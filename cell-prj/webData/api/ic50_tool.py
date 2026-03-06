# webData/api/ic50_tool.py
from __future__ import annotations

import sys
import tempfile
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import torch

# --- Paths ---
PROJECT_ROOT = Path(__file__).resolve().parents[3]
TRANSC_DIR = PROJECT_ROOT / "TransCDR" / "TransCDR-master" / "modify"
ALLCODES_DIR = TRANSC_DIR / "allCodes"
RESOURCE_DIR = TRANSC_DIR / "myPretrainedModel"
PRETRAINED_MODEL_PATH = PROJECT_ROOT / "TransCDR" / "TransCDR-master" / "result" / "Final_model" / "regression" / "model.pt"
DRUG_FILE = RESOURCE_DIR / "DB_alldrugs_smiles_129vec.csv"
GENE_FILE = RESOURCE_DIR / "geneid2vec128_978.csv"
TARGET_DSNAE_PATH = RESOURCE_DIR / "finetuned" / "target_dsnae.pt"
OTTER_DIR = PROJECT_ROOT / "cell-prj" / "webData" / "otter-knowledge-main"

if str(ALLCODES_DIR) not in sys.path:
    sys.path.append(str(ALLCODES_DIR))

from extract_fc_layers import load_transcdr_fc_weights  # noqa: E402
from utilis.data_utilis import read_drug_string, read_gene  # noqa: E402
from AE_part.dsn_AE import dsn_AE  # noqa: E402
from AE_part.Decoder import Decoder  # noqa: E402
from CE_part.CE_model import CE  # noqa: E402


def _ensure_numpy_core_aliases() -> None:
    if "numpy._core" in sys.modules:
        return
    import types
    import importlib

    numpy_core_module = types.ModuleType("numpy._core")
    sys.modules["numpy._core"] = numpy_core_module
    for submodule in ("multiarray", "overrides", "umath", "_multiarray_umath"):
        module_obj = importlib.import_module(f"numpy.core.{submodule}")
        sys.modules[f"numpy._core.{submodule}"] = module_obj
        setattr(numpy_core_module, submodule.split(".")[-1], module_obj)


class Datareader4p:
    def __init__(self, drug_file, gene_file, data_path, device, result_path, batch_size=32):
        self.device = device
        self.drug, self.drug_vec = read_drug_string(drug_file)
        self.gene = read_gene(gene_file, device)
        drug, cell, time_dose = self.read_data(data_path)
        feature_dict = self.trans_to_tensor(drug, cell, time_dose, self.drug_vec, result_path)
        self.ft = np.concatenate(
            (feature_dict["pert_time"], feature_dict["cell_id"], feature_dict["pert_idose"], feature_dict["drug"]),
            axis=1,
        )
        self.batch_num = max(1, len(self.ft) // batch_size)

    def read_data(self, adata_path):
        data = pd.read_csv(adata_path)
        drug = data.loc[:, "pert_id"]
        cell = data.iloc[:, 4:]
        time_dose = data.loc[:, ["pert_itime", "pert_idose"]]
        return np.asarray(drug), np.asarray(cell, dtype=np.float64), np.asarray(time_dose, dtype=np.float64)

    def trans_to_tensor(self, ft_drug, ft_cell, ft_time_dose, drug_vec, result_path):
        drug_feature = []
        for ft in ft_drug:
            drug_fp = drug_vec[ft]
            drug_feature.append(drug_fp)

        cell_id_feature_list = [np.array(ft, dtype=np.float64) for ft in ft_cell]

        feature_dict = dict()
        feature_dict["drug"] = np.asarray(drug_feature)

        ft_time = ft_time_dose[:, 0]
        ft_dose = ft_time_dose[:, 1]

        time_scaler = _load_scaler(result_path, "time_scaler.pkl")
        dose_scaler = _load_scaler(result_path, "dose_scaler.pkl")
        ft_time = time_scaler.transform(ft_time.reshape(-1, 1))
        ft_dose = np.log(ft_dose + 1e-6)
        ft_dose = dose_scaler.transform(ft_dose.reshape(-1, 1))

        feature_dict["pert_time"] = ft_time
        feature_dict["cell_id"] = np.asarray(cell_id_feature_list, dtype=np.float64)
        feature_dict["pert_idose"] = ft_dose

        return feature_dict

    def get_batch_data(self, batch_size):
        feature = torch.tensor(self.ft, dtype=torch.float64, device=self.device)
        for start_idx in range(0, feature.shape[0], batch_size):
            excerpt = slice(start_idx, start_idx + batch_size)
            output = dict()
            output["drug"] = feature[excerpt, 980:].clone().detach()
            output["pert_time"] = feature[excerpt, 0].clone().detach()
            output["cell_id"] = feature[excerpt, 1:979].clone().detach()
            output["pert_idose"] = feature[excerpt, 979].clone().detach()
            yield output


def _load_scaler(result_path: Path, name: str):
    import joblib

    scaler_path = Path(result_path) / name
    if not scaler_path.exists():
        raise FileNotFoundError(f"Scaler not found: {scaler_path}")
    return joblib.load(str(scaler_path))


class CodeToFCAdapter(torch.nn.Module):
    def __init__(self, num_gene, code_dim, target_dim, device):
        super().__init__()
        self.num_gene = num_gene
        self.code_dim = code_dim
        self.target_dim = target_dim
        self.slot_dim = int(np.ceil(target_dim / code_dim))
        self.project = torch.nn.Linear(
            num_gene, self.slot_dim, bias=False, dtype=torch.float64, device=device
        )
        self.post_linear = None
        if self.slot_dim * code_dim != target_dim:
            self.post_linear = torch.nn.Linear(
                self.slot_dim * code_dim, target_dim, bias=False, dtype=torch.float64, device=device
            )

    def forward(self, code):
        if code.dim() != 3:
            raise ValueError(f"Code tensor must be 3D (batch, gene, dim), got shape {code.shape}")
        x = code.transpose(1, 2)
        x = self.project(x)
        x = x.transpose(1, 2).contiguous().view(code.shape[0], -1)
        if self.post_linear is not None:
            x = self.post_linear(x)
        return x


class DSNAEToTransCDRPredictor(torch.nn.Module):
    def __init__(self, fc_layers, num_gene, code_dim, fc_input_dim, device):
        super().__init__()
        self.fc_layers = fc_layers.to(device).double()
        self.fc_layers.eval()
        for param in self.fc_layers.parameters():
            param.requires_grad = False
        self.adapter = CodeToFCAdapter(
            num_gene=num_gene, code_dim=code_dim, target_dim=fc_input_dim, device=device
        )
        self.fc_device = next(self.fc_layers.parameters()).device
        self.fc_dtype = next(self.fc_layers.parameters()).dtype

    def forward(self, code):
        adapted = self.adapter(code)
        adapted = adapted.to(device=self.fc_device, dtype=self.fc_dtype)
        return self.fc_layers(adapted)


class TransCDRPredictor:
    def __init__(self) -> None:
        _ensure_numpy_core_aliases()
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        fc_model, _ = load_transcdr_fc_weights(
            pretrained_model_path=str(PRETRAINED_MODEL_PATH),
            input_dim=1536,
        )
        self.fc_input_dim = fc_model.predictor[0].in_features

        params = {
            "drug_input_dim": 128,
            "gene_input_dim": 128,
            "drug_gene_emb_dim": 128,
            "CE_latent_dim": 256,
            "CE_output_dim": 128,
            "CE_drop": 0.2,
            "CE_n_layers": 1,
            "CE_n_heads": 2,
            "initializer": torch.nn.init.kaiming_uniform_,
            "device": self.device,
            "hid_dim": 1024,
            "gene_num": 978,
            "cell_id_input_dim": 978,
            "pert_time_emb_dim": 4,
            "cell_id_emb_dim": 128,
            "pert_dose_emb_dim": 4,
            "pt_lr": 0.00037,
            "pt_batch": 32,
            "pt_pert_epochs": 100,
            "pt_ctrl_epochs": 100,
            "ende_drop": 0.7,
            "ft_batch": 16,
            "hid_dim_list": [512, 512, 256, 256, 128],
            "latent_dim": 256,
            "ft_lr": 7e-05,
            "ft_epochs": 300,
            "gp": 10,
            "alpha": 1.0,
            "belta": 1.0,
            "conf_drop": 0.5,
            "ae_lr": 1.6e-05,
            "cc_lr": 0.0024,
            "tae_lr": 0.00017,
            "cc_output_dim": 128,
            "pred_drop": 0.5,
            "ende_epochs": 500,
            "ende_batch": 16,
            "critic_epochs": 500,
            "critic_batch": 16,
        }

        CE_model = CE(
            drug_input_dim=params["drug_input_dim"],
            gene_input_dim=params["gene_input_dim"],
            drug_gene_embed_dim=params["drug_gene_emb_dim"],
            device=self.device,
            hid_dim=params["hid_dim"],
            num_gene=params["gene_num"],
            CElatent_dim=params["CE_latent_dim"],
            drop=params["CE_drop"],
            CE_output_dim=params["CE_output_dim"],
            n_layers=params["CE_n_layers"],
            n_heads=params["CE_n_heads"],
            cell_id_input_dim=params["cell_id_input_dim"],
            cell_id_emb_dim=params["cell_id_emb_dim"],
            initializer=params["initializer"],
        )

        share_encoder = CE(
            drug_input_dim=params["drug_input_dim"],
            gene_input_dim=params["gene_input_dim"],
            drug_gene_embed_dim=params["drug_gene_emb_dim"],
            device=self.device,
            hid_dim=params["hid_dim"],
            num_gene=params["gene_num"],
            CElatent_dim=params["CE_latent_dim"],
            drop=params["CE_drop"],
            CE_output_dim=params["CE_output_dim"],
            n_layers=params["CE_n_layers"],
            n_heads=params["CE_n_heads"],
            cell_id_input_dim=params["cell_id_input_dim"],
            cell_id_emb_dim=params["cell_id_emb_dim"],
            initializer=params["initializer"],
        )

        share_decoder = Decoder(
            input_dim=params["CE_output_dim"] * 2,
            hid_dim_list=params["hid_dim_list"],
            latent_dim=share_encoder.linear_dim,
            drop=params["ende_drop"],
            device=self.device,
        )

        self.target_dsnae = dsn_AE(
            private_encoder=CE_model,
            share_encoder=share_encoder,
            decoder=share_decoder,
            **params,
        )
        self.target_dsnae.load_state_dict(torch.load(TARGET_DSNAE_PATH, map_location=self.device))
        self.target_dsnae.eval()

        code_dim = params["CE_output_dim"] * 2
        self.fc_predictor = DSNAEToTransCDRPredictor(
            fc_layers=fc_model,
            num_gene=params["gene_num"],
            code_dim=code_dim,
            fc_input_dim=self.fc_input_dim,
            device=self.device,
        )
        self.fc_predictor.eval()

    def predict(self, data_path: Path) -> float:
        data = Datareader4p(
            drug_file=str(DRUG_FILE),
            gene_file=str(GENE_FILE),
            data_path=str(data_path),
            device=self.device,
            result_path=str(RESOURCE_DIR),
            batch_size=1,
        )
        with torch.no_grad():
            for batch in data.get_batch_data(batch_size=1):
                x = batch["drug"]
                pert_time = batch["pert_time"]
                cell_id = batch["cell_id"]
                pert_idose = batch["pert_idose"]
                code, _ = self.target_dsnae.encode(x, data.gene, pert_time, cell_id, pert_idose)
                pred = self.fc_predictor(code).squeeze(-1)
                return float(pred.detach().cpu().double().numpy()[0])
        raise RuntimeError("No prediction output.")


class IC50Tool:
    def __init__(self) -> None:
        self._predictor = TransCDRPredictor()
        self._gene_ids = pd.read_csv(GENE_FILE, usecols=[0], header=None).iloc[:, 0].tolist()
        self._gene_lookup = {}
        for idx, gene_id in enumerate(self._gene_ids):
            key = str(int(gene_id))
            self._gene_lookup[key] = idx
            self._gene_lookup[f"GENE_{key}"] = idx

        self._smiles_to_id = {}
        self._smiles_to_vec_idx: Dict[str, int] = {}
        vectors: List[np.ndarray] = []
        with open(DRUG_FILE, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split(",")
                if len(parts) < 130:
                    continue
                drug_id, smiles = parts[0], parts[1]
                vec_strs = parts[2:130]
                if len(vec_strs) != 128:
                    continue
                vec = np.array(vec_strs, dtype=np.float32)
                vectors.append(vec)
                if smiles not in self._smiles_to_id:
                    self._smiles_to_id[smiles] = drug_id
                self._smiles_to_vec_idx[smiles] = len(vectors) - 1

        if vectors:
            self._drug_vectors = np.vstack(vectors)
            norms = np.linalg.norm(self._drug_vectors, axis=1, keepdims=True)
            norms[norms == 0] = 1.0
            self._drug_vectors_norm = self._drug_vectors / norms
        else:
            self._drug_vectors = np.zeros((0, 128), dtype=np.float32)
            self._drug_vectors_norm = self._drug_vectors

    def _compute_otter_embedding(self, smiles: str) -> List[float]:
        if str(OTTER_DIR) not in sys.path:
            sys.path.append(str(OTTER_DIR))
        from inference import get_embedding

        embedding = get_embedding(smiles)
        if embedding is None:
            raise ValueError("otter-knowledge failed to compute embedding.")
        if len(embedding) != 128:
            raise ValueError(f"Expected 128-d embedding, got {len(embedding)}.")
        return embedding

    def _append_drug_vector(self, drug_id: str, smiles: str, embedding: Sequence[float]) -> None:
        line = ",".join([drug_id, smiles] + [f"{float(x)}" for x in embedding]) + "\n"
        with open(DRUG_FILE, "a", encoding="utf-8") as f:
            f.write(line)

    def _update_vector_cache(self, smiles: str, embedding: Sequence[float]) -> None:
        vec = np.asarray(embedding, dtype=np.float32)
        norm = np.linalg.norm(vec)
        if norm == 0:
            norm = 1.0
        vec_norm = vec / norm
        if self._drug_vectors.size == 0:
            self._drug_vectors = vec[None, :]
            self._drug_vectors_norm = vec_norm[None, :]
            idx = 0
        else:
            self._drug_vectors = np.vstack([self._drug_vectors, vec])
            self._drug_vectors_norm = np.vstack([self._drug_vectors_norm, vec_norm])
            idx = self._drug_vectors.shape[0] - 1
        self._smiles_to_vec_idx[smiles] = idx

    def _similarity_confidence_from_embedding(self, embedding: Sequence[float]) -> Optional[float]:
        if self._drug_vectors_norm.size == 0:
            return None
        vec = np.asarray(embedding, dtype=np.float32)
        norm = np.linalg.norm(vec)
        if norm == 0:
            return None
        vec_norm = vec / norm
        sims = self._drug_vectors_norm @ vec_norm
        max_sim = float(np.max(sims))
        return (max_sim + 1.0) / 2.0

    def _resolve_drug_id(self, context: Dict[str, Any]) -> Tuple[str, Optional[List[float]]]:
        drug_id = context.get("drugId")
        if drug_id:
            return str(drug_id), None
        smiles = context.get("smiles")
        if not smiles:
            raise ValueError("Missing smiles; cannot map to drug ID.")
        match = self._smiles_to_id.get(smiles)
        if match is None:
            embedding = self._compute_otter_embedding(smiles)
            new_id = f"CUSTOM_{uuid.uuid4().hex[:10]}"
            self._append_drug_vector(new_id, smiles, embedding)
            self._smiles_to_id[smiles] = new_id
            return new_id, embedding
        return match, None

    def predict(self, expression: np.ndarray, context: Dict[str, Any]) -> Dict[str, Any]:
        if expression.shape[0] != len(self._gene_ids):
            raise ValueError(f"Expression vector length must be {len(self._gene_ids)}.")

        metadata = context.get("metadata") or {}
        time_val = metadata.get("timeHours")
        dose_val = metadata.get("doseUm")
        if time_val is None or dose_val is None:
            raise ValueError("Missing timeHours or doseUm in metadata.")

        smiles = context.get("smiles")
        drug_id, embedding = self._resolve_drug_id(context)
        cell_line = context.get("cellLineId") or "NA"

        row = {
            "pert_itime": float(time_val),
            "pert_id": drug_id,
            "cell_iname": str(cell_line),
            "pert_idose": float(dose_val),
        }
        for gene_id, value in zip(self._gene_ids, expression):
            row[str(int(gene_id))] = float(value)

        df = pd.DataFrame([row])
        tmp_dir = Path(tempfile.gettempdir()) / "transcdr_inputs"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        tmp_path = tmp_dir / f"transcdr_{uuid.uuid4().hex}.csv"
        try:
            df.to_csv(tmp_path, index=False)
            ic50_value = self._predictor.predict(tmp_path)
        finally:
            try:
                tmp_path.unlink(missing_ok=True)
            except Exception:
                pass

        confidence = None
        if smiles:
            if embedding is not None:
                confidence = self._similarity_confidence_from_embedding(embedding)
                if smiles not in self._smiles_to_vec_idx:
                    self._update_vector_cache(smiles, embedding)
            else:
                idx = self._smiles_to_vec_idx.get(smiles)
                if idx is not None and self._drug_vectors.size > 0:
                    confidence = self._similarity_confidence_from_embedding(self._drug_vectors[idx])
                else:
                    try:
                        confidence = self._similarity_confidence_from_embedding(self._compute_otter_embedding(smiles))
                    except Exception:
                        confidence = None

        ln_ic50 = float(ic50_value)
        ic50 = float(np.exp(ln_ic50))
        return {
            "lnIc50": ln_ic50,
            "ic50": ic50,
            "unit": metadata.get("ic50Unit", "uM"),
            "confidence": confidence,
            "toolVersion": "TransCDR",
            "drugId": drug_id,
        }


_predictor: IC50Tool | None = None


def expression_to_vector(expression: Sequence[Dict[str, float]]) -> np.ndarray:
    gene_ids = pd.read_csv(GENE_FILE, usecols=[0], header=None).iloc[:, 0].tolist()
    gene_lookup = {}
    for idx, gene_id in enumerate(gene_ids):
        key = str(int(gene_id))
        gene_lookup[key] = idx
        gene_lookup[f"GENE_{key}"] = idx

    vector = np.zeros(len(gene_ids), dtype=np.float64)
    for item in expression:
        key = item.get("geneId")
        if key is None:
            continue
        idx = gene_lookup.get(key)
        if idx is None:
            continue
        vector[idx] = float(item.get("expressionLevel", 0.0))
    return vector


def get_ic50_predictor() -> IC50Tool:
    global _predictor
    if _predictor is None:
        _predictor = IC50Tool()
    return _predictor
