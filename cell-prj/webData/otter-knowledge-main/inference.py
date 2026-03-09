import json
import os
import random
from typing import List, Optional

import numpy as np
import torch

from embeddings.morgan_fingerprint import MorganFingerprint


def set_seed(seed: int = 42) -> None:
    np.random.seed(seed)
    random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    os.environ["PYTHONHASHSEED"] = str(seed)
    print(f"Random seed set to: {seed}")


set_seed(42)

_net = None
_relation_map = None
_device = None
_initial_model = None
_model_loaded = False

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MODEL_DIR = os.path.join(BASE_DIR, "models")


def _build_initial_model():
    try:
        from embeddings.morgan_fingerprint import MorganGenerator  # type: ignore

        return MorganGenerator()
    except ImportError:
        print("MorganGenerator is unavailable; falling back to MorganFingerprint.")
        return MorganFingerprint()


def init_model(model_path: str = DEFAULT_MODEL_DIR, no_cuda: bool = False) -> bool:
    global _net, _relation_map, _device, _initial_model, _model_loaded

    try:
        _device = torch.device("cuda" if torch.cuda.is_available() and not no_cuda else "cpu")
        _initial_model = _build_initial_model()

        model_file = os.path.join(model_path, "model.pt")
        relation_map_path = os.path.join(model_path, "relation_map.json")

        if not os.path.exists(model_file):
            raise FileNotFoundError(f"Model file not found: {model_file}")
        if not os.path.exists(relation_map_path):
            raise FileNotFoundError(f"Relation map file not found: {relation_map_path}")

        _net = torch.load(model_file, map_location=torch.device(_device))
        if hasattr(_net, "eval"):
            _net.eval()

        with open(relation_map_path, encoding="utf-8") as f:
            _relation_map = json.load(f)

        _model_loaded = True
        print(f"Model initialized successfully on device: {_device}")
        return True
    except Exception as exc:
        print(f"Model initialization failed: {exc}")
        _model_loaded = False
        return False


def get_embedding(
    smiles: str,
    model_path: str = DEFAULT_MODEL_DIR,
    no_cuda: bool = False,
) -> Optional[List[float]]:
    global _net, _relation_map, _device, _initial_model, _model_loaded

    try:
        if not _model_loaded and not init_model(model_path, no_cuda):
            return None

        initial_embeddings = _initial_model.get_embedding([smiles])
        rel_id = _relation_map["smiles"]
        final_embedding: List[float] = []

        for init_embedding in initial_embeddings:
            nodes = {
                "morgan-fingerprint": {
                    "embeddings": init_embedding.unsqueeze(0).to(_device),
                    "node_indices": torch.tensor([0]).to(_device),
                },
                "Drug": {
                    "embeddings": [None],
                    "node_indices": torch.tensor([1]).to(_device),
                },
            }
            triples = torch.tensor([[0], [rel_id], [1]]).to(_device)

            with torch.no_grad():
                node_output_embeddings = _net.encoder(nodes, triples)
                output_embedding = node_output_embeddings[1]

            output_np = np.asarray(output_embedding.detach().cpu().numpy()).reshape(-1)
            init_np = np.asarray(init_embedding.detach().cpu().numpy()).reshape(-1)
            combined = np.concatenate([output_np, init_np])
            if combined.shape[0] < 128:
                combined = np.pad(combined, (0, 128 - combined.shape[0]))
            final_embedding.extend(combined[:128].tolist())

        return final_embedding or None
    except Exception as exc:
        print(f"Drug embedding failed: {exc}")
        return None
