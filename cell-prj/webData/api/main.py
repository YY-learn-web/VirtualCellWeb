from __future__ import annotations

import base64
import io
import json
import logging
import math
import os
import shutil
import sys
import traceback
import uuid
import urllib.parse
import urllib.request
from pathlib import Path
from threading import Lock
from time import perf_counter
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
from fastapi import BackgroundTasks, FastAPI, File, Form, HTTPException, Query, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

try:
    from scipy.stats import hypergeom
except Exception:
    hypergeom = None

LOGGER = logging.getLogger("virtual-cell-api")
if not LOGGER.handlers:
    logging.basicConfig(level=logging.INFO)

API_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = API_DIR.parents[2]
WEBDATA_DIR = API_DIR.parent

BATCH_UPLOAD_DIR = WEBDATA_DIR / "batch" / "uploads"
BATCH_RESULTS_DIR = WEBDATA_DIR / "batch" / "results"
SWEEP_EXPORT_DIR = WEBDATA_DIR / "batch" / "sweep"
IC50_EXPORT_DIR = WEBDATA_DIR / "batch" / "ic50"

METAB_BACKEND_DIR = PROJECT_ROOT / "MetabolomicsPred" / "backend"
METAB_RESULTS_DIR = API_DIR / "metabolomics_results"
METAB_UPLOAD_DIR = API_DIR / "metabolomics_uploads"
METAB_ASSETS_DIR = METAB_BACKEND_DIR / "app" / "assets"

ENRICH_LOCAL_DIR = API_DIR / "assets" / "enrichment"
ENRICH_LIBRARY_FILES = {
    "GO_Biological_Process_2021": "GO_Biological_Process_2021.gmt",
    "GO_Molecular_Function_2021": "GO_Molecular_Function_2021.gmt",
    "GO_Cellular_Component_2021": "GO_Cellular_Component_2021.gmt",
}

for folder in (
    BATCH_UPLOAD_DIR,
    BATCH_RESULTS_DIR,
    SWEEP_EXPORT_DIR,
    IC50_EXPORT_DIR,
    METAB_RESULTS_DIR,
    METAB_UPLOAD_DIR,
):
    folder.mkdir(parents=True, exist_ok=True)

if str(API_DIR) not in sys.path:
    sys.path.insert(0, str(API_DIR))
if str(METAB_BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(METAB_BACKEND_DIR))

ONNX_IMPORT_ERROR: Optional[str] = None
IC50_IMPORT_ERROR: Optional[str] = None
METAB_IMPORT_ERROR: Optional[str] = None

try:
    from onnx_model import load_onnx_model
except Exception as exc:
    load_onnx_model = None
    ONNX_IMPORT_ERROR = str(exc)

import pandas as pd

try:
    from ic50_tool import expression_to_vector, get_ic50_predictor
except Exception as exc:
    expression_to_vector = None
    get_ic50_predictor = None
    IC50_IMPORT_ERROR = str(exc)

try:
    from app.algorithms.scfea_analysis import generate_custom_diff_plot, run_downstream_analysis
    from app.algorithms.scfea_core import run_scfea_training
except Exception as exc:
    generate_custom_diff_plot = None
    run_downstream_analysis = None
    run_scfea_training = None
    METAB_IMPORT_ERROR = str(exc)

app = FastAPI(title="Virtual Cell Predictor API")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
app.mount("/api/metabolomics/results", StaticFiles(directory=str(METAB_RESULTS_DIR)), name="metabolomics_results")

_onnx_model = None
_onnx_model_error: Optional[str] = None
_onnx_model_lock = Lock()

_batch_result_index: Dict[str, Path] = {}
_sweep_cache: Dict[str, Dict[str, Any]] = {}
_download_index: Dict[str, Path] = {}

_metab_tasks: Dict[str, Dict[str, Any]] = {}
_metab_task_lock = Lock()

_entrez_to_symbol_cache: Optional[Dict[str, str]] = None
_gmt_cache: Dict[str, Dict[str, set[str]]] = {}


class RangeConfig(BaseModel):
    start: float
    end: float
    steps: int


class BatchSweepRequest(BaseModel):
    sweepVariable: str
    range: RangeConfig
    fixedParamValue: float
    smiles: str
    cellLineId: str = "PC3"


class ExpressionItem(BaseModel):
    geneId: str
    expressionLevel: float


class EnrichmentRequest(BaseModel):
    geneIds: List[str]
    library: str
    topN: int = 20


class StringNetworkRequest(BaseModel):
    genes: List[ExpressionItem]
    minExpression: float = 1.5
    species: int = 9606
    requiredScore: int = 400
    networkType: str = "functional"
    maxGenes: int = 200


class SingleMetabolomicsRequest(BaseModel):
    expression: List[ExpressionItem]
    epochs: int = 100
    imputation: bool = False


class KeyGeneRequest(BaseModel):
    expression: List[ExpressionItem]
    enrichment: List[Dict[str, Any]] = Field(default_factory=list)
    ppiGenes: List[str] = Field(default_factory=list)
    metabolomicsFlux: Dict[str, float] = Field(default_factory=dict)
    topN: int = 5


class BatchSnapshotRequest(BaseModel):
    smiles: str
    cellLineId: str
    sweepVariable: str
    fixedParamValue: float
    range: RangeConfig
    snapshotValue: float


def _safe_float(value: Any, default: float) -> float:
    try:
        return float(value)
    except Exception:
        return default


def _sweep_key(payload: Dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, ensure_ascii=True)


def _register_download(path: Path, filename: Optional[str] = None) -> str:
    name = filename or path.name
    _download_index[name] = path
    return f"/api/download/{name}"


def _ensure_onnx_model():
    global _onnx_model, _onnx_model_error
    if _onnx_model is not None:
        return _onnx_model
    with _onnx_model_lock:
        if _onnx_model is not None:
            return _onnx_model
        if load_onnx_model is None:
            detail = f"ASCEND model import failed: {ONNX_IMPORT_ERROR}"
            raise HTTPException(status_code=503, detail=detail)
        try:
            _onnx_model = load_onnx_model()
        except Exception as exc:
            _onnx_model_error = str(exc)
            raise HTTPException(status_code=503, detail=f"ASCEND model unavailable: {_onnx_model_error}") from exc
        return _onnx_model


def _predict_single(smiles: str, time_value: float, dose_value: float, cell_line_id: str) -> List[Tuple[str, float]]:
    model = _ensure_onnx_model()
    if hasattr(model, "predict_single"):
        return model.predict_single(smiles, float(time_value), float(dose_value), cell_line_id)
    if hasattr(model, "predict"):
        return model.predict(smiles, float(time_value), float(dose_value), cell_line_id)
    raise HTTPException(status_code=500, detail="ASCEND model wrapper does not expose predict function.")


def _run_sweep(payload: Dict[str, Any]) -> Dict[str, Any]:
    request = {
        "sweepVariable": str(payload.get("SweepVariable") or payload.get("sweepVariable") or "dose").lower(),
        "range": payload.get("range") or {},
        "fixedParamValue": payload.get("fixedParamValue", 24),
        "smiles": payload.get("smiles"),
        "cellLineId": payload.get("cellLineId") or "PC3",
    }
    start = _safe_float(request["range"].get("start"), 0.0)
    end = _safe_float(request["range"].get("end"), 100.0)
    steps = int(request["range"].get("steps", 6))
    if steps < 2:
        raise HTTPException(status_code=400, detail="range.steps must be >= 2.")
    if not request["smiles"]:
        raise HTTPException(status_code=400, detail="smiles is required.")
    if request["sweepVariable"] not in {"time", "dose"}:
        raise HTTPException(status_code=400, detail="sweepVariable must be 'time' or 'dose'.")

    fixed_value = request["fixedParamValue"]
    if isinstance(fixed_value, dict):
        if request["sweepVariable"] == "dose":
            fixed_value = fixed_value.get("time", 24)
        else:
            fixed_value = fixed_value.get("dose", 10)

    x_values = np.linspace(start, end, steps)
    expression_rows: List[Dict[str, float]] = []
    point_meta: List[Dict[str, float]] = []

    for x in x_values:
        if request["sweepVariable"] == "dose":
            time_value = float(fixed_value)
            dose_value = float(x)
        else:
            time_value = float(x)
            dose_value = float(fixed_value)
        predicted = _predict_single(request["smiles"], time_value, dose_value, request["cellLineId"])
        gene_map = {str(gene_id): float(value) for gene_id, value in predicted}
        expression_rows.append(gene_map)
        point_meta.append(
            {
                "xValue": float(x),
                "time": float(time_value),
                "dose": float(dose_value),
            }
        )

    matrix = pd.DataFrame(expression_rows, index=[f"point_{i:03d}" for i in range(len(expression_rows))]).fillna(0.0)
    if matrix.empty:
        raise HTTPException(status_code=500, detail="Sweep prediction returned empty result.")
    variability = matrix.var(axis=0).sort_values(ascending=False)
    top_genes = [str(g) for g in variability.head(5).index.tolist()]

    points = []
    for i, meta in enumerate(point_meta):
        values = {g: float(matrix.iloc[i][g]) for g in top_genes}
        points.append({"xValue": round(meta["xValue"], 6), "genes": values})

    return {
        "request": request,
        "points": points,
        "topGenes": top_genes,
        "matrix": matrix,
        "meta": point_meta,
    }


def _coerce_numeric_df(df: pd.DataFrame) -> pd.DataFrame:
    coerced = df.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
    return coerced.fillna(0.0)


def _extract_gene_columns(df: pd.DataFrame) -> List[str]:
    return [c for c in df.columns if str(c).startswith("GENE_")]


def _batch_result_file(file_id: str) -> Path:
    if file_id in _batch_result_index:
        return _batch_result_index[file_id]
    candidate = BATCH_RESULTS_DIR / f"result_{file_id}.csv"
    if candidate.exists():
        _batch_result_index[file_id] = candidate
        return candidate
    raise HTTPException(status_code=404, detail=f"Batch result file not found for file_id={file_id}.")


def _summary_pairs_from_matrix(df: pd.DataFrame, top_n: int = 500) -> List[List[Any]]:
    gene_cols = _extract_gene_columns(df)
    if not gene_cols:
        return []
    matrix = _coerce_numeric_df(df[gene_cols])
    summary = matrix.abs().mean(axis=0).sort_values(ascending=False).head(top_n)
    return [[str(gene), float(score)] for gene, score in summary.items()]


def _load_entrez_symbol_map() -> Dict[str, str]:
    global _entrez_to_symbol_cache
    if _entrez_to_symbol_cache is not None:
        return _entrez_to_symbol_cache

    mapping: Dict[str, str] = {}
    for name in ("geneinfo_beta.txt", "gene_info.txt"):
        path = METAB_ASSETS_DIR / name
        if not path.exists():
            continue
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            continue
        if df.empty:
            continue
        cols = {c.lower(): c for c in df.columns}
        gene_col = cols.get("gene_id") or cols.get("entrez_id") or cols.get("geneid")
        symbol_col = cols.get("gene_symbol") or cols.get("symbol")
        if not gene_col or not symbol_col:
            continue
        sub = df[[gene_col, symbol_col]].dropna()
        for _, row in sub.iterrows():
            entrez = str(row[gene_col]).strip()
            symbol = str(row[symbol_col]).strip()
            if entrez and symbol:
                mapping[entrez] = symbol
        if mapping:
            break
    _entrez_to_symbol_cache = mapping
    return mapping


def _gene_to_symbol(gene_id: str) -> Optional[str]:
    gene_id = str(gene_id).strip()
    if not gene_id:
        return None
    mapping = _load_entrez_symbol_map()
    if gene_id.startswith("GENE_"):
        entrez = gene_id[5:]
        return mapping.get(entrez)
    if gene_id.isdigit():
        return mapping.get(gene_id)
    return gene_id


def _normalize_gene_ids(gene_ids: Iterable[str]) -> Tuple[List[str], List[str]]:
    mapped: List[str] = []
    unmapped: List[str] = []
    for gene_id in gene_ids:
        symbol = _gene_to_symbol(gene_id)
        if symbol is None:
            unmapped.append(str(gene_id))
        else:
            mapped.append(symbol)
    uniq = list(dict.fromkeys(mapped))
    return uniq, unmapped


def _library_gmt_path(library: str) -> Path:
    filename = ENRICH_LIBRARY_FILES.get(library)
    if filename is None:
        raise HTTPException(status_code=400, detail=f"Unsupported enrichment library: {library}")
    path = ENRICH_LOCAL_DIR / filename
    if not path.exists():
        raise HTTPException(
            status_code=503,
            detail=f"Offline GO library file is missing: {path}. Please place the GMT file and retry.",
        )
    return path


def _load_gmt(path: Path) -> Dict[str, set[str]]:
    cache_key = str(path.resolve())
    if cache_key in _gmt_cache:
        return _gmt_cache[cache_key]
    terms: Dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = {g for g in parts[2:] if g}
            if genes:
                terms[term] = genes
    _gmt_cache[cache_key] = terms
    return terms


def _bh_adjust(p_values: List[float]) -> List[float]:
    if not p_values:
        return []
    arr = np.asarray(p_values, dtype=float)
    order = np.argsort(arr)
    ranked = arr[order]
    n = float(len(arr))
    adjusted = np.empty_like(ranked)
    running = 1.0
    for i in range(len(ranked) - 1, -1, -1):
        raw = ranked[i] * n / float(i + 1)
        running = min(running, raw)
        adjusted[i] = running
    out = np.empty_like(adjusted)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    return out.tolist()


def _offline_go_enrichment(gene_ids: List[str], library: str, top_n: int) -> Dict[str, Any]:
    symbols, unmapped = _normalize_gene_ids(gene_ids)
    if not symbols:
        raise HTTPException(status_code=400, detail="No valid genes after symbol normalization.")

    terms = _load_gmt(_library_gmt_path(library))
    if not terms:
        raise HTTPException(status_code=500, detail="Offline GO library is empty.")

    background = set().union(*terms.values())
    query = set([g for g in symbols if g in background])
    if not query:
        raise HTTPException(
            status_code=400,
            detail="No overlap between input genes and offline GO background. Check gene symbols.",
        )

    total_bg = len(background)
    total_query = len(query)
    rows = []
    p_values: List[float] = []
    for term, term_genes in terms.items():
        overlap = query.intersection(term_genes)
        if not overlap:
            continue
        k = len(overlap)
        n = len(term_genes)
        if hypergeom is not None:
            p_value = float(hypergeom.sf(k - 1, total_bg, n, total_query))
        else:
            # Conservative fallback when scipy is unavailable.
            p_value = float(min(1.0, max(1e-300, (n / max(total_bg, 1)) ** k)))
        p_values.append(p_value)
        rows.append(
            {
                "term": term,
                "pValue": p_value,
                "combinedScore": float(-math.log10(p_value + 1e-300) * k),
                "genes": sorted(overlap),
            }
        )

    if not rows:
        return {"success": True, "library": library, "data": [], "unmapped": unmapped}

    adjusted = _bh_adjust(p_values)
    for row, adj in zip(rows, adjusted):
        row["adjP"] = float(adj)
    rows.sort(key=lambda x: (x["adjP"], x["pValue"]))
    return {
        "success": True,
        "library": library,
        "data": rows[: max(1, int(top_n))],
        "unmapped": unmapped,
    }


def _fallback_ppi_image(genes: List[str]) -> str:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        payload = "PPI graph unavailable"
        b64 = base64.b64encode(payload.encode("utf-8")).decode("utf-8")
        return f"data:text/plain;base64,{b64}"

    n = max(1, len(genes))
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x = np.cos(theta)
    y = np.sin(theta)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_facecolor("#f8fafc")
    for i in range(n):
        j = (i + 1) % n
        ax.plot([x[i], x[j]], [y[i], y[j]], color="#94a3b8", linewidth=1, alpha=0.8)
    ax.scatter(x, y, s=220, c="#06b6d4", edgecolors="#0f172a", linewidths=0.4)
    for i, gene in enumerate(genes):
        ax.text(x[i], y[i], gene, fontsize=8, ha="center", va="center", color="#0f172a")
    ax.axis("off")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    return f"data:image/png;base64,{base64.b64encode(buf.getvalue()).decode('utf-8')}"


def _string_link(genes: List[str], species: int, required_score: int, network_type: str) -> str:
    identifiers = "\r".join(genes)
    encoded = urllib.parse.quote(identifiers, safe="")
    return (
        "https://string-db.org/cgi/network?identifiers="
        f"{encoded}&species={species}&required_score={required_score}&network_type={network_type}"
    )


def _update_metab_task(task_id: str, payload: Dict[str, Any]) -> None:
    with _metab_task_lock:
        current = _metab_tasks.get(task_id, {})
        current.update(payload)
        _metab_tasks[task_id] = current


def _ensure_metabolomics_available() -> None:
    if run_scfea_training is None or run_downstream_analysis is None:
        raise HTTPException(status_code=503, detail=f"Metabolomics service unavailable: {METAB_IMPORT_ERROR}")


def _process_metabolomics_task(task_id: str, job_id: str, input_file: Path, params: Dict[str, Any]) -> None:
    try:
        _update_metab_task(task_id, {"status": "TRAINING", "progress": 10, "message": "Training scFEA model..."})
        output_dir = METAB_RESULTS_DIR / job_id
        output_dir.mkdir(parents=True, exist_ok=True)

        flux_file, balance_file = run_scfea_training(
            input_file=str(input_file),
            output_dir=str(output_dir),
            assets_dir=str(METAB_ASSETS_DIR),
            sc_imputation=bool(params.get("imputation", False)),
            epochs=int(params.get("epochs", 100)),
        )

        _update_metab_task(task_id, {"status": "ANALYZING", "progress": 65, "message": "Running downstream analysis..."})
        analysis = run_downstream_analysis(
            flux_file=str(flux_file),
            balance_file=str(balance_file),
            output_dir=str(output_dir),
            assets_dir=str(METAB_ASSETS_DIR),
            n_clusters=int(params.get("n_clusters", 4)),
        )

        result = {
            "files": {
                "flux": str(flux_file),
                "balance": str(balance_file),
                **(analysis.get("files", {}) if isinstance(analysis, dict) else {}),
            },
            "images": analysis.get("images", {}) if isinstance(analysis, dict) else {},
        }
        _update_metab_task(
            task_id,
            {
                "status": "SUCCESS",
                "progress": 100,
                "message": "Metabolomics analysis completed.",
                "result": result,
            },
        )
    except Exception as exc:
        traceback.print_exc()
        _update_metab_task(
            task_id,
            {
                "status": "FAILURE",
                "progress": 0,
                "message": "Metabolomics analysis failed.",
                "error": str(exc),
            },
        )


def _expression_to_single_matrix(expression: List[ExpressionItem], sample_id: str = "sample_000") -> pd.DataFrame:
    row = {item.geneId: float(item.expressionLevel) for item in expression}
    df = pd.DataFrame([row], index=[sample_id])
    return _coerce_numeric_df(df)


def _key_gene_scores(payload: KeyGeneRequest) -> List[Dict[str, float]]:
    if not payload.expression:
        return []
    expr_abs = {item.geneId: abs(float(item.expressionLevel)) for item in payload.expression}
    expr_max = max(expr_abs.values()) if expr_abs else 1.0
    expr_norm = {gene: (value / expr_max if expr_max else 0.0) for gene, value in expr_abs.items()}

    enrich_weight_by_symbol: Dict[str, float] = {}
    for term in payload.enrichment:
        genes = term.get("genes") or []
        if not isinstance(genes, list):
            continue
        p = float(term.get("adjP") or term.get("pValue") or 1.0)
        weight = max(0.0, -math.log10(p + 1e-300))
        for gene in genes:
            key = str(gene)
            enrich_weight_by_symbol[key] = enrich_weight_by_symbol.get(key, 0.0) + weight
    enrich_max = max(enrich_weight_by_symbol.values()) if enrich_weight_by_symbol else 1.0

    ppi_symbols = set(_normalize_gene_ids(payload.ppiGenes)[0])
    flux_values = [abs(float(v)) for v in payload.metabolomicsFlux.values()]
    flux_strength = float(np.tanh(np.mean(flux_values))) if flux_values else 0.0

    scores: List[Dict[str, float]] = []
    for item in payload.expression:
        gene = item.geneId
        symbol = _gene_to_symbol(gene) or gene
        expression_score = float(expr_norm.get(gene, 0.0))
        enrichment_score = float(enrich_weight_by_symbol.get(symbol, 0.0) / enrich_max if enrich_max else 0.0)
        ppi_score = 1.0 if symbol in ppi_symbols else 0.0
        metabolomics_score = float(expression_score * flux_strength)
        final_score = 0.45 * expression_score + 0.25 * enrichment_score + 0.20 * ppi_score + 0.10 * metabolomics_score
        scores.append(
            {
                "gene": gene,
                "score": float(final_score),
                "expressionScore": expression_score,
                "enrichmentScore": enrichment_score,
                "ppiScore": ppi_score,
                "metabolomicsScore": metabolomics_score,
            }
        )
    scores.sort(key=lambda x: x["score"], reverse=True)
    top_n = max(1, int(payload.topN or 5))
    return scores[:top_n]


@app.get("/health")
def health() -> Dict[str, str]:
    return {"status": "ok"}


@app.post("/api/Analysis")
async def api_analysis(request: Dict[str, Any]) -> Dict[str, Any]:
    mode = str(request.get("AnalysisMode", "single")).lower()
    if mode == "single":
        return await api_analysis_single(request)
    if mode == "batch":
        method = str(request.get("BatchMethod", "sweep")).lower()
        if method == "sweep":
            return await api_batch_sweep(request)
        raise HTTPException(status_code=400, detail="Unsupported BatchMethod. Use 'sweep' or upload endpoint.")
    raise HTTPException(status_code=400, detail="AnalysisMode must be 'single' or 'batch'.")


@app.post("/api/Analysis/Single")
async def api_analysis_single(request: Dict[str, Any]) -> Dict[str, Any]:
    smiles = request.get("smiles")
    if not smiles:
        raise HTTPException(status_code=400, detail="smiles is required.")
    time_value = _safe_float(request.get("time"), 24.0)
    dose_value = _safe_float(request.get("dose"), 10.0)
    cell_line_id = str(request.get("cellLineId") or "PC3")

    start = perf_counter()
    predicted = _predict_single(smiles=smiles, time_value=time_value, dose_value=dose_value, cell_line_id=cell_line_id)
    elapsed_ms = int((perf_counter() - start) * 1000.0)

    return {
        "success": True,
        "mode": "single",
        "data": [[str(gene_id), float(value)] for gene_id, value in predicted],
        "metadata": {
            "inferenceTime": f"{elapsed_ms} ms",
            "modelVersion": "ASCEND-ONNX",
            "smiles": smiles,
            "time": time_value,
            "dose": dose_value,
            "cellLineId": cell_line_id,
        },
    }


@app.post("/api/Analysis/Batch/Sweep")
async def api_batch_sweep(request: Dict[str, Any]) -> Dict[str, Any]:
    cache_key = _sweep_key(
        {
            "sweepVariable": str(request.get("SweepVariable") or request.get("sweepVariable") or "dose").lower(),
            "range": request.get("range") or {},
            "fixedParamValue": request.get("fixedParamValue"),
            "smiles": request.get("smiles"),
            "cellLineId": request.get("cellLineId") or "PC3",
        }
    )
    result = _sweep_cache.get(cache_key)
    if result is None:
        result = _run_sweep(request)
        _sweep_cache[cache_key] = result

    return {
        "success": True,
        "type": "sweep",
        "data": result["points"],
        "topGenes": result["topGenes"],
    }


@app.post("/api/analysis/batch/upload")
async def api_batch_upload(file: UploadFile = File(...), cellLineId: Optional[str] = Form(default="PC3")) -> Dict[str, Any]:
    if not file.filename:
        raise HTTPException(status_code=400, detail="No file selected.")
    ext = Path(file.filename).suffix.lower()
    if ext not in {".csv", ".xlsx", ".xls"}:
        raise HTTPException(status_code=400, detail="Only CSV/XLS/XLSX files are supported.")

    file_id = uuid.uuid4().hex[:12]
    saved_input = BATCH_UPLOAD_DIR / f"input_{file_id}{ext}"
    with saved_input.open("wb") as f:
        shutil.copyfileobj(file.file, f)

    model = _ensure_onnx_model()
    if not hasattr(model, "predict_from_file"):
        raise HTTPException(status_code=500, detail="ASCEND model wrapper missing predict_from_file function.")

    try:
        result_df = model.predict_from_file(str(saved_input), file_id=file_id, cell_id=str(cellLineId or "PC3"))
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Batch upload prediction failed: {exc}") from exc

    output_file = BATCH_RESULTS_DIR / f"result_{file_id}.csv"
    if not output_file.exists():
        if isinstance(result_df, pd.DataFrame):
            output_file.parent.mkdir(parents=True, exist_ok=True)
            result_df.to_csv(output_file, index=False)
        else:
            raise HTTPException(status_code=500, detail="Batch prediction finished but result file is missing.")

    _batch_result_index[file_id] = output_file
    download_url = _register_download(output_file)
    processed_rows = int(result_df.shape[0]) if isinstance(result_df, pd.DataFrame) else 0

    return {
        "success": True,
        "type": "upload",
        "data": {
            "processedRows": processed_rows,
            "downloadUrl": download_url,
            "fileId": file_id,
        },
    }


@app.get("/api/download/{filename}")
async def api_download(filename: str) -> FileResponse:
    path = _download_index.get(filename)
    if path is None:
        candidates = [
            BATCH_RESULTS_DIR / filename,
            SWEEP_EXPORT_DIR / filename,
            IC50_EXPORT_DIR / filename,
        ]
        for candidate in candidates:
            if candidate.exists():
                path = candidate
                break
    if path is None or not path.exists():
        raise HTTPException(status_code=404, detail=f"File not found: {filename}")
    return FileResponse(str(path), filename=filename, media_type="application/octet-stream")


@app.get("/api/analysis/batch/export/{file_id}")
async def api_batch_export(file_id: str) -> FileResponse:
    file_path = _batch_result_file(file_id)
    return FileResponse(str(file_path), filename=file_path.name, media_type="text/csv")


@app.post("/api/analysis/batch/sweep/export")
async def api_sweep_export(request: Dict[str, Any]) -> FileResponse:
    cache_key = _sweep_key(
        {
            "sweepVariable": str(request.get("sweepVariable") or request.get("SweepVariable") or "dose").lower(),
            "range": request.get("range") or {},
            "fixedParamValue": request.get("fixedParamValue"),
            "smiles": request.get("smiles"),
            "cellLineId": request.get("cellLineId") or "PC3",
        }
    )
    result = _sweep_cache.get(cache_key)
    if result is None:
        result = _run_sweep(request)
        _sweep_cache[cache_key] = result

    export_name = f"sweep_{uuid.uuid4().hex[:12]}.csv"
    export_path = SWEEP_EXPORT_DIR / export_name
    matrix = result["matrix"].copy()
    matrix.index.name = "cell_id"
    matrix.to_csv(export_path)
    _register_download(export_path)
    return FileResponse(str(export_path), filename=export_name, media_type="text/csv")


def _ensure_ic50_predictor():
    if get_ic50_predictor is None or expression_to_vector is None:
        raise HTTPException(status_code=503, detail=f"IC50 service unavailable: {IC50_IMPORT_ERROR}")
    try:
        return get_ic50_predictor()
    except Exception as exc:
        raise HTTPException(status_code=503, detail=f"IC50 model failed to initialize: {exc}") from exc


@app.post("/api/ic50/predict")
async def api_ic50_predict(request: Dict[str, Any]) -> Dict[str, Any]:
    expression = request.get("expression") or []
    if not isinstance(expression, list) or len(expression) == 0:
        raise HTTPException(status_code=400, detail="expression is required.")

    predictor = _ensure_ic50_predictor()
    try:
        vector = expression_to_vector(expression)
        context = {
            "smiles": request.get("smiles"),
            "cellLineId": request.get("cellLineId"),
            "metadata": request.get("metadata") or {},
        }
        pred = predictor.predict(vector, context)
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"IC50 prediction failed: {exc}") from exc
    return {"success": True, "data": pred}


@app.post("/api/ic50/batch/predict")
async def api_ic50_batch_predict(request: Dict[str, Any]) -> Dict[str, Any]:
    file_id = str(request.get("fileId") or "").strip()
    if not file_id:
        raise HTTPException(status_code=400, detail="fileId is required.")
    preview_count = int(request.get("previewCount") or 10)
    predictor = _ensure_ic50_predictor()

    df = pd.read_csv(_batch_result_file(file_id))
    gene_cols = _extract_gene_columns(df)
    if not gene_cols:
        raise HTTPException(status_code=400, detail="No GENE_* columns found in uploaded batch file.")

    records: List[Dict[str, Any]] = []
    for _, row in df.iterrows():
        expression = [{"geneId": g, "expressionLevel": float(row[g])} for g in gene_cols]
        context = {
            "smiles": row.get("smiles"),
            "cellLineId": row.get("cellLineId") or row.get("cellLine") or request.get("cellLineId") or "PC3",
            "metadata": {
                "timeHours": float(row.get("time", 24)),
                "doseUm": float(row.get("dose", 10)),
            },
        }
        try:
            pred = predictor.predict(expression_to_vector(expression), context)
            records.append(
                {
                    "smiles": context["smiles"],
                    "time": context["metadata"]["timeHours"],
                    "dose": context["metadata"]["doseUm"],
                    "cellLineId": context["cellLineId"],
                    **pred,
                }
            )
        except Exception as exc:
            LOGGER.warning("IC50 batch row failed: %s", exc)

    out_df = pd.DataFrame(records)
    out_path = IC50_EXPORT_DIR / f"ic50_batch_{file_id}_{uuid.uuid4().hex[:8]}.csv"
    out_df.to_csv(out_path, index=False)
    download_url = _register_download(out_path)

    return {
        "success": True,
        "data": {
            "processedRows": int(len(records)),
            "preview": records[: max(1, preview_count)],
            "downloadUrl": download_url,
        },
    }


@app.post("/api/ic50/sweep/predict")
async def api_ic50_sweep_predict(request: Dict[str, Any]) -> Dict[str, Any]:
    preview_count = int(request.get("previewCount") or 10)
    predictor = _ensure_ic50_predictor()

    sweep_result = _run_sweep(
        {
            "sweepVariable": request.get("sweepVariable"),
            "range": request.get("range"),
            "fixedParamValue": request.get("fixedParamValue"),
            "smiles": request.get("smiles"),
            "cellLineId": request.get("cellLineId"),
        }
    )
    matrix = sweep_result["matrix"]
    meta = sweep_result["meta"]
    gene_cols = [str(c) for c in matrix.columns]

    records: List[Dict[str, Any]] = []
    for i in range(matrix.shape[0]):
        row = matrix.iloc[i]
        expression = [{"geneId": g, "expressionLevel": float(row[g])} for g in gene_cols]
        context = {
            "smiles": request.get("smiles"),
            "cellLineId": request.get("cellLineId") or "PC3",
            "metadata": {
                "timeHours": float(meta[i]["time"]),
                "doseUm": float(meta[i]["dose"]),
            },
        }
        try:
            pred = predictor.predict(expression_to_vector(expression), context)
            records.append(
                {
                    "smiles": context["smiles"],
                    "time": context["metadata"]["timeHours"],
                    "dose": context["metadata"]["doseUm"],
                    "cellLineId": context["cellLineId"],
                    **pred,
                }
            )
        except Exception as exc:
            LOGGER.warning("IC50 sweep point failed: %s", exc)

    out_df = pd.DataFrame(records)
    out_path = IC50_EXPORT_DIR / f"ic50_sweep_{uuid.uuid4().hex[:12]}.csv"
    out_df.to_csv(out_path, index=False)
    download_url = _register_download(out_path)

    return {
        "success": True,
        "data": {
            "processedRows": int(len(records)),
            "preview": records[: max(1, preview_count)],
            "downloadUrl": download_url,
        },
    }


@app.post("/api/metabolomics/analyze")
async def api_metabolomics_analyze(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(...),
    epochs: int = Form(100),
    imputation: bool = Form(False),
    n_clusters: int = Form(4),
) -> Dict[str, Any]:
    _ensure_metabolomics_available()

    if not file.filename:
        raise HTTPException(status_code=400, detail="No file selected.")
    ext = Path(file.filename).suffix.lower() or ".csv"
    job_id = uuid.uuid4().hex
    task_id = job_id
    saved_file = METAB_UPLOAD_DIR / f"{job_id}{ext}"
    with saved_file.open("wb") as f:
        shutil.copyfileobj(file.file, f)

    _update_metab_task(
        task_id,
        {
            "task_id": task_id,
            "status": "PENDING",
            "progress": 0,
            "message": "Task submitted.",
        },
    )
    background_tasks.add_task(
        _process_metabolomics_task,
        task_id,
        job_id,
        saved_file,
        {"epochs": epochs, "imputation": imputation, "n_clusters": n_clusters},
    )
    return {"job_id": job_id, "task_id": task_id, "message": "Metabolomics analysis started."}


@app.get("/api/metabolomics/status/{task_id}")
async def api_metabolomics_status(task_id: str) -> Dict[str, Any]:
    task = _metab_tasks.get(task_id)
    if task is None:
        return {
            "task_id": task_id,
            "status": "UNKNOWN",
            "progress": 0,
            "message": "Task not found.",
        }
    return {
        "task_id": task_id,
        "status": task.get("status", "UNKNOWN"),
        "progress": int(task.get("progress", 0)),
        "message": task.get("message", ""),
        "result": task.get("result"),
        "error": task.get("error"),
    }


@app.get("/api/metabolomics/diff_plot")
async def api_metabolomics_diff_plot(job_id: str, c1: int = 0, c2: int = 1) -> Dict[str, Any]:
    if generate_custom_diff_plot is None:
        raise HTTPException(status_code=503, detail=f"Diff-plot service unavailable: {METAB_IMPORT_ERROR}")
    output_dir = METAB_RESULTS_DIR / job_id
    if not output_dir.exists():
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")
    try:
        img_name, csv_name = generate_custom_diff_plot(str(output_dir), int(c1), int(c2))
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return {"status": "success", "image": img_name, "csv": csv_name}


@app.post("/api/enrichment/run")
async def api_enrichment_run(request: EnrichmentRequest) -> Dict[str, Any]:
    return _offline_go_enrichment(request.geneIds, request.library, request.topN)


@app.post("/api/enrichment/string_network")
async def api_string_network(request: StringNetworkRequest) -> Dict[str, Any]:
    filtered = [g for g in request.genes if float(g.expressionLevel) > float(request.minExpression)]
    if not filtered:
        raise HTTPException(status_code=400, detail="No genes pass the expression threshold.")
    filtered.sort(key=lambda x: x.expressionLevel, reverse=True)
    selected = filtered[: max(1, int(request.maxGenes))]
    symbols, _ = _normalize_gene_ids([g.geneId for g in selected])
    if not symbols:
        raise HTTPException(status_code=400, detail="No genes can be mapped to symbols for STRING query.")
    # Data URI image keeps this endpoint working offline and avoids API-rate issues.
    image = _fallback_ppi_image(symbols[: min(40, len(symbols))])
    link = _string_link(symbols[: min(200, len(symbols))], request.species, request.requiredScore, request.networkType)
    return {
        "success": True,
        "filteredCount": len(filtered),
        "mappedCount": len(symbols),
        "genesUsed": symbols,
        "image": image,
        "link": link,
    }


@app.get("/api/enrichment/batch/summary/{file_id}")
async def api_enrichment_batch_summary(file_id: str) -> Dict[str, Any]:
    df = pd.read_csv(_batch_result_file(file_id))
    return {"success": True, "data": _summary_pairs_from_matrix(df)}


@app.post("/api/enrichment/sweep/summary")
async def api_enrichment_sweep_summary(request: Dict[str, Any]) -> Dict[str, Any]:
    result = _run_sweep(
        {
            "sweepVariable": request.get("sweepVariable"),
            "range": request.get("range"),
            "fixedParamValue": request.get("fixedParamValue"),
            "smiles": request.get("smiles"),
            "cellLineId": request.get("cellLineId"),
        }
    )
    return {"success": True, "data": _summary_pairs_from_matrix(result["matrix"])}


@app.post("/api/workflow/single/metabolomics")
async def api_workflow_single_metabolomics(request: SingleMetabolomicsRequest) -> Dict[str, Any]:
    _ensure_metabolomics_available()
    if not request.expression:
        raise HTTPException(status_code=400, detail="expression is required.")

    job_id = f"single_{uuid.uuid4().hex[:12]}"
    output_dir = METAB_RESULTS_DIR / job_id
    output_dir.mkdir(parents=True, exist_ok=True)
    input_file = METAB_UPLOAD_DIR / f"{job_id}.csv"

    expression_df = _expression_to_single_matrix(request.expression, sample_id="single_0")
    expression_df.to_csv(input_file, index=True)

    try:
        flux_file, balance_file = run_scfea_training(
            input_file=str(input_file),
            output_dir=str(output_dir),
            assets_dir=str(METAB_ASSETS_DIR),
            sc_imputation=bool(request.imputation),
            epochs=int(request.epochs),
        )
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Single metabolomics workflow failed: {exc}") from exc

    flux_df = pd.read_csv(flux_file, index_col=0)
    balance_df = pd.read_csv(balance_file, index_col=0)
    flux = {str(k): float(v) for k, v in flux_df.iloc[0].to_dict().items()}
    balance = {str(k): float(v) for k, v in balance_df.iloc[0].to_dict().items()}
    return {
        "success": True,
        "jobId": job_id,
        "flux": flux,
        "balance": balance,
        "files": {
            "flux": f"/api/metabolomics/results/{job_id}/{Path(flux_file).name}",
            "balance": f"/api/metabolomics/results/{job_id}/{Path(balance_file).name}",
        },
    }


@app.post("/api/workflow/single/key_genes")
async def api_workflow_single_key_genes(request: KeyGeneRequest) -> Dict[str, Any]:
    return {"success": True, "topGenes": _key_gene_scores(request)}


@app.post("/api/workflow/batch/snapshot")
async def api_workflow_batch_snapshot(request: BatchSnapshotRequest) -> Dict[str, Any]:
    payload = {
        "sweepVariable": request.sweepVariable,
        "range": {"start": request.range.start, "end": request.range.end, "steps": request.range.steps},
        "fixedParamValue": request.fixedParamValue,
        "smiles": request.smiles,
        "cellLineId": request.cellLineId,
    }
    result = _run_sweep(payload)
    meta = result["meta"]
    matrix = result["matrix"]
    x = np.asarray([float(m["xValue"]) for m in meta], dtype=float)
    idx = int(np.abs(x - float(request.snapshotValue)).argmin())
    row = matrix.iloc[idx]
    expression = [{"geneId": str(g), "expressionLevel": float(v)} for g, v in row.to_dict().items()]
    cell_id = str(matrix.index[idx])
    return {
        "success": True,
        "pointIndex": idx,
        "cellId": cell_id,
        "xValue": float(meta[idx]["xValue"]),
        "time": float(meta[idx]["time"]),
        "dose": float(meta[idx]["dose"]),
        "expression": expression,
    }


@app.get("/api/workflow/batch/metabolomics_snapshot")
async def api_workflow_batch_metabolomics_snapshot(job_id: str = Query(...), point_index: int = Query(...)) -> Dict[str, Any]:
    flux_file = METAB_RESULTS_DIR / job_id / "predicted_flux.csv"
    if not flux_file.exists():
        raise HTTPException(status_code=404, detail=f"Metabolomics result not found for job_id={job_id}.")
    flux_df = pd.read_csv(flux_file, index_col=0)
    if point_index < 0 or point_index >= flux_df.shape[0]:
        raise HTTPException(status_code=400, detail="point_index out of range.")
    row = flux_df.iloc[int(point_index)]
    return {
        "success": True,
        "cellId": str(flux_df.index[int(point_index)]),
        "flux": {str(k): float(v) for k, v in row.to_dict().items()},
    }


@app.on_event("startup")
async def on_startup() -> None:
    if load_onnx_model is None:
        LOGGER.warning("ASCEND model import unavailable: %s", ONNX_IMPORT_ERROR)
        return
    try:
        _ensure_onnx_model()
        LOGGER.info("ASCEND model loaded.")
    except Exception as exc:
        LOGGER.warning("ASCEND model not loaded at startup: %s", exc)


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("main:app", host="0.0.0.0", port=8080, reload=True)
