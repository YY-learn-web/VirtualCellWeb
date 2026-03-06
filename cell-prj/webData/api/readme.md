# Virtual Cell API

FastAPI backend for virtual cell prediction and drug response analysis.

## 🚀 Quick Start

### Prerequisites
- Python 3.10+
- Conda (recommended) or pip

### Installation

#### Option 1: Using Conda (Recommended)

1. **Create and activate conda environment:**
```bash
conda create -n cell-prj-env python=3.10 -y
conda activate cell-prj-env
```

2. **Install PyTorch (CPU version):**
```bash
conda install pytorch==2.0.1 torchvision==0.15.2 torchaudio==2.0.2 cpuonly -c pytorch -y
```

3. **Install PyTorch Geometric dependencies:**
```bash
pip install torch-scatter==2.1.1 torch-sparse==0.6.17 --no-build-isolation --no-cache-dir
pip install torch-geometric==2.3.1
```

4. **Install scientific computing packages:**
```bash
conda install -c conda-forge rdkit -y
pip install "numpy<2.0" pandas transformers
```

5. **Install ONNX Runtime (支持算子):**
```bash
pip install onnxruntime>=1.20.0  # 支持 ReduceL2(13) 等算子
```

6. **Install FastAPI and web dependencies:**
```bash
pip install fastapi uvicorn python-multipart
```

#### Option 2: Using pip Only

```bash
# Create virtual environment
python -m venv cell-api-env
source cell-api-env/bin/activate  # On Windows: cell-api-env\Scripts\activate

# Install PyTorch
pip install torch==2.0.1 torchvision==0.15.2 torchaudio==2.0.2

# Install PyTorch Geometric
pip install torch-scatter==2.1.1 torch-sparse==0.6.17 --no-build-isolation
pip install torch-geometric==2.3.1

# Install ONNX Runtime
pip install onnxruntime>=1.20.0  # 支持 ReduceL2(13) 等算子

# Install other dependencies
pip install "numpy<2.0" pandas transformers rdkit-pypi fastapi uvicorn python-multipart
```

### 📋 Complete Requirements

The following packages are required with specific versions:

```txt
torch==2.0.1
torchaudio==2.0.2
torchvision==0.15.2
torch-geometric==2.3.1
torch-scatter==2.1.1
torch-sparse==0.6.17
numpy<2.0
pandas
transformers
rdkit
onnxruntime>=1.20.0
fastapi
uvicorn
python-multipart
```

## 🏃‍♂️ Running the Server

1. **Activate the environment:**
```bash
conda activate cell-prj-env
```

2. **Navigate to the API directory:**
```bash
cd /path/to/webData/api
```

3. **Start the server:**
```bash
python main.py
```

4. **Access the API:**
- **API Server**: http://localhost:8080
- **Interactive Docs**: http://localhost:8080/docs
- **OpenAPI Schema**: http://localhost:8080/openapi.json

## 📚 API Endpoints

### Single Molecule Analysis
- `POST /analyze` - Analyze single molecule
- **Body**: JSON with SMILES, time, dose, cell line

### Batch Analysis
- `POST /batch` - Batch analysis with parameter sweep
- **Body**: JSON with sweep parameters

### File Upload
- `POST /upload` - Upload CSV/Excel file for batch processing
- **Form Data**: File upload with processing options

### Health Check
- `GET /` - Server status and info

## 🔧 Configuration

The server is configured with:
- **CORS**: Enabled for all origins (development mode)
- **Port**: 8080 (default)
- **Host**: 0.0.0.0 (accessible from any interface)

## 🐛 Troubleshooting

### Common Issues

1. **NumPy Compatibility Error**
   ```bash
   pip install "numpy<2.0"
   ```

2. **torch-scatter/torch-sparse Build Failures**
   ```bash
   pip install torch-scatter torch-sparse --no-build-isolation --no-cache-dir
   ```

3. **Port Already in Use**
   ```bash
   # Find and kill process using port 8080
   lsof -ti:8080 | xargs kill -9
   ```

4. **RDKit Installation Issues**
   ```bash
   # Use conda-forge for RDKit
   conda install -c conda-forge rdkit
   ```

5. **ONNX Model Loading Failed (ReduceL2 Operator)**
   - **问题**: 模型使用了 `ReduceL2(13)` 算子，但旧版本的 onnxruntime 不支持
   - **解决方案**: 
     ```bash
     # 确保使用 Python 3.10+ 环境
     conda create -n cell-prj-env python=3.10 -y
     conda activate cell-prj-env
     
     # 安装支持算子的 onnxruntime 版本
     pip install onnxruntime>=1.20.0
     ```
   - **验证**: 检查 onnxruntime 版本和可用执行提供者
     ```python
     import onnxruntime as ort
     print(f"ONNX Runtime: {ort.__version__}")
     print(f"可用执行提供者: {ort.get_available_providers()}")
     ```

### Environment Verification

Verify all packages are correctly installed:

```python
import torch
import torch_geometric
import torch_scatter
import torch_sparse
import pandas
import numpy
import transformers
import rdkit
import fastapi
import onnxruntime as ort

print("✅ All packages imported successfully!")
print(f"PyTorch: {torch.__version__}")
print(f"PyTorch Geometric: {torch_geometric.__version__}")
print(f"NumPy: {numpy.__version__}")
print(f"ONNX Runtime: {ort.__version__}")
print(f"可用执行提供者: {ort.get_available_providers()}")
```

## 🤝 Development

### Project Structure
```
api/
├── model-app.py          # Main FastAPI application
├── README.md            # This file
└── requirements.txt     # Python dependencies (if needed)
```

### Adding New Endpoints

1. Define new route in `model-app.py`
2. Add appropriate Pydantic models for request/response
3. Update this README with new endpoint documentation

## 📄 License

This project is part of the Virtual Cell AI system.

## 🆘 Support

For issues and questions:
1. Check the troubleshooting section above
2. Verify all dependencies are installed with correct versions
3. Ensure Python 3.10+ is being used (required for ONNX operator support)
4. Check that conda/virtual environment is activated
5. For ONNX model issues, ensure onnxruntime>=1.20.0 is installed