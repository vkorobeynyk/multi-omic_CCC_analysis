#!/usr/bin/env bash
# Create a local conda env for the stLearn notebooks (09b, 10b) and register it as the
# Jupyter kernel "stlearn" (the kernel name both notebooks use).
#
# Usage:  bash setup_stlearn_env.sh
#
# stLearn 1.4.1 needs Python >= 3.12. torch/torchvision are pulled in by stLearn's imports
# but not used by the CCI functions, so the small CPU-only builds are installed.

set -euo pipefail

ENV_NAME="stlearn"
PY_VER="3.12"
STLEARN_VER="1.4.1"   # version the notebooks were tested with

# make `conda activate` work inside a script
eval "$(conda shell.bash hook)"

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "!! conda env '$ENV_NAME' already exists -- remove it first: conda env remove -n $ENV_NAME" >&2
  exit 1
fi

conda create -y -n "$ENV_NAME" -c conda-forge python="$PY_VER" pip
conda activate "$ENV_NAME"

# CPU-only torch first, so stLearn doesn't pull the multi-GB CUDA build
pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu

pip install "stlearn==${STLEARN_VER}" ipykernel nbclient nbformat jupyter_client

# register as Jupyter kernel "stlearn" (overwrites an existing kernelspec with that name,
# e.g. one pointing to the Singularity container)
python -m ipykernel install --user --name stlearn --display-name "stlearn"

# quick check
python - <<'PY'
import stlearn as st, scanpy as sc, anndata as ad, numba, numpy as np, pandas as pd
print("stlearn", st.__version__, "| scanpy", sc.__version__, "| anndata", ad.__version__,
      "| numba", numba.__version__, "| numpy", np.__version__, "| pandas", pd.__version__)
print("connectomeDB2020_lit mouse pairs:", len(st.tl.cci.load_lrs(["connectomeDB2020_lit"], species="mouse")))
PY

echo "Done. Kernel registered:"
jupyter kernelspec list | grep -i stlearn || true
