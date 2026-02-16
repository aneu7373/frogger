FROM python:3.11-slim

WORKDIR /app

# Basic OS packages needed for downloading micromamba and running tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl bzip2 \
  && rm -rf /var/lib/apt/lists/*

# --- Install micromamba (lightweight conda) ---
ENV MAMBA_ROOT_PREFIX=/opt/mamba
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | \
    tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

# --- Create a small conda env just for CryptKeeper native deps ---
# Provides:
# - transtermhp (binary)
# - viennarna (includes Python 'RNA' bindings in this env)
RUN micromamba create -y -n ckdeps -c conda-forge -c bioconda \
    python=3.11 viennarna transtermhp \
  && micromamba clean -a -y

# Put ckdeps on PATH so CryptKeeper can find transterm, and Python can import RNA
ENV PATH=/opt/mamba/envs/ckdeps/bin:$PATH
ENV PYTHONPATH=/opt/mamba/envs/ckdeps/lib/python3.11/site-packages:$PYTHONPATH

# Provide a `transterm` command name (CryptKeeper expects exact name `transterm`)
RUN if command -v transterm >/dev/null 2>&1; then \
      echo "transterm already present"; \
    elif command -v transtermhp >/dev/null 2>&1; then \
      printf '#!/bin/sh\nexec transtermhp "$@"\n' > /usr/local/bin/transterm && \
      chmod +x /usr/local/bin/transterm && \
      echo "Created /usr/local/bin/transterm wrapper -> transtermhp"; \
    else \
      echo "ERROR: neither transterm nor transtermhp found in ckdeps env"; \
      exit 1; \
    fi

# Python deps for codonopt (pip)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Fail-fast checks
RUN python - <<'PY'
import shutil
tt = shutil.which("transterm")
assert tt, "ERROR: transterm not found on PATH"
print("transterm OK:", tt)

import RNA
print("ViennaRNA (RNA module) import OK")

import cryptkeeper
print("CryptKeeper import OK")
PY

# Copy app
COPY codonopt/ codonopt/
COPY README.md .

ENTRYPOINT ["python", "-m", "codonopt.main"]
