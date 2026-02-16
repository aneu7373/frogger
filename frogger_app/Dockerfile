FROM python:3.11-slim

WORKDIR /app

# Install runtime deps
COPY frogger_app/requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

# Copy sources into image
COPY frogger_app/frogger/ /app/frogger/
COPY codonopt/ /app/codonopt/

# Ensure both are importable as Python packages
# (Prefer pip install if codonopt is a proper package; otherwise fallback to PYTHONPATH)
RUN python -c "import pathlib; p=pathlib.Path('/app/codonopt'); \
    has_setup=(p/'setup.py').exists() or (p/'pyproject.toml').exists(); \
    print('codonopt package installable:', has_setup)"

# Try install codonopt if possible
RUN pip install --no-cache-dir /app/codonopt || true

# Make sure python can still find the copied source trees
ENV PYTHONPATH="/app:/app/codonopt"

ENTRYPOINT ["python", "-m", "frogger.cli"]
