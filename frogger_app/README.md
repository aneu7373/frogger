FROGGER (FRagmented Orfs, Golden Gate Excision Ready)

This container bundles:
- codonopt v1.2.1 (vendored into /app/codonopt at build time)
- Frogger pipeline (in /app/frogger)

Build (from workspace root):
  docker build -t frogger:0.1 -f frogger/Dockerfile .

Run:
  docker run --rm -v "$PWD/in:/in" -v "$PWD/out:/out" frogger:0.1 \
    run --config /in/pipeline.yaml --outdir /out

Outputs are written into --outdir.
