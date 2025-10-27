docker build -t h5mu2frag .
docker run --rm \
  -v "$(pwd)/data":/app/data \
  h5mu2frag:latest \
  python /app/h5mu_to_pseudofragments.py \
  --h5mu /app/data/HT.h5mu --out /app/data/HT_pseudo.tsv

