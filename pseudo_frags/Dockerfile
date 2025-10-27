# Use official lightweight Python base
FROM python:3.10-slim

# Install system dependencies for hdf5 and scipy
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libhdf5-dev \
        libatlas-base-dev \
        liblapack-dev \
        gfortran && \
    rm -rf /var/lib/apt/lists/*

# Create work directory
WORKDIR /app

# Copy script
COPY h5mu_to_pseudofragments.py /app/

# Install Python dependencies
RUN pip install --no-cache-dir \
    muon \
    numpy \
    scipy \
    h5py \
    anndata

# Default command
ENTRYPOINT ["python", "h5mu_to_pseudofragments.py"]
