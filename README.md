better with R 4.4.0 (seurat object) 857s to build   

GitHub Multi-Omic Demo  data: https://zenodo.org/records/15162323 Datasets for ShinyCell2 Example Applications 


# 📊 ShinyCell2 Docker

This project provides a **Dockerized Shiny Server** with [ShinyCell2](https://github.com/the-ouyang-lab/ShinyCell2) preinstalled.  
It allows researchers to deploy interactive single-cell visualization dashboards directly from **Seurat** or **AnnData (`.h5ad`)** objects.

---

## 🔧 Features

- Shiny Server preconfigured with **ShinyCell2**
- R packages: `Seurat`, `ggplot2`, `ggpubr`, `ggdendro`, `bslib`
- Python support (`anndata`, `scanpy`) for `.h5ad` input
- Runs Shiny apps under `/srv/shiny-server/`
- Accessible from browser via `http://localhost:3838`

---

## 🚀 Usage

### 1. Build Docker image
```bash
docker build -t shinycell2 -f Dockfile .



