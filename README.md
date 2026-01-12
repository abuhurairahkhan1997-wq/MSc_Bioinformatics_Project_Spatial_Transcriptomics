# MSc_Bioinformatics_Project_Spatial_Transcriptomics
This repository contains all my Python and R programming code related to my dissertation for my MSc taught program in Bioinformatics.

🧬 Integrated Spatial Transcriptomics Analysis Pipeline


A unified Python + R workflow for 10x Visium spatial transcriptomics data
This repository contains a complete analysis pipeline for 10x Visium spatial transcriptomics, integrating Python (Scanpy, Squidpy) and R (BayesSpace, clusterProfiler) with optional LLM‑assisted cell‑type annotation workflows.
The pipeline supports preprocessing, spatial clustering, hybrid clustering, marker discovery, pathway enrichment, and biological validation using established cartilage‑zone markers.


🚀 Key Features


Gene-expression clustering: Leiden, Louvain, GMM, K‑means, hierarchical
Spatial clustering: BayesSpace (default, k‑means, mclust initialisation)
Hybrid clustering: image + gene expression
Marker gene discovery using ranked DE (Wilcoxon)
LLM‑ready cluster annotation prompts (JSON + TXT)
Pathway enrichment using clusterProfiler
Biological validation via PRG4 (superficial) & COL10A1 (hypertrophic)
Clean modular scripts for reproducibility


🗂 Repository Structure


<img width="910" height="703" alt="image" src="https://github.com/user-attachments/assets/e8cabc87-71da-4f03-a22c-e31877417f32" />






🔧 Installation

Python Environment

<img width="657" height="107" alt="image" src="https://github.com/user-attachments/assets/8197371e-5cd8-4a27-9b47-3c22ee53e679" />



R Environment

<img width="732" height="265" alt="image" src="https://github.com/user-attachments/assets/ee956707-6129-472d-97ac-e5388065edc3" />




📘 Workflow Summary

1. Preprocessing & Leiden Clustering (Python)

Load Visium data with Squidpy
QC, normalisation, log‑transform
Identify highly variable genes (HVGs)
PCA → neighbourhood graph → Leiden sweep (0.2–1.6)
Generate spatial PDFs of cluster maps
Compute ARI vs expert annotations

2. BayesSpace Clustering (R)

Build a SpatialExperiment object
Run BayesSpace (q = 5, d = 15)
Compare default, k-means, and mclust initialisation
Export cluster assignments

3. Marker Gene Discovery & LLM Prompt Generation (Python)

DE testing per cluster (Wilcoxon)
Filter top positive markers (logFC-based)
Generate GPTCelltype‑format prompts
Export to .txt and .json

4. Pathway Validation & Zone‑Specific Marker QC (R)

Validate PRG4 (superficial) and COL10A1 (hypertrophic)
Define background = all expressed genes
Run GO:BP enrichment with BH correction
Output spatial feature plots + dotplots




📊 Clustering Performance (ARI Scores)


<img width="459" height="452" alt="image" src="https://github.com/user-attachments/assets/48b698d9-f9c6-485f-a075-92d8b7128546" />





🤝 Contributing

Pull requests, discussions, and suggestions are welcome.
For major changes, please open an issue first to discuss the proposal.




📜 License

This project is released under the MIT License.
