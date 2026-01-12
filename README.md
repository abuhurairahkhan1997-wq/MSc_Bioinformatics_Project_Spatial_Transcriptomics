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

<img width="843" height="503" alt="image" src="https://github.com/user-attachments/assets/b6d79e52-41d9-4fdc-96fb-54a098b264b0" />








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


<<img width="560" height="550" alt="image" src="https://github.com/user-attachments/assets/012ce60b-2fbd-4b8d-811b-d1d50d931f3e" />






🤝 Contributing

Pull requests, discussions, and suggestions are welcome.
For major changes, please open an issue first to discuss the proposal.




📜 License

This project is released under the MIT License.
