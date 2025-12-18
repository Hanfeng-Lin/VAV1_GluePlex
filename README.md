# VAV1_GluePlex
Code for GluePlex analysis in the publication "Beyond the G-Loop: VAV1 Degradation Unveils an RT-Loop Degron for CRBN Molecular Glues"

This repository contains a set of scripts to automate the Boltz prediction and post-analysis workflow for protein-ligand docking.

## Workflow Overview

1.  **Batch Prediction**: Generate and run Boltz prediction scripts for multiple clusters in parallel.
2.  **Post-Analysis**: Collect results, align structures, extract coordinates and scores, and visualize the data.

## Files

-   `run_batch_clusters.py`: Main script to generate and execute Boltz prediction scripts.
    -   Customizable parameters: Protein sequences, ligand SMILES, diffusion samples.
    -   Parallel execution on available GPUs.
-   `post_analysis.py`: Script to process the output of the batch run.
    -   Collects PDB and JSON files.
    -   Aligns all structures to a reference (Chain A).
    -   Extracts Chain B coordinates and confidence scores (`iptm`, `pair_chains_iptm`).
    -   Generates `atom.csv` with aligned coordinates and scores.
    -   Runs UMAP visualization.
-   `umap_coord_energy.py`: Script to perform UMAP dimensionality reduction on the extracted coordinates and visualize the results.
-   `CRBN_vav1_template1.yaml`: Example YAML configuration file for Boltz.
-   `cluster*_1.pdb.cif`: Input CIF files for the batch prediction.

## Installation

To set up the environment for this workflow, run the following commands:

```bash
conda create -n boltz_analysis python=3.12
conda activate boltz_analysis
pip install "boltz[cuda]" -U
pip install umap-learn
pip install umap-learn[plot]
```

## Usage

### 1. Run Batch Predictions

Edit `run_batch_clusters.py` to set your desired protein sequences, ligand SMILES, and other parameters. Then run:

```bash
python run_batch_clusters.py
```

This will generate `run_boltz2_clusterX.py` scripts and execute them. Results will be saved in `CRBN_VAV1_template_noMSA_20runs/combined`.

### 2. Run Post-Analysis

After the batch predictions are complete, run the post-analysis script:

```bash
python post_analysis.py
```

This script will:
1.  Create a `results` directory and copy all PDB and JSON files there.
2.  Align structures and extract data to `atom.csv`.
3.  Run `umap_coord_energy.py` to generate UMAP plots.

### 3. View Results

-   `atom.csv`: Contains the aligned coordinates and confidence scores for all models.
-   `UMAP_*.html` / `UMAP_*.png`: UMAP visualization plots showing the clustering of docked poses colored by confidence scores.
