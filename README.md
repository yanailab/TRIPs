# TRIPs

### Dependencies
Analysis was done with Python v3.8.6 with the following libraries:
- numpy v1.23.3
- pandas v1.3.3
- scanpy v1.7.1 (https://scanpy.readthedocs.io)
- scvi v0.9.0 (https://scvi-tools.org/)

Code has been tested only on these library versions. No non-standard hardware is required and installation only takes minutes. We recommend using conda to create an specific environment - here all installations were done using anaconda3 (https://www.anaconda.com/download).

### Guide to notebooks
There are two example workflows (see Table S1 from the manuscript for details):
- Ecoli_D1: *E. coli* MG1655 grown in LB
- Saureus_D5: *S. aureus* USA300 LAC grown in TSB

For each of these we have the following notebooks:
- `initial_processing.ipynb`: Code for doing initial data import, filtering, and scVI denoising. This also has the code to generate global correlation patterns. The whole workbook takes ~1 hr to run.
- `cycle_analysis.ipynb`: Code for doing cell cycle analysis including assigning cell and gene angles, finding cycle variable genes, and aligning cell angle by the predicted point of replication. All analysis takes ~10 min to run (with additional time for `origin_angle_circular_model.R` (see below).
- `promoter_distance_analysis.ipynb`: Determining relationships between expression patterns and distance from the transcriptional start site (see Figure 4 in the manuscript). Takes <5 min to run.
- `trip_analysis.ipynb`: Defining Transcription-Replication Interaction Profiles (TRIPs) and performing clustering on these. Takes <5 min to run.

In addition there is the `origin_angle_circular_model.R` file. This is a script that must be run separately in RStudio as part of the `cycle_analysis.ipynb` workflow (see notebook for details). This takes ~5 min to run.

  There are also the following folders:
  - `count_matrices`: Raw count matrices for each library as output of the PETRI-seq processing pipeline (https://tavazoielab.c2b2.columbia.edu/PETRI-seq/).
  - `outputs`: Output files from the scripts above.
  - `samples`: Actual files used within the analysis presented for the manuscript. These are frequently used as inputs in the Jupyter notebooks to aid comparison to manuscript figures.
  - `reference`: Gene annotation files.

In order to make the analysis work, you **must run the initial steps of the `initial_processing.ipynb` notebook** to generate the AnnData file (`*_adata.h5ad`). This file was too large to upload directly to Github (>500 MB), so you'll need to generate it yourself. Please contact the lab if you have any issues doing so. Once you have done this, you should be able to run the notebooks individually without issues.
