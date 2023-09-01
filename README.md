# TRIPs

### Dependencies
Analysis was done with Python v3.8.6 with the following packages:
- numpy v1.23.3
- pandas v1.3.3
- scanpy v1.7.1 (https://scanpy.readthedocs.io)
- scvi v0.9.0 (https://scvi-tools.org/)

### Guide to notebooks
There are two example workflows (see Table S1 from the manuscript for details):
- Ecoli_D1: *E. coli* MG1655 grown in LB
- Saureus_D5: *S. aureus* USA300 LAC grown in TSB

For each of these we have the following notebooks:
- `initial_processing.ipynb`: Code for doing initial data import, filtering, and scVI denoising. This also has the code to generate global correlation patterns.
- `cycle_analysis.ipynb`: Code for doing cell cycle analysis including assigning cell and gene angles, finding cycle variable genes, and aligning cell angle by the predicted point of replication.
- `promoter_distance_analysis.ipynb`: Determining relationships between expression patterns and distance from the transcriptional start site (see Figure 4 in the manuscript).
- `trip_analysis.ipynb`: Defining Transcription-Replication Interaction Profiles (TRIPs) and performing clustering on these.

In addition there is the `origin_angle_circular_model.R` file. This is a script that must be run separately in RStudio as part of the `cycle_analysis.ipynb` workflow (see notebook for details).

  There are also the following folders:
  - `count_matrices`: Raw count matrices for each library as output of the PETRI-seq processing pipeline (https://tavazoielab.c2b2.columbia.edu/PETRI-seq/).
  - `outputs`: Output files from the scripts above.
  - `samples`: Actual files used within the analysis presented for the manuscript. These are frequently used as inputs in the Jupyter notebooks to aid comparison to manuscript figures.
  - `reference`: Gene annotation files.
