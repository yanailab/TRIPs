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
- `cycle_analysis.ipynb`
