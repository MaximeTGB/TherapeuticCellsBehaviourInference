# TherapeuticCellsBehaviourInference
Files and Script relating to the publication "Mathematical Modelling with Bayesian Inference to Quantitatively Characterise Therapeutic Cell Behaviour in Nerve Tissue Engineering" published in the Journal of the Royal Society interface

This repository contains 3 main Folder:

-The Figure folder, which contains 1) a folder containing a .pdf of each Figure in the manuscript and 2) a folder containing the source files (.tex) to create each individual panel in each figure. Initial compilation made on Overleaf.

-The Data Folder that contains, for each cell type explored in the manuscript the measures of oxygen, glucose, VEGF and viable cell density, including error bar where available. It also contains a folder for the acellular measure of oxygen. We note that oxygen measure files vary from the other in their organisation, as the first column represents time (in hours), the second the average measure over the repeat and the third one the associated spread. For all other files the first column is the oxygen level at the air/media interface. Units for oxygen are in %O2, units for VEGF are in pg/ml, units for cell density is in MCell/ml. Glucose data are dimensionless as they represents a fraction. 

-Python Script folder that contains the code used to simulate cell-solute interactions (see comment in run_case.py)
