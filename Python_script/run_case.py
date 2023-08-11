#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 16:39:19 2023

@author: maxime
"""

from design_parameters import createDesignParameters
from compute_cell_solute_interaction import computeCellSoluteInterations
from post_processing_functions import compareFields

#This script runs cell-solute simulation for a given cell type in [F7, CTX, dADSC] 
# and for a given VEGF_model in [model_1, model_2, model_3, model_4].
# It is possible to compare simulations with experiment (in the same spirit as Figure 7, 8 and 9 in the manuscript)
# for specie in [oxygen, glucose, VEGF, cell]. note that for dADSC specie='glucose' is not available due to the lack of glucose measurement
cell_type='F7'
VEGF_model='model_4'
specie='oxygen'  
design_parameters=createDesignParameters(cell_type,VEGF_model)
predictions=computeCellSoluteInterations(design_parameters)
compareFields(design_parameters['experimental_data'],predictions,cell_type,specie)

