#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 17:07:30 2023

@author: maxime
"""
import numpy as np

#Definitions of metablic functions used to decribe cell-solute interations following equations 2.14 to 2.18 and Table 3)
def create_metabolism_function_oxygen(parameters): 
    def metabolism_function_oxygen(concentrations):
        oxygen=concentrations[0]
        density=concentrations[3]
        metabolism=-parameters['oxygen_consumption_rate']*density*oxygen/(parameters['oxygen_concentration_half']+oxygen)
        return metabolism
    return metabolism_function_oxygen

def create_metabolism_function_glucose(parameters): 
    def metabolism_function_glucose(concentrations):
                oxygen=concentrations[0]
                glucose=concentrations[1]
                density=concentrations[3]
                metabolism=-parameters['glucose_consumption_rate']*density*glucose/(parameters['glucose_concentration_half']+glucose)*(
                            (1+parameters['glucose_anaerobic_factor']*(parameters['oxygen_concentration_half']/(oxygen+parameters['oxygen_concentration_half']))))
                return metabolism
    return metabolism_function_glucose


def create_metabolism_function_cell(parameters): 
    def growth_death_function(concentrations):
        oxygen=concentrations[0]
        glucose=concentrations[1]
        density=concentrations[3]
        growth=parameters['cell_growth_rate']*density*(1-density/parameters['cell_maximum_density'])*oxygen/(parameters['oxygen_concentration_half']+oxygen)*glucose/(parameters['glucose_concentration_half']+glucose)
        death=-density*(parameters['cell_death_rate_baseline']+parameters['cell_death_rate_oxygen']*(parameters['oxygen_concentration_half']/(oxygen+parameters['oxygen_concentration_half']))+(
                                                      +parameters['cell_death_rate_glucose']*(parameters['glucose_concentration_half']/(glucose+parameters['glucose_concentration_half']) )))
        growth_death=growth+death
        return growth_death
    return growth_death_function


def create_metabolism_function_VEGF(parameters):
    
    if parameters['VEGF_model']=="model_1":
        def metabolism_function_VEGF(concentrations):
            oxygen=concentrations[0]
            density=concentrations[3]
            threshold=oxygen<parameters['VEGF_hypoxic_threshold']
            metabolism_VEGF=density*(parameters['VEGF_secretion_rate']*(1-oxygen/parameters['VEGF_hypoxic_threshold'])*threshold)
            return metabolism_VEGF
        
    elif parameters['VEGF_model']=="model_2":
        def metabolism_function_VEGF(concentrations):
            oxygen=concentrations[0]
            density=concentrations[3]
            metabolism_VEGF=np.zeros(len(oxygen))
            idx=0
            for ox, den in zip(oxygen,density):
                if ox>parameters['VEGF_hyperoxic_threshold']:
                   metabolism_VEGF[idx]=den*parameters['VEGF_secretion_rate']
                elif ox>parameters['VEGF_hypoxic_threshold']:
                   metabolism_VEGF[idx]=den*(parameters['VEGF_secretion_rate']*(1+parameters['VEGF_factor']*((parameters['VEGF_hyperoxic_threshold']-ox)/(parameters['VEGF_hyperoxic_threshold']-parameters['VEGF_hypoxic_threshold']))**parameters['VEGF_exponent']))
                else:
                   metabolism_VEGF[idx]=den*parameters['VEGF_secretion_rate']*(1+parameters['VEGF_factor'])
                idx=idx+1
            return metabolism_VEGF
    elif parameters['VEGF_model']=="model_3":
        def metabolism_function_VEGF(concentrations):
            oxygen=concentrations[0]
            density=concentrations[3]
            
            alpha=(parameters['VEGF_secretion_rate_1']+
                   parameters['VEGF_secretion_rate_2']*parameters['cell_initial_concentration_gel']+
                   parameters['VEGF_secretion_rate_3']*parameters['cell_initial_concentration_gel']**2
                   )
            
            V=(parameters['VEGF_factor_1']+
               parameters['VEGF_factor_2']*parameters['cell_initial_concentration_gel'])
            metabolism_VEGF=density*alpha*(1+V/2+V/2*np.tanh(parameters['VEGF_exponent']*(1-oxygen/parameters['VEGF_hypoxic_threshold'])) )
            return metabolism_VEGF
    elif parameters['VEGF_model']=="model_4":
        def metabolism_function_VEGF(concentrations):
            oxygen=concentrations[0]
            density=concentrations[3]
            metabolism_VEGF=density*(parameters['VEGF_secretion_rate_alpha']*oxygen/parameters['VEGF_hypoxic_threshold']+(parameters['VEGF_secretion_rate_beta']*np.exp(-oxygen/parameters['VEGF_hypoxic_threshold']-density/parameters['VEGF_space_threshold'])))
            return metabolism_VEGF
    return metabolism_function_VEGF
    
    