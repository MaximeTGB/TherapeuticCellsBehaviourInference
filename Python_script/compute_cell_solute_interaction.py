#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 12:38:41 2022

@author: maxime
"""

from well_model_manager import WellModelManager
from cell_solute_interaction_functions import create_metabolism_function_oxygen
from cell_solute_interaction_functions import create_metabolism_function_glucose
from cell_solute_interaction_functions import create_metabolism_function_cell
from cell_solute_interaction_functions import create_metabolism_function_VEGF


# This function performs simulations associated with the experimental design described in Table 2 in the manuscript
def computeCellSoluteInterations(design_parameters):
    def setupMesh(parameters):
         mesh_info={'axial_length':parameters['axial_length'],
                    'top_transverse_length':parameters['top_transverse_length'],
                       'bottom_transverse_length':parameters['bottom_transverse_length'],
                       'gel_length':parameters['gel_length'],
                       'nb_gel_cell':parameters['nb_gel_cell'],
                       'nb_media_cell':parameters['nb_media_cell'],
                       'mesh_id':'test'
                        }
         wM.setupMesh(mesh_info)
            
    def setupOxygen(parameters):
        metabolism_function_oxygen=create_metabolism_function_oxygen(parameters)
        oxygen_info={
                        'boundary_label':parameters['oxygen_bc_label'],
                         'boundary_concentration':parameters['oxygen_bc_concentration'],
                         'diffusion_coefficient_media': parameters['oxygen_diffusion_media'],
                         'diffusion_coefficient_gel': parameters['oxygen_diffusion_gel'],
                         'partition_coefficient': parameters['oxygen_partition_coefficient'],
                         'degradation_rate':parameters['oxygen_degradation_rate'],
                         'initial_concentration_media':parameters['oxygen_initial_concentration_media'],
                         'initial_concentration_gel':parameters['oxygen_initial_concentration_gel'],
                         'reaction_function':metabolism_function_oxygen,
                         'cell_solute_id': 'oxygen',
                         'dimension':parameters['oxygen_dimension']
                        }   
        wM.setupNewCellSoluteManager(oxygen_info)

    def setupGlucose(parameters):
        metabolism_function_glucose=create_metabolism_function_glucose(parameters)
        
        glucose_info={
                     'boundary_label':parameters['glucose_bc_label'],
                     'boundary_concentration':parameters['glucose_bc_concentration'],
                     'diffusion_coefficient_media': parameters['glucose_diffusion_media'],
                     'diffusion_coefficient_gel': parameters['glucose_diffusion_gel'],
                     'partition_coefficient': parameters['glucose_partition_coefficient'],
                     'degradation_rate':parameters['glucose_degradation_rate'],
                     'initial_concentration_media':parameters['glucose_initial_concentration_media'],
                     'initial_concentration_gel':parameters['glucose_initial_concentration_gel'],
                     'reaction_function':metabolism_function_glucose,
                     'cell_solute_id': 'glucose',
                     'dimension':parameters['glucose_dimension']
                    }  
        wM.setupNewCellSoluteManager(glucose_info)        

    def setupVEGF(parameters):
        metabolism_function_VEGF=create_metabolism_function_VEGF(parameters)
     
        VEGF_info={
                     'boundary_label':parameters['VEGF_bc_label'],
                     'boundary_concentration':parameters['VEGF_bc_concentration'],
                     'diffusion_coefficient_media': parameters['VEGF_diffusion_media'],
                     'diffusion_coefficient_gel': parameters['VEGF_diffusion_gel'],
                     'partition_coefficient': parameters['VEGF_partition_coefficient'],
                     'degradation_rate':parameters['VEGF_degradation_rate'],
                     'initial_concentration_media':parameters['VEGF_initial_concentration_media'],
                     'initial_concentration_gel':parameters['VEGF_initial_concentration_gel'],
                     'reaction_function':metabolism_function_VEGF,
                     'cell_solute_id': 'VEGF',
                     'dimension':parameters['VEGF_dimension']
                    } 
        wM.setupNewCellSoluteManager(VEGF_info)
        
    def setupDensity(parameters):    
        growth_death_function=create_metabolism_function_cell(parameters)
            
        density_info={
                     'boundary_label':parameters['cell_bc_label'],
                     'boundary_concentration':parameters['cell_bc_concentration'],
                     'diffusion_coefficient_media': parameters['cell_diffusion_media'],
                     'diffusion_coefficient_gel': parameters['cell_diffusion_gel'],
                     'partition_coefficient': parameters['cell_partition_coefficient'],
                     'degradation_rate':parameters['cell_degradation_rate'],
                     'initial_concentration_media':parameters['cell_initial_concentration_media'],
                     'initial_concentration_gel':parameters['cell_initial_concentration_gel'],
                     'reaction_function':growth_death_function,
                     'cell_solute_id': 'cell',
                     'dimension':parameters['cell_dimension']
                    } 
        wM.setupNewCellSoluteManager(density_info)        
     
    
    predicted_field={'oxygen':{},
                     'glucose':{},
                     'VEGF':{},
                     'cell':{}}
    for initial_cell_idx, initial_cell in zip(design_parameters['experimental_design']['initial_cell_density'],design_parameters['experimental_design']['initial_cell_density'].values()):    
        predicted_field['oxygen'][initial_cell_idx]=[]        
        predicted_field['glucose'][initial_cell_idx]=[]  
        predicted_field['VEGF'][initial_cell_idx]=[]  
        predicted_field['cell'][initial_cell_idx]=[]  
                  
        for oxygen_air_idx, oxygen_air in zip(design_parameters['experimental_design']['oxygen_concentration_air_media_interface'],design_parameters['experimental_design']['oxygen_concentration_air_media_interface'].values()):
            print(initial_cell/1e12, oxygen_air/4.21e-4)
            print(initial_cell_idx, oxygen_air_idx)
            
            parameters=design_parameters
            
            parameters['oxygen_initial_concentration_media']=design_parameters['oxygen_initial_concentration']
            parameters['oxygen_initial_concentration_gel']=design_parameters['oxygen_initial_concentration']
            parameters['oxygen_bc_label']='Dirichlet'
            parameters['oxygen_bc_concentration']=oxygen_air
            parameters['oxygen_dimension']=4.21e-4
            
            parameters['glucose_initial_concentration_media']=design_parameters['glucose_initial_concentration']
            parameters['glucose_initial_concentration_gel']=design_parameters['glucose_initial_concentration']/20
            parameters['glucose_bc_label']='Neumann'
            parameters['glucose_bc_concentration']=-1
            parameters['glucose_dimension']=4.5
            
            
            parameters['VEGF_initial_concentration_media']=0
            parameters['VEGF_initial_concentration_gel']=0
            parameters['VEGF_bc_label']='Neumann'
            parameters['VEGF_bc_concentration']=-1
            parameters['VEGF_dimension']=100e-9
            
            parameters['cell_initial_concentration_media']=0
            parameters['cell_initial_concentration_gel']=initial_cell
            parameters['cell_bc_label']='Neumann'
            parameters['cell_bc_concentration']=-1
            parameters['cell_dimension']=77e12
            
            wM=WellModelManager()
            setupMesh(parameters) 
            setupOxygen(parameters)
            setupGlucose(parameters)
            setupVEGF(parameters)
            setupDensity(parameters)
            wM.setTimestep(parameters['timestep'])
            wM.setPhysicalTime(parameters['physical_time'])
            wM.solve()
            oxygen_concentration_gel=wM.cell_solute_list['oxygen'].average_concentration_gel
            glucose_concentration_media=wM.cell_solute_list['glucose'].average_concentration_media[-1]
            VEGF_concentration_media=wM.cell_solute_list['VEGF'].average_concentration_media[-1]
            cell_concentration_gel=wM.cell_solute_list['cell'].average_concentration_gel[-1]
            t=0
            predicted_oxygen_level=[]
            for oxygen_concentration in oxygen_concentration_gel:
                if t%1800==0:
                    predicted_oxygen_level.append(oxygen_concentration/4.21e-4)
                t+=parameters['timestep']
            
            predicted_field['oxygen'][initial_cell_idx].append({oxygen_air_idx:predicted_oxygen_level})        
            predicted_field['glucose'][initial_cell_idx].append({oxygen_air_idx:glucose_concentration_media/design_parameters['glucose_initial_concentration_media']})
            predicted_field['VEGF'][initial_cell_idx].append({oxygen_air_idx:VEGF_concentration_media/1e-9})
            predicted_field['cell'][initial_cell_idx].append({oxygen_air_idx:cell_concentration_gel/1e12})

    return predicted_field




