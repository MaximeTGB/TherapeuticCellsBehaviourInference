#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 16:39:19 2023

@author: maxime
"""
import csv

#function that parameterises the cell solute model, values used here corresponds to the ones reported in Table 5 and Table 6 in the manuscript
def createDesignParameters(cell_type, vegf_model):
    design_parameter={}
    if cell_type=='F7':
       #mesh
       design_parameter['timestep']=180
       design_parameter['physical_time']=24*3600
       design_parameter['axial_length']=5.7e-3
       design_parameter['top_transverse_length']=3.45e-3
       design_parameter['bottom_transverse_length']=3.35e-3
       design_parameter['gel_length']=200e-6
       design_parameter['nb_gel_cell']=1
       design_parameter['nb_media_cell']=9
       #Oxygen
       design_parameter['oxygen_diffusion_media']=1.96e-9
       design_parameter['oxygen_diffusion_gel']=4.51e-10
       design_parameter['oxygen_degradation_rate']=0
       design_parameter['oxygen_partition_coefficient']=0.98
       design_parameter['oxygen_consumption_rate']=11.89e-20
       design_parameter['oxygen_concentration_half']=0.95*4.21e-4
       design_parameter['oxygen_initial_concentration']=12.54*4.21e-4
       
       #Glucose
       design_parameter['glucose_diffusion_media']=9.56e-10
       design_parameter['glucose_diffusion_gel']=2.70e-10
       design_parameter['glucose_degradation_rate']=0
       design_parameter['glucose_partition_coefficient']=1.4
       design_parameter['glucose_consumption_rate']=7.82e-18
       design_parameter['glucose_concentration_half']=0.55
       design_parameter['glucose_anaerobic_factor']=0.88
       design_parameter['glucose_initial_concentration']=4.5
       #Cell
       design_parameter['cell_diffusion_media']=0
       design_parameter['cell_diffusion_gel']=0
       design_parameter['cell_degradation_rate']=0
       design_parameter['cell_partition_coefficient']=0
       design_parameter['cell_death_rate_baseline']=0.35e-6
       design_parameter['cell_death_rate_oxygen']=0.78e-6
       design_parameter['cell_death_rate_glucose']=1.71e-6
       design_parameter['cell_growth_rate']=0.78e-6
       design_parameter['cell_maximum_density']=51.87e12
                        
       #VEGF
       design_parameter['VEGF_diffusion_media']=1.49e-10
       design_parameter['VEGF_diffusion_gel']=4.96e-11
       design_parameter['VEGF_degradation_rate']=7.05e-5
       design_parameter['VEGF_partition_coefficient']=1.07
       
       if vegf_model=='model_1':
           design_parameter['VEGF_model']='model_1'
           design_parameter['VEGF_secretion_rate']=0.05e-22
           design_parameter['VEGF_hypoxic_threshold']=1.33*4.21e-4  
           
       elif vegf_model=='model_2':
           design_parameter['VEGF_model']='model_2'
           design_parameter['VEGF_secretion_rate']=2.53e-24
           design_parameter['VEGF_factor']=0.63
           design_parameter['VEGF_exponent']=3.88
           design_parameter['VEGF_hypoxic_threshold']=0.99*4.21e-4
           design_parameter['VEGF_hyperoxic_threshold']=4.77*4.21e-4
           
       elif vegf_model=='model_3':
           design_parameter['VEGF_model']='model_3'
           design_parameter['VEGF_secretion_rate_1']=2.61e-24
           design_parameter['VEGF_secretion_rate_2']=0.04e-38
           design_parameter['VEGF_secretion_rate_3']=-0.01e-51
           design_parameter['VEGF_factor_1']=0.55
           design_parameter['VEGF_factor_2']=0.83e-14
           design_parameter['VEGF_exponent']=4.05
           design_parameter['VEGF_hypoxic_threshold']=1.08*4.21e-4
           
       elif vegf_model=='model_4':
           design_parameter['VEGF_model']='model_4'
           design_parameter['VEGF_secretion_rate_alpha']=0.19e-24
           design_parameter['VEGF_secretion_rate_beta']=0.07e-22
           design_parameter['VEGF_space_threshold']=51.44e12
           design_parameter['VEGF_hypoxic_threshold']=1.08*4.21e-4
      #Initial Boundary Conditions
       design_parameter['experimental_design']={'initial_cell_density':{'20':19.41e12,
                                                               '31':28.59e12,
                                                               '60':51.0e12},
                                       'oxygen_concentration_air_media_interface':{'1':1.07*4.21e-4,
                                                                                   '3':2.90*4.21e-4,
                                                                                   '7':6.99*4.21e-4,
                                                                                   '19':19.01*4.21e-4}
                                       }
       design_parameter['experimental_filenames']={'oxygen':['../Data/Cellular/F7/O2_concentration_60MCellml_1_percent.txt',
                                                             '../Data/Cellular/F7/O2_concentration_60MCellml_3_percent.txt',
                                                             '../Data/Cellular/F7/O2_concentration_60MCellml_7_percent.txt'],
                                                   
                                                   'glucose':['../Data/Cellular/F7/fraction_glucose_remaining_20MCellml_O2_mean_std.txt',
                                                              '../Data/Cellular/F7/fraction_glucose_remaining_31MCellml_O2_mean_std.txt',
                                                              '../Data/Cellular/F7/fraction_glucose_remaining_60MCellml_O2_mean_std.txt'],
                                                   
                                                   'VEGF':['../Data/Cellular/F7/VEGF_concentration_20MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/F7/VEGF_concentration_31MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/F7/VEGF_concentration_60MCellml_O2_mean_std.txt'
                                                           ],
                                                   
                                                   'cell':['../Data/Cellular/F7/viable_cell_density_20MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/F7/viable_cell_density_31MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/F7/viable_cell_density_60MCellml_O2_mean_std.txt',]}
       
       design_parameter['experimental_data']={'oxygen':[],'glucose':[],'VEGF':[],'cell':[],}
       for specie in design_parameter['experimental_filenames']:
           for filename in design_parameter['experimental_filenames'][specie]:
               experimental_data=loadMeasure(filename,specie)
               design_parameter['experimental_data'][specie].append(experimental_data)
               
    if cell_type=='CTX':
       #mesh
       design_parameter['timestep']=180
       design_parameter['physical_time']=24*3600
       design_parameter['axial_length']=5.7e-3
       design_parameter['top_transverse_length']=3.45e-3
       design_parameter['bottom_transverse_length']=3.35e-3
       design_parameter['gel_length']=200e-6
       design_parameter['nb_gel_cell']=1
       design_parameter['nb_media_cell']=9
       #Oxygen
       design_parameter['oxygen_diffusion_media']=1.96e-9
       design_parameter['oxygen_diffusion_gel']=4.51e-10
       design_parameter['oxygen_degradation_rate']=0
       design_parameter['oxygen_partition_coefficient']=0.98
       design_parameter['oxygen_consumption_rate']=4.95e-20
       design_parameter['oxygen_concentration_half']=0.87*4.21e-4
       design_parameter['oxygen_initial_concentration']=12.54*4.21e-4
       
       #Glucose
       design_parameter['glucose_diffusion_media']=9.56e-10
       design_parameter['glucose_diffusion_gel']=2.70e-10
       design_parameter['glucose_degradation_rate']=0
       design_parameter['glucose_partition_coefficient']=1.4
       design_parameter['glucose_consumption_rate']=2.22e-18
       design_parameter['glucose_concentration_half']=1.63
       design_parameter['glucose_anaerobic_factor']=4.61
       design_parameter['glucose_initial_concentration']=4.5
       #Cell
       design_parameter['cell_diffusion_media']=0
       design_parameter['cell_diffusion_gel']=0
       design_parameter['cell_degradation_rate']=0
       design_parameter['cell_partition_coefficient']=0
       design_parameter['cell_death_rate_baseline']=8.78e-6
       design_parameter['cell_death_rate_oxygen']=2.08e-6
       design_parameter['cell_death_rate_glucose']=0.72e-6
       design_parameter['cell_growth_rate']=0
       design_parameter['cell_maximum_density']=1e12
                        
       #VEGF
       design_parameter['VEGF_diffusion_media']=1.49e-10
       design_parameter['VEGF_diffusion_gel']=4.96e-11
       design_parameter['VEGF_degradation_rate']=7.05e-5
       design_parameter['VEGF_partition_coefficient']=1.07
       
       if vegf_model=='model_1':
           design_parameter['VEGF_model']='model_1'
           design_parameter['VEGF_secretion_rate']=0.57e-22
           design_parameter['VEGF_hypoxic_threshold']=1.38*4.21e-4  
           
       elif vegf_model=='model_2':
           design_parameter['VEGF_model']='model_2'
           design_parameter['VEGF_secretion_rate']=10.61e-24
           design_parameter['VEGF_factor']=2.43
           design_parameter['VEGF_exponent']=3.67
           design_parameter['VEGF_hypoxic_threshold']=1.01*4.21e-4
           design_parameter['VEGF_hyperoxic_threshold']=2.20*4.21e-4
           
       elif vegf_model=='model_3':
           design_parameter['VEGF_model']='model_3'
           design_parameter['VEGF_secretion_rate_1']=15.91e-24
           design_parameter['VEGF_secretion_rate_2']=-0.24e-38
           design_parameter['VEGF_secretion_rate_3']=-0.70e-51
           design_parameter['VEGF_factor_1']=3.24
           design_parameter['VEGF_factor_2']=-0.35e-14
           design_parameter['VEGF_exponent']=3.86
           design_parameter['VEGF_hypoxic_threshold']=0.83*4.21e-4
           
       elif vegf_model=='model_4':
           design_parameter['VEGF_model']='model_4'
           design_parameter['VEGF_secretion_rate_alpha']=1.12e-24
           design_parameter['VEGF_secretion_rate_beta']=2.54e-22
           design_parameter['VEGF_space_threshold']=26.67e12
           design_parameter['VEGF_hypoxic_threshold']=1.12*4.21e-4
      #Initial Boundary Conditions
       design_parameter['experimental_design']={'initial_cell_density':{'20':33.41e12,
                                                               '31':52.59e12,
                                                               '60':92.0e12},
                                       'oxygen_concentration_air_media_interface':{'1':1.11*4.21e-4,
                                                                                   '3':2.71*4.21e-4,
                                                                                   '7':7.14*4.21e-4,
                                                                                   '19':19.00*4.21e-4}
                                       }
       design_parameter['experimental_filenames']={'oxygen':['../Data/Cellular/CTX/O2_concentration_60MCellml_1_percent.txt',
                                                             '../Data/Cellular/CTX/O2_concentration_60MCellml_3_percent.txt',
                                                             '../Data/Cellular/CTX/O2_concentration_60MCellml_7_percent.txt'],
                                                   
                                                   'glucose':['../Data/Cellular/CTX/fraction_glucose_remaining_20MCellml_O2_mean_std.txt',
                                                              '../Data/Cellular/CTX/fraction_glucose_remaining_31MCellml_O2_mean_std.txt',
                                                              '../Data/Cellular/CTX/fraction_glucose_remaining_60MCellml_O2_mean_std.txt'],
                                                   
                                                   'VEGF':['../Data/Cellular/CTX/VEGF_concentration_20MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/CTX/VEGF_concentration_31MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/CTX/VEGF_concentration_60MCellml_O2_mean_std.txt'
                                                           ],
                                                   
                                                   'cell':['../Data/Cellular/CTX/viable_cell_density_20MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/CTX/viable_cell_density_31MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/CTX/viable_cell_density_60MCellml_O2_mean_std.txt',]}
       
       design_parameter['experimental_data']={'oxygen':[],'glucose':[],'VEGF':[],'cell':[],}
       for specie in design_parameter['experimental_filenames']:
           for filename in design_parameter['experimental_filenames'][specie]:
               experimental_data=loadMeasure(filename,specie)
               design_parameter['experimental_data'][specie].append(experimental_data)
     
    if cell_type=='dADSC':
       design_parameter['timestep']=180
       design_parameter['physical_time']=24*3600
       design_parameter['axial_length']=5.7e-3
       design_parameter['top_transverse_length']=4.03e-3
       design_parameter['bottom_transverse_length']=3.31e-3
       design_parameter['gel_length']=180e-6
       design_parameter['nb_gel_cell']=1
       design_parameter['nb_media_cell']=9
       #Oxygen
       design_parameter['oxygen_diffusion_media']=1.96e-9
       design_parameter['oxygen_diffusion_gel']=4.51e-10
       design_parameter['oxygen_degradation_rate']=0
       design_parameter['oxygen_partition_coefficient']=0.98
       design_parameter['oxygen_consumption_rate']=0.95e-20
       design_parameter['oxygen_concentration_half']=0.68*4.21e-4
       design_parameter['oxygen_initial_concentration']=19.2*4.21e-4
       
       #Glucose
       design_parameter['glucose_diffusion_media']=9.56e-10
       design_parameter['glucose_diffusion_gel']=2.70e-10
       design_parameter['glucose_degradation_rate']=0
       design_parameter['glucose_partition_coefficient']=1.4
       design_parameter['glucose_consumption_rate']=5.31e-18
       design_parameter['glucose_concentration_half']=0.87
       design_parameter['glucose_anaerobic_factor']=3.94
       design_parameter['glucose_initial_concentration']=4.5
       #Cell
       design_parameter['cell_diffusion_media']=0
       design_parameter['cell_diffusion_gel']=0
       design_parameter['cell_degradation_rate']=0
       design_parameter['cell_partition_coefficient']=0
       design_parameter['cell_death_rate_baseline']=7.9e-6
       design_parameter['cell_death_rate_oxygen']=4.48e-6
       design_parameter['cell_death_rate_glucose']=1.78e-6
       design_parameter['cell_growth_rate']=28.22e-6
       design_parameter['cell_maximum_density']=144e12
                        
       #VEGF
       design_parameter['VEGF_diffusion_media']=1.49e-10
       design_parameter['VEGF_diffusion_gel']=4.96e-11
       design_parameter['VEGF_degradation_rate']=7.05e-5
       design_parameter['VEGF_partition_coefficient']=1.07
       
       if vegf_model=='model_1':
           design_parameter['VEGF_model']='model_1'
           design_parameter['VEGF_secretion_rate']=1.17e-22
           design_parameter['VEGF_hypoxic_threshold']=7.5*4.21e-4  
           
       elif vegf_model=='model_2':
           design_parameter['VEGF_model']='model_2'
           design_parameter['VEGF_secretion_rate']=62.11e-24
           design_parameter['VEGF_factor']=1.50
           design_parameter['VEGF_exponent']=3.28
           design_parameter['VEGF_hypoxic_threshold']=5.44*4.21e-4
           design_parameter['VEGF_hyperoxic_threshold']=11.90*4.21e-4
           
       elif vegf_model=='model_3':
           design_parameter['VEGF_model']='model_3'
           design_parameter['VEGF_secretion_rate_1']=57.03e-24
           design_parameter['VEGF_secretion_rate_2']=2.66e-38
           design_parameter['VEGF_secretion_rate_3']=0.12e-51
           design_parameter['VEGF_factor_1']=1.92
           design_parameter['VEGF_factor_2']=0.1e-14
           design_parameter['VEGF_exponent']=2.78
           design_parameter['VEGF_hypoxic_threshold']=5.88*4.21e-4
           
       elif vegf_model=='model_4':
           design_parameter['VEGF_model']='model_4'
           design_parameter['VEGF_secretion_rate_alpha']=25.11e-24
           design_parameter['VEGF_secretion_rate_beta']=3.32e-22
           design_parameter['VEGF_space_threshold']=184e12
           design_parameter['VEGF_hypoxic_threshold']=5.94*4.21e-4
       #Initial Boundary Conditions
       design_parameter['experimental_design']={'initial_cell_density':{'39':39e12,
                                                                       '77':77e12,
                                                                       '154':154e12,
                                                                       '231':231e12,
                                                                       '385':385e12
                                                                       },
                                       'oxygen_concentration_air_media_interface':{'1':1*4.21e-4,
                                                                                   '3':3*4.21e-4,
                                                                                   '5':5*4.21e-4,
                                                                                   '10':10*4.21e-4,
                                                                                   '16':16*4.21e-4}
                                       }
       design_parameter['experimental_filenames']={'oxygen':['../Data/Cellular/dADSC/O2_concentration_39MCellml_3_percent.txt',
                                                             '../Data/Cellular/dADSC/O2_concentration_39MCellml_5_percent.txt',
                                                             '../Data/Cellular/dADSC/O2_concentration_77MCellml_5_percent.txt',
                                                             '../Data/Cellular/dADSC/O2_concentration_154MCellml_3_percent.txt',
                                                             '../Data/Cellular/dADSC/O2_concentration_231MCellml_5_percent.txt'],
                                                   
                                                   'VEGF':['../Data/Cellular/dADSC/VEGF_concentration_39MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/dADSC/VEGF_concentration_77MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/dADSC/VEGF_concentration_154MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/dADSC/VEGF_concentration_231MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/dADSC/VEGF_concentration_385MCellml_O2_mean_std.txt'],
                                                   
                                                   'cell':['../Data/Cellular/dADSC/viable_cell_density_39MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/dADSC/viable_cell_density_77MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/dADSC/viable_cell_density_154MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/dADSC/viable_cell_density_231MCellml_O2_mean_std.txt',
                                                           '../Data/Cellular/dADSC/viable_cell_density_385MCellml_O2_mean_std.txt']}

                        
       design_parameter['experimental_data']={'oxygen':[],'glucose':[],'VEGF':[],'cell':[],}
       for specie in design_parameter['experimental_filenames']:
           for filename in design_parameter['experimental_filenames'][specie]:
               experimental_data=loadMeasure(filename,specie)
               design_parameter['experimental_data'][specie].append(experimental_data)
            
        
    return design_parameter


# function that load the experimental data
def loadMeasure(filename,field):
    if field=='oxygen':
        time=[]
        measured_oxygen=[]
        errorbar_oxygen=[]
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                while '' in row:
                    row.remove('')
                time.append(float(row[0]))
                measured_oxygen.append(float(row[1]))
                if len(row)>2:
                    errorbar_oxygen.append(float(row[2]))
                else:
                    errorbar_oxygen.append(0.5)
        return [time, measured_oxygen, errorbar_oxygen]
    else:
        oxygen_level=[]
        measured_field=[]
        errorbar=[]
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                while '' in row:
                    row.remove('')
                oxygen_level.append(float(row[0]))
                measured_field.append(float(row[1]))
                errorbar.append(float(row[2]))
        return [oxygen_level, measured_field, errorbar]

