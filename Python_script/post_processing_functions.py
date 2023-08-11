#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 10:09:46 2023

@author: maxime
"""
import matplotlib.pyplot as plt
import numpy as np

# This function create figure displayed predicted and experimental field for a specific cell type and a specific specie
def compareFields(measured_data,predicted_data,cell_type,specie):
      if specie=='oxygen':
        if cell_type=='F7' or cell_type=='CTX':
            keys_to_fetch={'60':['1','3','7']}
        else:   
            keys_to_fetch={'39':['3','5'],
                           '77':['5'],
                           '154':['3'],
                           '231':['5']
                          }
        k=0
        oxygen_measured=measured_data[specie] 
        for cell_density_key, oxygen_level_key_list in keys_to_fetch.items():
            for oxygen_predicted_values in predicted_data[specie][cell_density_key]:
                for oxygen_key, oxygen_values in oxygen_predicted_values.items():
                    if oxygen_key in oxygen_level_key_list:
                       oxygen_measured_values=oxygen_measured[k][1]
                       time=oxygen_measured[k][0]
                       oxygen_error_bar=oxygen_measured[k][2]
                       fig, ax=plt.subplots()
                       ax.errorbar(time,oxygen_measured_values,yerr=oxygen_error_bar)
                       ax.plot(time,oxygen_values[0:len(time)],'r-*')
                       plt.xlabel('time')
                       plt.ylabel(specie)
                       plt.title(cell_density_key+'_'+oxygen_key)
                       plt.legend(['simulation','experiment'])
                       k=k+1    
           
      else:       
           field_measured=measured_data[specie] 
           field_predicted=predicted_data[specie]
           for field_predicted_item, field_measured_item in zip(field_predicted.items(),field_measured):
               initial_cell_density_key=field_predicted_item[0]
               field_predicted_item_dict_list=field_predicted_item[1]
               
               oxygen_level_keys=[]
               field_predicted_values=[]
               for field_predicted_item_dict in field_predicted_item_dict_list:
                   oxygen_level_keys.append([key for key in field_predicted_item_dict.keys()][0])
                   field_predicted_values.append([value for value in field_predicted_item_dict.values()][0])
               values={'experiment':[field_measured_item[1],field_measured_item[2]] ,
                      'simulation':[field_predicted_values, 0*np.arange(len(field_predicted_values))]
                      }
              
               x = np.arange(len(oxygen_level_keys))  # the label locations
               width = 0.25  # t
               multiplier = 0
               
               fig, ax = plt.subplots()
               colors=['blue','red']
               for item, color in zip(values.items(),colors):
                attribute=item[0]
                value=item[1][0]
                yerr=item[1][1]
                offset = width * multiplier
                ax.bar(x + offset, value, width, label=attribute, color=color, yerr=yerr)
                multiplier += 1
               
               ax.set_ylabel(specie)
               ax.set_xticks(x)
               ax.set_xticklabels(oxygen_level_keys)
               ax.set_title(initial_cell_density_key)
               ax.legend(loc='upper right')