#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

#This class manages the each specie individually (oxygen, glucose, VEGF, cell)
class WellCellSoluteManager(object):
    def __init__(self,cell_solute_info,mesh):
        self.mesh=mesh
        self.id=cell_solute_info['cell_solute_id']
        self.dimension=cell_solute_info['dimension']
        self.diffusion_coefficient_media=cell_solute_info['diffusion_coefficient_media']
        self.diffusion_coefficient_gel=cell_solute_info['diffusion_coefficient_gel']
        self.degradation_rate=cell_solute_info['degradation_rate']
        self.partition_coefficient=cell_solute_info['partition_coefficient']
        self.initial_concentration_media=cell_solute_info['initial_concentration_media']
        self.average_concentration_media=[self.initial_concentration_media]
        self.initial_concentration_gel=cell_solute_info['initial_concentration_gel']
        self.average_concentration_gel=[self.initial_concentration_gel]
        self.reaction=cell_solute_info['reaction_function']
        if 'reaction_derivative_function' in cell_solute_info:
            self.reaction_derivative=cell_solute_info['reaction_derivative_function']
        self.boundary_label=cell_solute_info['boundary_label']
        concentration_mesh=[]
        for idx in range(self.mesh.nb_media_cell):
            concentration_mesh.append(self.initial_concentration_media)
        for idx in range(self.mesh.nb_gel_cell):
            concentration_mesh.append(self.initial_concentration_gel)
            
        self.concentration_mesh=np.array(concentration_mesh)         
        self.boundary_concentration=cell_solute_info['boundary_concentration']
        self.computeEffectivePermeability()
        if(self.diffusion_coefficient_media>0 and self.diffusion_coefficient_gel>0):
            self.defineStationaryFunction()
            self.defineStationaryDerivative()
            self.computeNumericalSchemeCoefficient()
        else:
            self.coefficient_media_left=[ [0,0] for idx in range(self.mesh.nb_media_cell) ]
            self.coefficient_media_right=[ [0,0] for idx in range(self.mesh.nb_media_cell) ]
            self.coefficient_gel_left=[ [0,0] for idx in range(self.mesh.nb_gel_cell) ]
            self.coefficient_gel_right=[ [0,0] for idx in range(self.mesh.nb_gel_cell) ]
    #------------------
        
    #------------------    
    def computeEffectivePermeability(self):
        self.effective_permeability=2*self.diffusion_coefficient_gel/self.mesh.gel_cell_length
    #------------------    
        
    #------------------
    def defineStationaryFunction(self):
        def stationaryFunction1(position,transverse_length,slope,diffusion_coefficient, degradation_rate):
            A=transverse_length/slope
            B=degradation_rate/diffusion_coefficient
            if degradation_rate>0:
               value=np.exp(-np.sqrt(B)*position)/(A+position)
            else:
               value=-1/(A+position)
            return value
        
        def stationaryFunction2(position,transverse_length,slope,diffusion_coefficient, degradation_rate):
            A=transverse_length/slope
            B=degradation_rate/diffusion_coefficient
            if degradation_rate>0:
               value=np.exp(np.sqrt(B)*position)/(2*np.sqrt(B)*(A+position))
            else:
               value=1
            return value
        
        self.stationary_function_1=stationaryFunction1
        self.stationary_function_2=stationaryFunction2
    #------------------
        
    #------------------
    def defineStationaryDerivative(self):
        def stationaryDerivative1(position,transverse_length,slope,diffusion_coefficient, degradation_rate):
            A=transverse_length/slope
            B=degradation_rate/diffusion_coefficient
            if degradation_rate>0:
               value=-(1+np.sqrt(B)*(A+position))*np.exp(-np.sqrt(B)*position)/(A+position)**2
            else:
               value=1/(A+position)**2
            return value
        
        def stationaryDerivative2(position,transverse_length,slope,diffusion_coefficient, degradation_rate):
            A=transverse_length/slope
            B=degradation_rate/diffusion_coefficient
            if degradation_rate>0:
               value=(np.sqrt(B)*(A+position)-1)*np.exp(np.sqrt(B)*position)/(2*np.sqrt(B)*(A+position)**2)
            else:
               value=0
            return value
        
        self.stationary_derivative_1=stationaryDerivative1
        self.stationary_derivative_2=stationaryDerivative2
    #------------------    
        
    #------------------ 
    def computeNumericalSchemeCoefficient(self):
        coefficient_media_left=[]
        coefficient_media_right=[]
        for idx in range(self.mesh.nb_media_cell):
            if idx==0:
                if self.boundary_label=='Dirichlet':
                    coefficient_media_left.append(self.computeCoefficientAirMediaInterface())
                else:
                    coefficient_media_left.append([0,0])
                coefficient_media_right.append(self.computeCoefficientCentral(self.mesh.media_cell_length,self.mesh.media_cell_transverse_length_list[idx],
                                                                              self.diffusion_coefficient_media))
            elif idx==self.mesh.nb_media_cell-1:
                coefficient_media_left.append(self.computeCoefficientCentral(self.mesh.media_cell_length,self.mesh.media_cell_transverse_length_list[idx-1],
                                                                                     self.diffusion_coefficient_media))
                coefficient_media_right.append(self.computeCoefficientMediaGelInterface())
            else:
                
                coefficient_media_left.append(self.computeCoefficientCentral(self.mesh.media_cell_length,self.mesh.media_cell_transverse_length_list[idx-1],
                                                                                     self.diffusion_coefficient_media))
                
                coefficient_media_right.append(self.computeCoefficientCentral(self.mesh.media_cell_length,self.mesh.media_cell_transverse_length_list[idx],
                                                                                      self.diffusion_coefficient_media))
    
        coefficient_gel_left=[]
        coefficient_gel_right=[]           
        for idx in range(self.mesh.nb_gel_cell):
            if idx==0:
                coefficient_gel_left.append(self.computeCoefficientMediaGelInterface())
                coefficient_gel_right.append(self.computeCoefficientCentral(self.mesh.gel_cell_length,self.mesh.gel_cell_transverse_length_list[idx],
                                                                                    self.diffusion_coefficient_gel))
            elif idx==self.mesh.nb_gel_cell-1:
                coefficient_gel_left.append(self.computeCoefficientCentral(self.mesh.gel_cell_length,self.mesh.gel_cell_transverse_length_list[idx-1],
                                                                                   self.diffusion_coefficient_gel))
                coefficient_gel_right.append([0,0])
            else:
                coefficient_gel_left.append(self.computeCoefficientCentral(self.mesh.gel_cell_length,self.mesh.gel_cell_transverse_length_list[idx-1],
                                                                                   self.diffusion_coefficient_gel))
                coefficient_gel_right.append(self.computeCoefficientCentral(self.mesh.gel_cell_length,self.mesh.gel_cell_transverse_length_list[idx],
                                                                                    self.diffusion_coefficient_gel))
        self.coefficient_media_left=coefficient_media_left
        self.coefficient_media_right=coefficient_media_right
        self.coefficient_gel_left=coefficient_gel_left
        self.coefficient_gel_right=coefficient_gel_right
    #------------------ 
        
    #------------------
    def computeCoefficientCentral(self,cell_length,top_transverse_length,diffusion_coefficient):
        w11=self.stationary_function_1(0,top_transverse_length,
                                         self.mesh.well_wall_slope,diffusion_coefficient,
                                         self.degradation_rate)
        
        w12=self.stationary_function_2(0,top_transverse_length,
                                       self.mesh.well_wall_slope,diffusion_coefficient,
                                       self.degradation_rate)
        
        w21=self.stationary_function_1(cell_length,top_transverse_length,
                                       self.mesh.well_wall_slope,diffusion_coefficient,
                                       self.degradation_rate)
        
        w22=self.stationary_function_2(cell_length,top_transverse_length,
                                       self.mesh.well_wall_slope,diffusion_coefficient,
                                       self.degradation_rate)
                
        d1=self.stationary_derivative_1(cell_length/2,top_transverse_length,
                                        self.mesh.well_wall_slope,diffusion_coefficient,
                                        self.degradation_rate)
                
        d2=self.stationary_derivative_2(cell_length/2,top_transverse_length,
                                        self.mesh.well_wall_slope,diffusion_coefficient,self.degradation_rate)
               
        coefficient_central= [ (w22*d1-w21*d2)/(w11*w22-w12*w21),
                                (w11*d2-w12*d1)/(w11*w22-w12*w21) ] 
        return coefficient_central
    #------------------
        
    #------------------
    def computeCoefficientAirMediaInterface(self):
        w11=self.stationary_function_1(0,self.mesh.top_transverse_length,self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
        
        w12=self.stationary_function_2(0,self.mesh.top_transverse_length,self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
        
        w21=self.stationary_function_1(self.mesh.media_cell_length/2,self.mesh.top_transverse_length,self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)

        
        w22=self.stationary_function_2(self.mesh.media_cell_length/2,self.mesh.top_transverse_length,self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
                
        d1=self.stationary_derivative_1(0,self.mesh.top_transverse_length,self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
                
        d2=self.stationary_derivative_2(0,self.mesh.top_transverse_length,self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
               
        coefficient_air_media_interface= [ (w22*d1-w21*d2)/(w11*w22-w12*w21),
                                          (w11*d2-w12*d1)/(w11*w22-w12*w21) ]  
        return coefficient_air_media_interface
    #------------------
        
    #------------------
    def computeCoefficientMediaGelInterface(self):
        w11=self.stationary_function_1(0,self.mesh.media_cell_transverse_length_list[self.mesh.nb_media_cell-1],self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
        
        w12=self.stationary_function_2(0,self.mesh.media_cell_transverse_length_list[self.mesh.nb_media_cell-1],self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
        
        sigma=self.diffusion_coefficient_media/self.effective_permeability
        
        w21p1=sigma*self.stationary_derivative_1(self.mesh.media_cell_length/2,self.mesh.media_cell_transverse_length_list[self.mesh.nb_media_cell-1],self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
        
        w21p2=self.stationary_function_1(self.mesh.media_cell_length/2,self.mesh.media_cell_transverse_length_list[self.mesh.nb_media_cell-1],self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
        
        w21=w21p1+w21p2

        w21p1=sigma*self.stationary_derivative_1(self.mesh.media_cell_length/2,self.mesh.media_cell_transverse_length_list[self.mesh.nb_media_cell-1],self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
        
        w21p2=self.stationary_function_1(self.mesh.media_cell_length/2,self.mesh.media_cell_transverse_length_list[self.mesh.nb_media_cell-1],self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)

        
        w22p1=sigma*self.stationary_derivative_2(self.mesh.media_cell_length/2,self.mesh.media_cell_transverse_length_list[self.mesh.nb_media_cell-1],self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
        
        w22p2=self.stationary_function_2(self.mesh.media_cell_length/2,self.mesh.media_cell_transverse_length_list[self.mesh.nb_media_cell-1],self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
        
        w22=w22p1+w22p2
                
        f1=self.stationary_function_1(self.mesh.media_cell_length/2,self.mesh.media_cell_transverse_length_list[self.mesh.nb_media_cell-1],self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
                
        f2=self.stationary_function_2(self.mesh.media_cell_length/2,self.mesh.media_cell_transverse_length_list[self.mesh.nb_media_cell-1],self.mesh.well_wall_slope,self.diffusion_coefficient_media, self.degradation_rate)
               
        coefficient_media_gel_interface= [ (w22*f1-w21*f2)/(w11*w22-w12*w21),
                                          (w11*f2-w12*f1)/(w11*w22-w12*w21)-1]  
        return coefficient_media_gel_interface               
    
    #------------------ 
    def computeAverageConcentrationMedia(self):
        average_concentration_media=0
        for idx in range(self.mesh.nb_media_cell):
            average_concentration_media+=self.concentration_mesh[idx]*self.mesh.media_cell_volume_list[idx]
        average_concentration_media=average_concentration_media/self.mesh.media_volume
        self.average_concentration_media.append(average_concentration_media)        
    #------------------
        
    #------------------
    def computeAverageConcentrationGel(self):
        average_concentration_gel=0
        for idx in range(self.mesh.nb_gel_cell):
            g_idx=self.mesh.nb_media_cell+idx
            average_concentration_gel+=self.concentration_mesh[g_idx]*self.mesh.gel_cell_volume_list[idx]
        average_concentration_gel=average_concentration_gel/self.mesh.gel_volume
        self.average_concentration_gel.append(average_concentration_gel)      
    #------------------
    


        
        
        
