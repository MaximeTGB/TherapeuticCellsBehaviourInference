#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from well_mesh import WellMesh
from well_cell_solute_manager import WellCellSoluteManager
from scipy.sparse.linalg import spsolve
from scipy import sparse
import math as mth
import numpy as np
class WellModelManager(object):
    #------------------
    def __init__(self):
        self.id='default'
        self.solving_method='explicit'
        self.scheme='compact'
        self.zero_coupling=0
        self.cell_solute_list={}
    #------------------ 
        
    #------------------
    def setupMesh(self,mesh_info):
        self.mesh=WellMesh(mesh_info)
    #------------------
        
    #------------------
    def setupNewCellSoluteManager(self,cell_solute_info):
        cell_solute=WellCellSoluteManager(cell_solute_info,self.mesh)
        self.cell_solute_list[cell_solute.id]=cell_solute
    #------------------
    
    #------------------
    def setTimestep(self,timestep):
        self.timestep=timestep
        for cell_solute in self.cell_solute_list.values():
            cell_solute.timestep=timestep
    #------------------
        
    #------------------
    def setPhysicalTime(self,physical_time):
        self.physical_time=physical_time
        self.nb_time_iteration=int(physical_time/self.timestep)        
        for cell_solute in self.cell_solute_list.values():
            cell_solute.nb_time_iteration=self.nb_time_iteration
    #------------------

    #------------------
    def solve(self):
        self.dimension=np.array([cell_solute.dimension for cell_solute in self.cell_solute_list.values()])
        self.computeTransportMatrix()
        for time_iteration in range(self.nb_time_iteration):
            self.solution=[cs.concentration_mesh for cs in self.cell_solute_list.values()]
            self.computeTransportRightHandSide(time_iteration)
            self.computeReactionRightHandSide()
            self.computeRightHandSide()
            for cell_solute in self.cell_solute_list.values():
               cell_solute.concentration_mesh=spsolve(cell_solute.transport_matrix,cell_solute.right_hand_side)
               cell_solute.concentration_mesh=np.where(cell_solute.concentration_mesh<0, 0, cell_solute.concentration_mesh)
               cell_solute.computeAverageConcentrationMedia()
               cell_solute.computeAverageConcentrationGel()
    #------------------  
    
    #------------------ 
    def computeTransportRightHandSide(self,time_iteration):
        for cell_solute in self.cell_solute_list.values():
            transport_rhs=cell_solute.concentration_mesh
            if cell_solute.boundary_label=='Dirichlet':
               prefactor=cell_solute.diffusion_coefficient_media*self.timestep*self.mesh.media_cell_top_surface_list[0]/self.mesh.media_cell_volume_list[0]
               transport_rhs[0]=transport_rhs[0]-prefactor*cell_solute.coefficient_media_left[0][0]*cell_solute.boundary_concentration
               #print(cell_solute.boundary_concentration)
            cell_solute.transport_right_hand_side=transport_rhs
    #------------------         

    #------------------
    def computeReactionRightHandSide(self):
          for cell_solute in self.cell_solute_list.values():
              reaction_right_hand_side=self.timestep*cell_solute.reaction(self.solution)
              cell_solute.reaction_right_hand_side=reaction_right_hand_side
              
    #------------------ 

    #------------------
    def computeRightHandSide(self):
        for cell_solute in self.cell_solute_list.values():
            cell_solute.right_hand_side=cell_solute.transport_right_hand_side+cell_solute.reaction_right_hand_side
    #------------------
            
    #------------------            
    def computeTransportMatrix(self):
        for cell_solute in self.cell_solute_list.values():
            transport_matrix=np.zeros((len(cell_solute.concentration_mesh),len(cell_solute.concentration_mesh)))
            for idx in range(self.mesh.nb_media_cell):
                prefactor_top=cell_solute.diffusion_coefficient_media*self.timestep*self.mesh.media_cell_top_surface_list[idx]/self.mesh.media_cell_volume_list[idx]
                prefactor_bottom=cell_solute.diffusion_coefficient_media*self.timestep*self.mesh.media_cell_bottom_surface_list[idx]/self.mesh.media_cell_volume_list[idx]
                
                transport_matrix[idx][idx]=1+self.timestep*cell_solute.degradation_rate
                if idx==0:
                    if cell_solute.boundary_label=='Dirichlet':
                        transport_matrix[idx][idx]=transport_matrix[idx][idx]+prefactor_top*cell_solute.coefficient_media_left[idx][1]
                    transport_matrix[idx][idx]=transport_matrix[idx][idx]-prefactor_bottom*cell_solute.coefficient_media_right[idx][0]
                    transport_matrix[idx][idx+1]=-prefactor_bottom*cell_solute.coefficient_media_right[idx][1]
                elif idx==self.mesh.nb_media_cell-1:
                    transport_matrix[idx][idx-1]=prefactor_top*cell_solute.coefficient_media_left[idx][0]
                    transport_matrix[idx][idx]=transport_matrix[idx][idx]+prefactor_top*cell_solute.coefficient_media_left[idx][1]
                    
                    prefactor_gel_media_interface=cell_solute.effective_permeability*self.timestep*self.mesh.media_cell_bottom_surface_list[idx]/self.mesh.media_cell_volume_list[idx]
                    transport_matrix[idx][idx]=transport_matrix[idx][idx]+prefactor_gel_media_interface*cell_solute.coefficient_media_right[idx][0]
                    transport_matrix[idx][idx+1]=prefactor_gel_media_interface*cell_solute.partition_coefficient*cell_solute.coefficient_media_right[idx][1]
                else:
                    transport_matrix[idx][idx-1]=prefactor_top*cell_solute.coefficient_media_left[idx][0]
                    transport_matrix[idx][idx]=transport_matrix[idx][idx]+prefactor_top*cell_solute.coefficient_media_left[idx][1]
                    transport_matrix[idx][idx]=transport_matrix[idx][idx]-prefactor_bottom*cell_solute.coefficient_media_right[idx][0]
                    transport_matrix[idx][idx+1]=-prefactor_bottom*cell_solute.coefficient_media_right[idx][1]
            
            for idx in range(self.mesh.nb_gel_cell):
                g_idx=self.mesh.nb_media_cell+idx
                transport_matrix[g_idx][g_idx]=1+self.timestep*cell_solute.degradation_rate
                prefactor_top=cell_solute.diffusion_coefficient_gel*self.timestep*self.mesh.gel_cell_top_surface_list[idx]/self.mesh.gel_cell_volume_list[idx]
                prefactor_bottom=cell_solute.diffusion_coefficient_gel*self.timestep*self.mesh.gel_cell_bottom_surface_list[idx]/self.mesh.gel_cell_volume_list[idx]
                if self.mesh.nb_gel_cell==1:
                    prefactor_gel_media_interface=cell_solute.effective_permeability*self.timestep*self.mesh.gel_cell_top_surface_list[idx]/self.mesh.gel_cell_volume_list[idx]
                    transport_matrix[g_idx][g_idx]=transport_matrix[g_idx][g_idx]-prefactor_gel_media_interface*cell_solute.partition_coefficient*cell_solute.coefficient_gel_left[idx][1]
                    transport_matrix[g_idx][g_idx-1]=-prefactor_gel_media_interface*cell_solute.coefficient_gel_left[idx][0]
                else:
                    if idx==0:
                        prefactor_gel_media_interface=cell_solute.effective_permeability*self.timestep*self.mesh.gel_cell_top_surface_list[idx]/self.mesh.gel_cell_volume_list[idx]
                        transport_matrix[g_idx][g_idx]=transport_matrix[g_idx][g_idx]-prefactor_gel_media_interface*cell_solute.partition_coefficient*cell_solute.coefficient_gel_left[idx][1]
                        transport_matrix[g_idx][g_idx-1]=-prefactor_gel_media_interface*cell_solute.coefficient_gel_left[idx][0]
                        transport_matrix[g_idx][g_idx]=transport_matrix[g_idx][g_idx]-prefactor_bottom*cell_solute.coefficient_gel_right[idx][0]
                        transport_matrix[g_idx][g_idx+1]=-prefactor_bottom*cell_solute.coefficient_gel_right[idx][1]
                    elif idx==self.mesh.nb_gel_cell-1:
                        transport_matrix[g_idx][g_idx-1]=prefactor_top*cell_solute.coefficient_gel_left[idx][0]
                        transport_matrix[g_idx][g_idx]=transport_matrix[g_idx][g_idx]+prefactor_top*cell_solute.coefficient_gel_left[idx][1]
                    else:
                        transport_matrix[g_idx][g_idx-1]=prefactor_top*cell_solute.coefficient_gel_left[idx][0]
                        transport_matrix[g_idx][g_idx]=transport_matrix[g_idx][g_idx]+prefactor_top*cell_solute.coefficient_gel_left[idx][1]
                        transport_matrix[g_idx][g_idx]=transport_matrix[g_idx][g_idx]-prefactor_bottom*cell_solute.coefficient_gel_right[idx][0]
                        transport_matrix[g_idx][g_idx+1]=-prefactor_bottom*cell_solute.coefficient_gel_right[idx][1]

                        
                            
            csr_transport_matrix=sparse.csr_matrix(transport_matrix)
            cell_solute.transport_matrix=csr_transport_matrix
    #------------------

    #------------------
#########
    #------------------           
        
        
      
        
