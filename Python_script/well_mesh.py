#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math as mth
import numpy as np

#This class manage the creation of the VF mesh
class WellMesh(object):
    """ Create the mesh""" 
    #------------------
    def __init__(self,mesh_info):
        self.id=mesh_info['mesh_id']
        self.axial_length=mesh_info['axial_length']
        self.top_transverse_length=mesh_info['top_transverse_length']
        self.bottom_transverse_length=mesh_info['bottom_transverse_length']
        self.gel_length=mesh_info['gel_length']
        self.media_length=self.axial_length-self.gel_length
        self.nb_gel_cell=mesh_info['nb_gel_cell']
        self.gel_cell_length=self.gel_length/self.nb_gel_cell
        self.nb_media_cell=mesh_info['nb_media_cell']
        self.media_cell_length=self.media_length/self.nb_media_cell
       
        
        self.computeWellWallSlope()
        self.computeGelVolume()
        self.computeGelCell()
        self.computeMediaVolume()
        self.computeMediaCell()
    #------------------
        
    #------------------ 
    def computeWellWallSlope(self):
        self.well_wall_slope=(self.bottom_transverse_length-self.top_transverse_length)/self.axial_length
    #------------------
        
    #------------------ 
    def computeGelVolume(self):
        media_bottom_transverse_length=self.top_transverse_length+self.well_wall_slope*self.media_length
        self.gel_volume=-1./3.*mth.pi/self.well_wall_slope*(media_bottom_transverse_length**3-self.bottom_transverse_length**3)
    #------------------
    
    #------------------ 
    def computeMediaVolume(self):
        media_bottom_transverse_length=self.top_transverse_length+self.well_wall_slope*self.media_length
        self.media_volume=-1./3.*mth.pi/self.well_wall_slope*(self.top_transverse_length**3-media_bottom_transverse_length**3)
    #------------------
        
    #------------------
    def computeMediaCell(self):
        media_cell_volume_list=[]
        media_cell_top_transverse_length_list=[]
        media_cell_bottom_transverse_length_list=[]
        media_cell_transverse_length_list=[]
        media_cell_top_surface_list=[]
        media_cell_bottom_surface_list=[]
        media_cell_position_list=[]
        media_cell_position_top_list=[]
        media_cell_position_bottom_list=[]
        for idx in range(self.nb_media_cell):
            media_cell_position=(idx+1/2)*self.media_cell_length
            media_cell_position_top=idx*self.media_cell_length
            media_cell_position_bottom=(idx+1)*self.media_cell_length
            media_cell_top_transverse_length=self.top_transverse_length+self.well_wall_slope*(media_cell_position_top)
            media_cell_bottom_transverse_length=self.top_transverse_length+self.well_wall_slope*(media_cell_position_bottom)
            media_cell_transverse_length=self.top_transverse_length+self.well_wall_slope*(media_cell_position)
            
            media_cell_volume=-1./3.*mth.pi/self.well_wall_slope*(media_cell_top_transverse_length**3-media_cell_bottom_transverse_length**3)
            
            media_cell_top_surface=mth.pi*media_cell_top_transverse_length**2
            media_cell_bottom_surface=mth.pi*media_cell_bottom_transverse_length**2
            
            media_cell_volume_list.append(media_cell_volume)
            media_cell_top_transverse_length_list.append(media_cell_top_transverse_length)
            media_cell_bottom_transverse_length_list.append(media_cell_bottom_transverse_length)
            media_cell_transverse_length_list.append(media_cell_transverse_length)
            media_cell_top_surface_list.append(media_cell_top_surface)
            media_cell_bottom_surface_list.append(media_cell_bottom_surface)
            
            media_cell_position_list.append(media_cell_position)
            media_cell_position_top_list.append(media_cell_position_top)
            media_cell_position_bottom_list.append(media_cell_position_bottom)
        self.media_cell_volume_list=np.array(media_cell_volume_list)
        self.media_cell_top_transverse_length_list=np.array(media_cell_top_transverse_length_list)
        self.media_cell_bottom_transverse_length_list=np.array(media_cell_bottom_transverse_length_list)
        self.media_cell_transverse_length_list=np.array(media_cell_transverse_length_list)
        
        self.media_cell_top_surface_list=np.array(media_cell_top_surface_list)
        self.media_cell_bottom_surface_list=np.array(media_cell_bottom_surface_list)
        
        self.media_cell_position_list=np.array(media_cell_position_list)
        self.media_cell_position_top_list=np.array(media_cell_position_top_list)
        self.media_cell_position_bottom_list=np.array(media_cell_position_bottom_list)
    #------------------     
    
    #------------------
    def computeGelCell(self):
        gel_cell_volume_list=[]

        gel_cell_top_transverse_length_list=[]
        gel_cell_bottom_transverse_length_list=[]
        gel_cell_transverse_length_list=[]        
        
        gel_cell_top_surface_list=[]
        gel_cell_bottom_surface_list=[]
        gel_cell_position_list=[]
        gel_cell_position_top_list=[]
        gel_cell_position_bottom_list=[]
        for idx in range(self.nb_gel_cell):
            gel_cell_position=self.media_length+(idx+1/2)*self.gel_cell_length
            gel_cell_position_top=self.media_length+idx*self.gel_cell_length
            gel_cell_position_bottom=self.media_length+(idx+1)*self.gel_cell_length
            
            gel_cell_top_transverse_length=self.top_transverse_length+self.well_wall_slope*(self.media_length+idx*self.gel_cell_length)
            gel_cell_bottom_transverse_length=self.top_transverse_length+self.well_wall_slope*(self.media_length+(idx+1)*self.gel_cell_length)
            gel_cell_transverse_length=self.top_transverse_length+self.well_wall_slope*(gel_cell_position)

            gel_cell_volume=-1./3.*mth.pi/self.well_wall_slope*(gel_cell_top_transverse_length**3-gel_cell_bottom_transverse_length**3)
            
            gel_cell_top_surface=mth.pi*gel_cell_top_transverse_length**2
            gel_cell_bottom_surface=mth.pi*gel_cell_bottom_transverse_length**2
            
            gel_cell_volume_list.append(gel_cell_volume)
            gel_cell_top_transverse_length_list.append(gel_cell_top_transverse_length)
            gel_cell_bottom_transverse_length_list.append(gel_cell_bottom_transverse_length)
            gel_cell_transverse_length_list.append(gel_cell_transverse_length)
            
            gel_cell_top_surface_list.append(gel_cell_top_surface)
            gel_cell_bottom_surface_list.append(gel_cell_bottom_surface)
            
            gel_cell_position_list.append(gel_cell_position)
            gel_cell_position_top_list.append(gel_cell_position_top)
            gel_cell_position_bottom_list.append(gel_cell_position_bottom)
            
        self.gel_cell_volume_list=np.array(gel_cell_volume_list)
        self.gel_cell_top_transverse_length_list=np.array(gel_cell_top_transverse_length_list)
        self.gel_cell_bottom_transverse_length_list=np.array(gel_cell_bottom_transverse_length_list)
        self.gel_cell_transverse_length_list=np.array(gel_cell_transverse_length_list)
        
        self.gel_cell_top_surface_list=np.array(gel_cell_top_surface_list)
        self.gel_cell_bottom_surface_list=np.array(gel_cell_bottom_surface_list)
        
        self.gel_cell_position_list=np.array(gel_cell_position_list)
        self.gel_cell_position_top_list=np.array(gel_cell_position_top_list)
        self.gel_cell_position_bottom_list=np.array(gel_cell_position_bottom_list)
    #------------------    
    
            
        

