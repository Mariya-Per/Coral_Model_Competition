# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 12:06:38 2022

@author: perepely
"""

import numpy as np
from tqdm import tqdm
import os
from netCDF4 import Dataset # may need for the output writing

from coral_model.core  import Coral, Massive, Branching
from coral_model.hydrodynamics import Delft3D,Transect #,Reef0D,Reef1D
from coral_model.environment import Constants, Environment

# canopy object

class Canopy:
    """ Canopy object, representing the space available for colonisation"""
    
    def __init__(self):
        """
        description of parameters that are required for canopy initialisation
        
        Canopy describes and governs space on the reef. So, for initialisation 
        I need the x, y and bathymetry? 
        """
        
        # carrying capacity
        self._carrying_capacity = None
        
        # competition process        
        self.competition = None

        @property 
        def dc_rep(self):
            """Representative coral diameter; weighted average of all coral diameters."""
            return dc_rep
        
        @setter
        def dc_rep (self, Massive, Branching):
            return (Massive.pop_den * Massive.dc + Branching.pop_den * Branching.dc) / (Massive.pop_den + Branching.pop_den)
        
        @property
        def veg_den(self):
            return Massive.pop_den + Branching.pop_den
        
        # do I need to have property and then setter for it separately? and I need to make sure at every timestep it is returned to hydrodynamics
        
        @property
        def total_pop_state(self):
            return total_pop_state
        
        @property
        def dc_rep_matrix(self):
            """Reshaped representative coral diameter."""
            return RESHAPE.variable2matrix(self.dc_rep, 'space')
        
        
        # TODO: resolve the shape of the carrying capacity
        def carrying_capacity(self):
            """Carrying capacity."""
            if self._carrying_capacity is None:
                carrying_capacity = np.ones(np.array(self.volume).shape) # should not be related to coral shape # also, it is actually set somewhere from hydrodynamics domain
                carrying_capacity[self.volume == 0.] = 0. 
                return carrying_capacity
    
            return self._carrying_capacity
    
        @cover.setter
        def carrying_capacity(self, carrying_capacity):
            """
            :param carrying_capacity: carrying capacity [m2 m-2]
            :type carrying_capacity: float, list, tuple, numpy.ndarray
            """
            carrying_capacity = RESHAPE.variable2array(carrying_capacity) # Check if this has to be here!
            if not Coral.volume.shape == carrying_capacity.shape:
                raise ValueError(
                    f'Shapes do not match: '
                    f'{self.volume.shape} =/= {carrying_capacity.shape}'
                )
    
            if sum(Coral.volume[carrying_capacity == 0.]) > 0. :
                print(
                    f'WARNING: Coral volume present where the carrying capacity is zero. This is unrealistic.'
                )
    
            self._carrying_capacity = carrying_capacity
            
            # What is the point of this? In the property it is already set that carrying capacity took the dimensions from coral volume. So why do we compare it to the coral volume again? 
    
        @property
        def living_cover(self):
            """Living coral cover based on population states."""
            if self.total_pop_state is not None:
                return self.total_pop_state.sum(axis=2)
            
        @property
        def coral_cover(self):
            """Total coral cover in the cell"""
            self.
        
    # call for the coral groups (or the output alresady from loop.py)
    
    

        
    def coral_potential (self, coral, year):
        
        """ This function calls for all the coral development processes and calculates 
        potential biomass production by the end of the year""" 
        
       # Can I actually call for something, that was calculated in the part of loop.exec function?
        
        # Or I call for 
        massive = Massive(inputs)
        
        branching = Branching(inputs)
        
        Massive.p0[i]
        Branching.p0[i]
        
        # does it make sense to have recruitment inside Coral? Because it depends on the free space available
        # Or I can still code it in there, like def recruitment (space_available) and when calling, then do:   x = Massive.recruitment(carrying_capacity - total_cover)


        2. both coral outputs by the end of the year are run through the space-checker
        
        
        3. enable competition, recruitment
        
        4. recalculate new biomass dimentions for both groups
        
        5. write output file and return info to the next timestep
        
    def spatial_competition(self, bbbb):
        
        """ Contains all the processes and functions, related to spatial competition
        between coral functional groups."""
              
                
        area_massive = Massive.planar_area * Massive.pop_den
            
        area_branching = Branching.planar_area * Branching.pop_den
            
        total_occupied = area_massive + area_branching
            
                
        if total_occupied >= 1.00 * self.carrying_capacity: # can I also have just if total_occupied >= self.carrying_capacity
            print(f'WARNING: Potential total coral cover higher than space available. Spatial competition enabled.')
            total_occupied *= self.spatial_mortality()
        else:
            free_space = self.carrying_capacity - total_occupied
            Massive.recruitment_update(free_space)
            Branching.recruitment_update(free_space)
        return dvndkbv d
    
    def spatial_mortality(self):
        
        

