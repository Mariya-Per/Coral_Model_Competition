# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 12:06:38 2022

@author: perepely
"""

import numpy as np
from tqdm import tqdm
import os
from netCDF4 import Dataset # may need for the output writing

from coral_model.core  import Group_parameters, Coral, Massive, Branching, Light, Flow, Temperature, Photosynthesis, \
    PopulationStates, Calcification, Morphology, \
    Dislodgement, Recruitment
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
        
        @property 
        
        def group1_const(self):
            dfbfbb
            return fdbfdbg
        
        def group2_const(self):
            fsdgbfdgbf
            return fsbfgb 
        
        
        def dc_rep(self):
            """Representative coral diameter; weighted average of base and plate diameters."""
            return (self.bc * (self.hc - self.tc) + self.dc * self.tc) / self.hc
        
    
        def dc_rep_matrix(self):
            """Reshaped representative coral diameter."""
            return RESHAPE.variable2matrix(self.dc_rep, 'space')
        
           
        def cover(self):
            """Carrying capacity.""" # it should be most likely out of here! or we define here that all the space is available for potential settling and growth and then murder them in Canopy file
            if self._cover is None:
                cover = np.ones(np.array(self.volume).shape)
                cover[self.volume == 0.] = 0. 
                return cover
    
            return self._cover
    
        @cover.setter
        def cover(self, carrying_capacity):
            """
            :param carrying_capacity: carrying capacity [m2 m-2]
            :type carrying_capacity: float, list, tuple, numpy.ndarray
            """
            carrying_capacity = RESHAPE.variable2array(carrying_capacity) # Check if this has to be here!
            if not self.volume.shape == carrying_capacity.shape:
                raise ValueError(
                    f'Shapes do not match: '
                    f'{self.volume.shape} =/= {carrying_capacity.shape}'
                )
    
            if sum(self.volume[carrying_capacity == 0.]) > 0. :
                print(
                    f'WARNING: Coral volume present where the carrying capacity is zero. This is unrealistic.'
                )
    
            self._cover = carrying_capacity
    
        @property
        def living_cover(self):
            """Living coral cover based on population states."""
            if self.pop_states is not None:
                return self.pop_states.sum(axis=2)
        
        # define those parameters here
        
        # Class canopy 
        
        
        # initiate class objects - processes that will be run here as =None
        
        
        # print smth to show that initialisation of the class worked
        
        # list of properties
        
        # list of setters
        
        # subclasses with processes as classes with their functions
    
    
    massive_const = self.read_coral_parameters(file='massive_input.txt')    
    
    branching_const = self.read_coral_parameters(file='branching_input.txt')
    
        
    def coral_potential (self, coral, year):
        
        """ This function calls for all the coral development processes and calculates 
        potential biomass production by the end of the year""" 
        
        # total laight available to the coral colony
        
        self.lme = Light(
            constants=self.constants, cor_const = self.cor_const,
            light_in=time_series_year(self.environment.light, year),
            lac=time_series_year(self.environment.light_attenuation, year),
            depth=self.hydrodynamics.water_depth)
        
        self.lme.rep_light(coral)
        

        # flow micro-environment - will be disabled most likely
        self.fme = Flow(
            constants=self.constants,
            u_current = current_vel, 
            u_wave = wave_vel, 
            h = self.hydrodynamics.water_depth, 
            peak_period = wave_per)
                    
        self.fme.velocities(coral, in_canopy=self.constants.fme) # check if actually needed
        self.fme.thermal_boundary_layer(coral) # will be disabled
                # thermal micro-environment
        self.tme = Temperature(
            constants = self.constants,
            temperature = time_series_year(self.environment.temp_kelvin, year))
                    
        self.tme.coral_temperature(coral) # also will be disabled
        
        # coral physiology
        self.phd = Photosynthesis(
            constants = self.constants,
            cor_const = self.cor_const,
            light_in = time_series_year(self.environment.light, year),
            first_year=True if i == 0 else False )
                
        self.phd.photo_rate(coral, self.environment, year)
                # population states
        self.ps = PopulationStates(cor_const = self.cor_const)
        self.ps.pop_states_t(coral)
                # calcification
        self.cr = Calcification(constants = self.constants)
        self.cr.calcification_rate(coral, time_series_year(self.environment.aragonite, year))
                  
        # morphology

                # morphological development
        self.mor = Morphology(
                    cor_const = self.cor_const,
                    calc_sum = coral.calc.sum(axis=1),
                    light_in = time_series_year(self.environment.light, year)
                )
        self.mor.update(coral)
        

    # storm damage
        if self.environment.storm_category is not None:
            tt=self.environment.storm_category
            yr=year
            stormcat = int(tt['stormcat'].values[tt.index==yr])
            if stormcat > 0:
            
            # update hydrodynamic model
                current_vel, wave_vel, wave_per = self.hydrodynamics.update(coral, stormcat)
                        # storm flow environment
                sfe = Flow(constants = self.constants,
                       u_current =current_vel, 
                       u_wave = wave_vel,
                       h = self.hydrodynamics.water_depth,
                       peak_period= wave_per)
                sfe.velocities(coral, in_canopy=self.constants.fme)
                        # storm dislodgement criterion
            self.sdc = Dislodgement(self.constants)
            sdc.update(coral)

# recruitment

                # recruitment
    rec = Recruitment(self.constants)
    rec.update(coral)

                # # export results
              
                # map-file
    self.output.update_map(coral, year)
                # his-file
    # self.output.update_his(coral, self.environment.dates[self.environment.dates.dt.year == year)


        2. both coral outputs by the end of the year are run through the space-checker
        
        
        3. enable competition, recruitment
        
        4. recalculate new biomass dimentions for both groups
        
        5. write output file and return info to the next timestep
        
    def coral_competition():
        
        """ Contains all the processes and functions, related to spatial competition
        between coral functional groups."""
        
            
            
        def space_checker (self,canopy, massive, branching):
            """ Checks whether there is space available for coral development or 
            spatial mortality should be enabled."""
            
            # Should I call here just (self, canopy, coral), or specify the call towards the 2 coral groups I have? 
            
            def area_occupied (self, coral):
                """ Total area (planar view), occupied by a single coral group. """
                
                # same here, do I need to specify the call, or as both groups inherit the properties, I can
                # keep these functions generic? 
                
                planar_area = coral.planar_area 
                
                pop_den = coral.pop_den 
                
                area_occ = planar_area * pop_den
                
                return area_occ
                
                
            area_massive = self.area_occupied(Massive)
            
            area_branching = self.area_occupied(Branching)
            
            total_occupied = area_massive + area_branching
            
            if total_occupied <= size_cell: # size_cell stands for total space available for occupation
                
            if any(p.sum(axis=1) > 1.0001 * coral.cover):
                slot_1 = np.arange(len(coral.cover))[p.sum(axis=1) > 1.0001 * coral.cover]
                slot_2 = p[p.sum(axis=1) > 1.0001 * coral.cover]
                slot_3 = coral.cover[p.sum(axis=1) > 1.0001 * coral.cover]
                print(
                    f'WARNING: Total population than carrying capacity at {slot_1}. '
                    f'\n\tPT = {slot_2}; K = {slot_3}' )
        
                
                pass
            else:
                print ("Potential coral cover higher than space available. Spatial competition enabled.")
                call for competition function     