# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 12:06:38 2022

@author: perepely
"""

import numpy as np
from tqdm import tqdm
import os
from netCDF4 import Dataset # may need for the output writing

from coral_model.core  import Group_parameters, Coral, Massive, Branching


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
        
        
        # define those parameters here
        
        # Class canopy 
        
        
        # initiate class objects - processes that will be run here as =None
        
        
        # print smth to show that initialisation of the class worked
        
        # list of properties
        
        # list of setters
        
        # subclasses with processes as classes with their functions
    
    
    massive_const = self.read_coral_parameters(file='massive_input.txt')    
    
    branching_const = self.read_coral_parameters(file='branching_input.txt')
    
    
        
            1. inner loop - coral physiology
        lme (light micro-environment) is affected by presence of branching coral, because then
        less light can go through the angle or reflecton off the sediment - side-effect LAC    
        
    def coral_potential (self, coral, coral_constants, year):
        
        """ This function calls for all the coral development processes and calculates 
        potential biomass production by the end of the year""" 
        
        self.lme = Light(
                    constants=self.constants,
                    light_in=time_series_year(self.environment.light, year),
                    lac=time_series_year(self.environment.light_attenuation, year),
                    depth=self.hydrodynamics.water_depth
                )
                lme.rep_light(coral)
        
        
        
    def exec(self, coral, duration=None):
        """Execute simulation.

        :param coral: coral animal
        :param duration: simulation duration [yrs], defaults to None

        :type coral: Coral
        :type duration: int, optional
        """
        # auto-set duration based on environmental time-series
        if duration is None:
            duration = int(self.environment.dates.iloc[-1].year - self.environment.dates.iloc[0].year)
        years = range(int(self.environment.dates.iloc[0].year), int(self.environment.dates.iloc[0].year + duration))

        with tqdm(range((int(duration)))) as progress:
            for i in progress:

                # # environment
                progress.set_postfix(inner_loop='coral environment')
                # light micro-environment
                lme = Light(
                    constants=self.constants,
                    light_in=time_series_year(self.environment.light, years[i]),
                    lac=time_series_year(self.environment.light_attenuation, years[i]),
                    depth=self.hydrodynamics.water_depth
                )
                lme.rep_light(coral)
                # flow micro-environment
                fme = Flow(
                    constants=self.constants,
                    u_current = current_vel, 
                    u_wave = wave_vel, 
                    h = self.hydrodynamics.water_depth, 
                    peak_period = wave_per
                    )
                fme.velocities(coral, in_canopy=self.constants.fme)
                fme.thermal_boundary_layer(coral)
                # thermal micro-environment
                tme = Temperature(
                    constants = self.constants,
                    temperature = time_series_year(self.environment.temp_kelvin, years[i])
                    )
                tme.coral_temperature(coral)

                # # physiology
                progress.set_postfix(inner_loop='coral physiology')
                # photosynthetic dependencies
                phd = Photosynthesis(
                    constants = self.constants,
                    light_in = time_series_year(self.environment.light, years[i]),
                    first_year=True if i == 0 else False
                )
                phd.photo_rate(coral, self.environment, years[i])
                # population states
                ps = PopulationStates(constants = self.constants)
                ps.pop_states_t(coral)
                # calcification
                cr = Calcification(constants = self.constants)
                cr.calcification_rate(
                    coral, time_series_year(self.environment.aragonite, years[i])
                )
                # # morphology
                progress.set_postfix(inner_loop='coral morphology')
                # morphological development
                mor = Morphology(
                    constants = self.constants,
                    calc_sum = coral.calc.sum(axis=1),
                    light_in = time_series_year(self.environment.light, years[i])
                )
                mor.update(coral)

                # # storm damage
                if self.environment.storm_category is not None:
                    tt=self.environment.storm_category
                    yr=years[i]
                    stormcat = int(tt['stormcat'].values[tt.index==yr])
                    if stormcat > 0:
                        progress.set_postfix(inner_loop='storm damage')
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
                        sdc = Dislodgement(self.constants)
                        sdc.update(coral)

                # # recruitment
                progress.set_postfix(inner_loop='coral recruitment')
                # recruitment
                rec = Recruitment(self.constants)
                rec.update(coral)

                # # export results
                progress.set_postfix(inner_loop='export results')
                # map-file
                self.output.update_map(coral, years[i])
                # his-file
                self.output.update_his(coral, self.environment.dates[self.environment.dates.dt.year == years[i]])


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
                pass
            else:
                print ("Potential coral cover higher than space available. Spatial competition enabled.")
                call for competition function     