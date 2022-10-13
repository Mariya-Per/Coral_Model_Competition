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
        TODO
        return spatial_mortality
    
    def light_competition(self):
        if Branching.hc > Massive.hc:
            has to output the light-attenuation coefficient that branching shade on massive (value)
            if no shade, that coef = 0:
                else formula
       return lac_shading_branching 
            
    def otput_canopy(self):
        make it write 1 map file but with layers for massive and branching
        
        massive = ncfile.createGroup('Massive')
        branching = ncfile.createGroup('Branching')
        total = ncfile.createGroup ('Total_coral')
        for grp in ncfile.groups.items():
            print(grp)
        
       TODO: should I have the Output writing here or in utils? Probably makes more sense to have it here, because
       in canopy file the "group" and total outputs will be generated.
       
        def initiate_map(self, coral):
        """Initiate mapping output file in which annual output covering the whole model domain is stored.

        :param coral: coral animal
        :type coral: Coral
        """
        if self._map_output is not None and any(self._map_output.values()):
            self._map_data = Dataset(self.file_name_map, 'w', format='NETCDF4')
            self._map_data.description = 'Mapped simulation data of the CoralCompetitionModel.'

            # dimensions
            self._map_data.createDimension('time', None)
            self._map_data.createDimension('nmesh2d_face', self.space)

            # variables
            t = self._map_data.createVariable('time', int, ('time',))
            t.long_name = 'year'
            t.units = 'years since 0 B.C.'

            x = self._map_data.createVariable('nmesh2d_x', 'f8', ('nmesh2d_face',))
            x.long_name = 'x-coordinate'
            x.units = 'm'

            y = self._map_data.createVariable('nmesh2d_y', 'f8', ('nmesh2d_face',))
            y.long_name = 'y-coordinate'
            y.units = 'm'

            t[:] = self.first_year
            x[:] = self.xy_coordinates[:, 0]
            y[:] = self.xy_coordinates[:, 1]
            
            # groups
            grp1 = ncfile.createGroup('Massive')
            grp2 = ncfile.createGroup('Branching')
            grp3 = ncfile.createGroup ('Total_coral')
            for grp in ncfile.groups.items():
                print(grp)
            

            # initial conditions
            if self._map_output['lme']:
                light_set = self._map_data.createVariable('Iz', 'f8', ('time', 'nmesh2d_face'))
                light_set.long_name = 'annual mean representative light-intensity'
                light_set.units = 'micro-mol photons m-2 s-1'
                light_set[:, :] = 0
            if self._map_output['fme']:
                flow_set = self._map_data.createVariable('ucm', 'f8', ('time', 'nmesh2d_face'))
                flow_set.long_name = 'annual mean in-canopy flow'
                flow_set.units = 'm s-1'
                flow_set[:, :] = 0
            if self._map_output['tme']:
                temp_set = self._map_data.createVariable('Tc', 'f8', ('time', 'nmesh2d_face'))
                temp_set.long_name = 'annual mean coral temperature'
                temp_set.units = 'K'
                temp_set[:, :] = 0

                low_temp_set = self._map_data.createVariable('Tlo', 'f8', ('time', 'nmesh2d_face'))
                low_temp_set.long_name = 'annual mean lower thermal limit'
                low_temp_set.units = 'K'
                low_temp_set[:, :] = 0

                high_temp_set = self._map_data.createVariable('Thi', 'f8', ('time', 'nmesh2d_face'))
                high_temp_set.long_name = 'annual mean upper thermal limit'
                high_temp_set.units = 'K'
                high_temp_set[:, :] = 0
            if self._map_output['pd']:
                pd_set = self._map_data.createVariable('PD', 'f8', ('time', 'nmesh2d_face'))
                pd_set.long_name = 'annual sum photosynthetic rate'
                pd_set.units = '-'
                pd_set[:, :] = 0
            if self._map_output['ps']:
                pt_set = self._map_data.createVariable('PT', 'f8', ('time', 'nmesh2d_face'))
                pt_set.long_name = 'total living coral population at the end of the year'
                pt_set.units = '-'
                pt_set[:, :] = coral.living_cover

                ph_set = self._map_data.createVariable('PH', 'f8', ('time', 'nmesh2d_face'))
                ph_set.long_name = 'healthy coral population at the end of the year'
                ph_set.units = '-'
                ph_set[:, :] = coral.living_cover

                pr_set = self._map_data.createVariable('PR', 'f8', ('time', 'nmesh2d_face'))
                pr_set.long_name = 'recovering coral population at the end of the year'
                pr_set.units = '-'
                pr_set[:, :] = 0

                pp_set = self._map_data.createVariable('PP', 'f8', ('time', 'nmesh2d_face'))
                pp_set.long_name = 'pale coral population at the end of the year'
                pp_set.units = '-'
                pp_set[:, :] = 0

                pb_set = self._map_data.createVariable('PB', 'f8', ('time', 'nmesh2d_face'))
                pb_set.long_name = 'bleached coral population at the end of the year'
                pb_set.units = '-'
                pb_set[:, :] = 0
            if self._map_output['calc']:
                calc_set = self._map_data.createVariable('calc', 'f8', ('time', 'nmesh2d_face'))
                calc_set.long_name = 'annual sum calcification rate'
                calc_set.units = 'kg m-2 yr-1'
                calc_set[:, :] = 0
            if self._map_output['md']:
                dc_set = self._map_data.createVariable('dc', 'f8', ('time', 'nmesh2d_face'))
                dc_set.long_name = 'coral plate diameter'
                dc_set.units = 'm'
                dc_set[0, :] = coral.dc

                hc_set = self._map_data.createVariable('hc', 'f8', ('time', 'nmesh2d_face'))
                hc_set.long_name = 'coral height'
                hc_set.units = 'm'
                hc_set[0, :] = coral.hc

                bc_set = self._map_data.createVariable('bc', 'f8', ('time', 'nmesh2d_face'))
                bc_set.long_name = 'coral base diameter'
                bc_set.units = 'm'
                bc_set[0, :] = coral.bc

                tc_set = self._map_data.createVariable('tc', 'f8', ('time', 'nmesh2d_face'))
                tc_set.long_name = 'coral plate thickness'
                tc_set.units = 'm'
                tc_set[0, :] = coral.tc

                ac_set = self._map_data.createVariable('ac', 'f8', ('time', 'nmesh2d_face'))
                ac_set.long_name = 'coral axial distance'
                ac_set.units = 'm'
                ac_set[0, :] = coral.ac
                
                pa_set = self._map_data.createVariable('PA', 'f8', ('time', 'nmesh2d_face'))
                pa_set.long_name = 'planar area'
                pa_set.units = 'm2'
                pa_set[0, :] = coral.planar_area     
                
                SA_set = self._map_data.createVariable('SA', 'f8', ('time', 'nmesh2d_face'))
                SA_set.long_name = 'surface area'
                SA_set.units = 'm2'
                SA_set[0, :] = coral.surface_area

                vc_set = self._map_data.createVariable('Vc', 'f8', ('time', 'nmesh2d_face'))
                vc_set.long_name = 'coral volume'
                vc_set.units = 'm3'
                vc_set[0, :] = coral.volume
            self._map_data.close()

    def update_map(self, coral, year):
        """Write data as annual output covering the whole model domain.

        :param coral: coral animal
        :param year: simulation year

        :type coral: Coral
        :type year: int
        """
        if self._map_output is not None and any(self._map_output.values()):
            self._map_data = Dataset(self.file_name_map, mode='a')

            i = int(year - self.first_year)
            self._map_data['time'][i] = year
            if self._map_output['lme']:
                self._map_data['Iz'][-1, :] = coral.light[:, -1]
            if self._map_output['fme']:
                self._map_data['ucm'][-1, :] = coral.ucm
            if self._map_output['tme']:
                self._map_data['Tc'][-1, :] = coral.temp[:, -1]
                self._map_data['Tlo'][-1, :] = coral.Tlo if len(DataReshape.variable2array(coral.Tlo)) > 1 else coral.Tlo * np.ones(self.space)
                self._map_data['Thi'][-1, :] = coral.Thi if len(DataReshape.variable2array(coral.Thi)) > 1 else coral.Thi * np.ones(self.space)
            if self._map_output['pd']:
                self._map_data['PD'][-1, :] = coral.photo_rate.mean(axis=1)
            if self._map_output['ps']:
                self._map_data['PT'][-1, :] = coral.pop_states[:, -1, :].sum(axis=1)
                self._map_data['PH'][-1, :] = coral.pop_states[:, -1, 0]
                self._map_data['PR'][-1, :] = coral.pop_states[:, -1, 1]
                self._map_data['PP'][-1, :] = coral.pop_states[:, -1, 2]
                self._map_data['PB'][-1, :] = coral.pop_states[:, -1, 3]
            if self._map_output['calc']:
                self._map_data['calc'][-1, :] = coral.calc.sum(axis=1)
            if self._map_output['md']:
                self._map_data['dc'][-1, :] = coral.dc
                self._map_data['hc'][-1, :] = coral.hc
                self._map_data['bc'][-1, :] = coral.bc
                self._map_data['tc'][-1, :] = coral.tc
                self._map_data['ac'][-1, :] = coral.ac
                self._map_data['PA'][-1, :] = coral.planar_area
                self._map_data['SA'][-1, :] = coral.surface_area
                self._map_data['Vc'][-1, :] = coral.volume

            self._map_data.close()
      
        
        

