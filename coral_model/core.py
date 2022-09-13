"""
coral_model - core

@author: Gijs G. Hendrickx
@contributor: Peter M.J. Herman
"""

import numpy as np
import pandas as pd
from scipy.optimize import newton

from coral_model.utils import DataReshape, CoralOnly

# # data formatting -- to be reformatted in the model simulation
RESHAPE = DataReshape()

# coral object
class Coral:
    """Coral object, representing one coral type."""

    def __init__(self, constants, cor_type,dc, hc, bc, tc, ac, species_constant):
        """
        :param dc: diameter coral plate [m]
        :param hc: coral height [m]
        :param bc: diameter coral base [m]
        :param tc: thickness coral plate [m]
        :param ac: axial distance corals [m]
        :param species_constant: species constant [-]

        :type dc: float, list, tuple, numpy.ndarray
        :type hc: float, list, tuple, numpy.ndarray
        :type bc: float, list, tuple, numpy.ndarray
        :type tc: float, list, tuple, numpy.ndarray
        :type ac: float, list, tuple, numpy.ndarray
        :type species_constant: float
        """
# =============================================================================
#         With initialisation set minimum dimensions of corals that can already
#        reproduce - data available in CTDb in terms of planar area
# =============================================================================
        self.constants = constants

        self.dc = RESHAPE.variable2array(dc)
        self.hc = RESHAPE.variable2array(hc)
        self.bc = RESHAPE.variable2array(bc)
        self.tc = RESHAPE.variable2array(tc)
        self.ac = RESHAPE.variable2array(ac)

        self.Csp = species_constant #make into the constants list

        self._cover = None  # check what does it mean exactly here 

        # initiate environmental working objects
        # > light micro-environment
        self.light = None
        #self.Bc = None
        # > flow micro-environment
        self.ucm = None
        self.um = None
        self.delta_t = None
        # > thermal micro-environment
        self.dTc = None
        self.temp = None
        # > photosynthesis
        self.photo_rate = None
        self.Tlo = None
        self.Thi = None
        # > population states
        self.pop_states = None
        self.p0 = None
        # np.array([
        #     self.cover, np.zeros(self.cover.shape), np.zeros(self.cover.shape), np.zeros(self.cover.shape),
        # ])  # What was this and why is it commented out? - creating the p0 array ith the population states as 0s, but why is it disabled
        # > calcification
        self.calc = None

    def __repr__(self):
        """Development representation."""
        return f'Morphology({self.dc}, {self.hc}, {self.bc}, {self.bc}, {self.ac})'

    def __str__(self):
        """Print representation."""
        return (f'Coral morphology with: dc = {self.dc} m; hc = {self.hc} ;'
                f'bc = {self.bc} m; tc = {self.tc} m; ac = {self.ac} m')

    @property
    def dc_rep(self):
        """Representative coral diameter; weighted average of base and plate diameters."""
        return (self.bc * (self.hc - self.tc) + self.dc * self.tc) / self.hc

    @property
    def rf(self):
        """Form ratio: height-to-diameter ratio."""
        return self.hc / self.dc

    @property
    def rp(self):
        """Plate ratio: base-to-diameter ratio."""
        return self.bc / self.dc

    @property
    def rs(self):
        """Spacing ratio: diameter-to-axial distance ratio."""
        return self.dc / self.ac
    
    @property
    def planar_area(self):
        """Planar area of the coral.
        Can be assigned from:
            1) the map - check how to do it
            2) dependant on diameter """
    
        # for now I will make only the diameter dependency
        planar_area = self.constants.dc_pa_coef * ((self.dc)** self.constants.dc_pa_exp) 
        return planar_area 
    
    @property
    def surface_area(self):
        """Active area of the coral. Basically, its surface area"""
        
        SA = self.constants.sa_dc_coef * (self.dc ** self.constants.sa_dc_exp)
        return SA

    
    @property
    # changed the volume function - assigned the output coral_volume
    
    #change recalculation based on the planar area - volume or diameter volume relationship - Check this statement,
    # maybe it is better to have both and one updates/recalculates another? As for branching coral we have volume
    # between branches, filled with water, so cannot be directly related to just geomentry
    
    def volume(self):
        """Coral volume."""
        # coral_volume = .25 * np.pi * ((self.hc - self.tc) * self.bc ** 2 + self.tc * self.dc ** 2)
        coral_volume = self.constants.pa_vol_coef * ((self.planar_area)**self.constants.pa_vol_exp)
        return coral_volume

    @volume.setter # what is the difference? And which volume does it call then? 
    # also function update morphology does not update coral volume - was disabled last year, subject of check - it does
    
    # NOTE: 
    # What was not changed but needs a lot of attention and causes the numerical crashes - spacing 
    # of corals ac.
    # If we do probabilistic occupation of canopy and coexistance of branching and massive coral volume under,
    # this ac has to be adjusted to the space of coral "stem" base diameter with 95% confidence interval - as corals
    # of different species do not like other branches/stems too close and release some substances to discourage overgrowth
    
    def volume(self, coral_volume):
        """
        :param coral_volume: coral volume [m3]
        :type coral_volume: float, int, list, tuple, np.ndarray
        """
        self.update_morphology(coral_volume, rf=self.rf, rp=self.rp, rs=self.rs) 
        
        # NOTE: Do not mix the general volume and volume update function due to competition.
        # This volume here basically updates the volume based on the ratios - that are stable now, but were dynamic due to flow and competition

    def update_morphology(self, coral_volume, rf, rp, rs):
        """Update the coral morphology based on updated coral volume and morphological ratios.

        :param coral_volume: coral volume [m3]
        :param rf: form ratio [-]
        :param rp: plate ratio [-]
        :param rs: spacing ratio [-]

        :type coral_volume: float, numpy.ndarray
        :type rf: float, numpy.ndarray
        :type rp: float, numpy.ndarray
        :type rs: float, numpy.ndarray
        """

        def vc2dc(coral_volume, rf, rp):
            """Coral volume to coral plate diameter."""
            dc = ((4. * coral_volume) /
                  (np.pi * rf * rp * (1. + rp - rp ** 2))) ** (1. / 3.)
            return dc

        def vc2hc(coral_volume, rf, rp):
            """Coral volume to coral height."""
            hc = ((4. * coral_volume * rf ** 2) / 
                  (np.pi * rp * (1. + rp - rp ** 2))) ** (1. / 3.)
            return hc

        def vc2bc(coral_volume, rf, rp):
            """Coral volume > diameter of the base."""
            bc = ((4. * coral_volume * rp ** 2) / 
                  (np.pi * rf * (1. + rp - rp ** 2))) ** (1. / 3.)
            return bc

        def vc2tc(coral_volume, rf, rp):
            """Coral volume > thickness of the plate."""
            tc = ((4. * coral_volume * rf ** 2 * rp ** 2) / 
                  (np.pi * (1. + rp - rp ** 2))) ** (1. / 3.)
            return tc

        def vc2ac(coral_volume, rf, rp, rs):
            """Coral volume > axial distance."""
            ac = (1. / rs) * ((4. * coral_volume) / 
                              (np.pi * rf * rp * (1. + rp - rp ** 2))) ** (1. / 3.)
            return ac

        # # update morphology
        self.dc = vc2dc(coral_volume, rf, rp)
        self.hc = vc2hc(coral_volume, rf, rp)
        self.bc = vc2bc(coral_volume, rf, rp)
        self.tc = vc2tc(coral_volume, rf, rp)
        self.ac = vc2ac(coral_volume, rf, rp, rs)

    @property
    def dc_matrix(self):
        """Reshaped coral plate diameter."""
        return RESHAPE.variable2matrix(self.dc, 'space') # check what exactly is "space" here and does it correlate with canopy?

    @property
    def hc_matrix(self):
        """Reshaped coral height."""
        return RESHAPE.variable2matrix(self.hc, 'space')

    @property
    def bc_matrix(self):
        """Reshaped coral base diameter."""
        return RESHAPE.variable2matrix(self.bc, 'space')

    @property
    def tc_matrix(self):
        """Reshaped coral plate thickness."""
        return RESHAPE.variable2matrix(self.tc, 'space')

    @property
    def ac_matrix(self):
        """Reshaped axial distance."""
        return RESHAPE.variable2matrix(self.ac, 'space')
    
    @property 
    def planar_area_matrix(self):
        """Reshaped planar area."""
        return RESHAPE.variable2matrix(self.planar_area, 'space')
    
    @property 
    def surface_area_matrix(self):
        """Reshaped planar area."""
        return RESHAPE.variable2matrix(self.surface_area, 'space')
    
    @property
    def dc_rep_matrix(self):
        """Reshaped representative coral diameter."""
        return RESHAPE.variable2matrix(self.dc_rep, 'space')

    @property
    def as_vegetation_density(self):
        """Translation from coral morphology to (vegetation) density."""
        #would it really work with having coral volumes co-existing in the same canopy space? 

        def function(dc_rep, ac):
            return (2 * dc_rep) / (ac ** 2)

        return CoralOnly().in_space(
            coral=self,
            function=function,
            args=(self.dc_rep, self.ac)
        )

    @property
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

    def initiate_spatial_morphology(self, cover=None):
        """Initiate the morphology based on the on set of morphological dimensions and the coral cover. This method
        contains a catch that it can only be used to initiate the morphology, and cannot overwrite existing spatial
        heterogeneous morphology definitions.

        :param cover: custom coral cover definition, defaults to None
        :type cover: None, numpy.ndarray
        """
        if cover is not None:
            cover = RESHAPE.variable2array(cover)   # ! Check the function! 
            if not cover.shape[0] == RESHAPE.space:
                msg = f'Spatial dimension of cover does not match: {cover.shape} =/= {RESHAPE.space}.'
                raise ValueError(msg)
        else:
            cover = np.ones(RESHAPE.space)

        self.p0 = np.array([
            cover,
            np.zeros(cover.shape),
            np.zeros(cover.shape),
            np.zeros(cover.shape),
            ]).transpose()
        

        self.dc = cover * self.dc
        self.hc = cover * self.hc
        self.bc = cover * self.bc
        self.tc = cover * self.tc
        self.ac = cover * self.ac


# # biophysical processes

class Massive(Coral):
    """ Class of Massive corals. Inherits all the initialisation properties and functions of parent super-class Coral.
    Here some class-specific functions can be added. """
    def __init__(self):
            super().__init__()
            

class Branching(Coral):
    """ Class of Branching corals. Inherits all the properties of class Coral"""
    def __init__(self):
            super().__init__()            
            

class Light:
    """Light micro-environment."""

    def __init__(self, constants, cor_typ, light_in, lac, depth):
        """Light micro-environment.

        :param light_in: incoming light-intensity at the water-air interface [u mol photons m-2 s-1]
        :param lac: light-attenuation coefficient [m-1]
        :param depth: water depth [m]

        :type light_in: float, list, tuple, numpy.ndarray
        :type lac: float, list, tuple, numpy.ndarray
        :type depth: float, list, tuple, numpy.ndarray
        """
        self.I0 = RESHAPE.variable2matrix(light_in, 'time')
        self.Kd = RESHAPE.variable2matrix(lac, 'time') # also there for turbidity proxy
        self.h = RESHAPE.variable2matrix(depth, 'space')                                          # CHECK!
        self.constants = constants
        
      
    

    # def biomass(self, constants, coral):
    #     """Coral biomass; as surface.
        
    #     :param coral: coral animal
    #     :type coral: Coral
    #     """
    #     self.constants = constants
    #     # base_section = self.base_light(coral)
    #     # coral.Bc = np.pi * (
    #     #         .25 * coral.dc_matrix ** 2 +
    #     #         coral.dc_matrix * coral.tc_matrix +
    #     #         coral.bc_matrix * base_section
    #     # )
        
    #     coral.Bc = self.constants.sa_dc_coef * (coral.dc_matrix ** self.constants.sa_dc_exp)
    #     return coral.Bc

    def base_light(self, coral):
        """Section of coral base receiving light.

        :param coral: coral animal
        :type coral: Coral
        """
        # # spreading of light
        theta = self.light_spreading(coral)

        # # coral base section
        base_section = coral.hc_matrix - coral.tc_matrix - (
                (coral.dc_matrix - coral.bc_matrix) / (2. * np.tan(.5 * theta))
        )
        # no negative lengths
        base_section[base_section < 0] = 0

        return base_section

    def light_spreading(self, coral):
        """Spreading of light as function of depth.

        :param coral: coral animal
        :type coral: Coral
        """
        return self.constants.theta_max * np.exp(
            -self.Kd * (self.h - coral.hc_matrix + coral.tc_matrix))
    
     

    def side_correction(self, coral):
        """Correction of the light-intensity on the sides of the coral object.

        :param coral: coral animal
        :type coral: Coral
        """
        # # spreading of light
        theta = self.light_spreading(coral)

        # # correction factor
        return np.sin(.5 * theta)
        
    def rep_light(self,coral):
        
        base_section = self.base_light(coral)
        # # light catchment per coral section
        # top of plate
        # top = .25 * np.pi * coral.dc_matrix ** 2 * self.I0 * np.exp(
        #    -self.Kd * (self.h - coral.hc_matrix))
        
        top = coral.planar_area_matrix * self.I0 * np.exp(
            -self.Kd * (self.h - coral.hc_matrix))     # formula is the same, 
        # but I introduced planar area earlier
        
        # For the rest, according to Kaniewska, 2014 - there is within colony light attenuation
        # different for branches direction
        
        # at the same time, massive corals do not have the branches and much extra SA,
        # So within colony LAC is not needed for them 
        
        # If I introduce within-colony Light attemuation, do I need side correction in this case? 
        # Probably yes, but we do not know how branches are aligned
        
        # ToDo: surface_area is actually biomass in this model, still has to be changed to the equation of dependency on planar area or volume
        
        side1 = ((coral.surface_area_matrix - coral.planar_area_matrix)/2.0) * self.I0*np.exp(
            -self.Kd*(self.h - coral.hc_matrix))* self.side_correction(coral) *np.exp(-self.constants.coral_Kd_1*(coral.tc_matrix/2.))
        
        side2 = ((coral.surface_area_matrix - coral.planar_area_matrix)/2.0) * self.I0*np.exp(
            -self.Kd*(self.h - coral.hc_matrix))* self.side_correction(coral) *np.exp(-self.constants.coral_Kd_2*(coral.tc_matrix/2.))
        
        
        
        # side1 = (coral.surface_area_matrix - coral.planar_area_matrix - base_section)* 0.5 * self.I0 * np.exp(
        #     -1.0 *self.Kd * (self.h - coral.hc_matrix)) * np.exp(-1.*self.constants.coral_Kd_1 * coral.tc_matrix)* self.side_correction(coral)
        
        # side2 = (coral.surface_area_matrix - coral.planar_area_matrix - base_section)* 0.5 * self.I0 * np.exp(
        #     -1.0 * self.Kd * (self.h - coral.hc_matrix)) * np.exp(-1.0* self.constants.coral_Kd_2 * coral.tc_matrix)* self.side_correction(coral)
        
        # # side of plate
        # side = (np.pi * coral.dc_matrix * self.I0) / self.Kd * (
        #         np.exp(-self.Kd * (self.h - coral.hc_matrix)) -
        #         np.exp(-self.Kd * (self.h - coral.hc_matrix + coral.tc_matrix))
        # ) * self.side_correction(coral)
        
        # # side of base
        # base = (np.pi * coral.bc_matrix * self.I0) / self.Kd * (
        #         np.exp(-self.Kd * (self.h - base_section)) -
        #         np.exp(-self.Kd * self.h) 
        # ) * self.side_correction(coral)
        # total
        total = top + side1 + side2

        def averaged_light(total_light, surface_area):
            """Averaged light-intensity."""
            return total_light / surface_area
        
        # Still good and a validation point for average light intensity dispersion 
        # over the colony (Kaniewska 2014, colony averaged PAR)

        coral.light = CoralOnly().in_spacetime(
            coral=coral,
            function=averaged_light,
            args=(total, coral.surface_area),
            no_cover_value=self.I0 * np.exp(-self.Kd * self.h)) # check the args, so uses total from rep_light
        



class Flow:
    """Flow micro-environment."""

    def __init__(self, constants, u_current, u_wave, h, peak_period):
        """
        :param u_current: current flow velocity [m s-1]
        :param u_wave: wave flow velocity [m s-1]
        :param h: water depth [m]
        :param peak_period: peak wave period [s]

        :type u_current: float, list, tuple, numpy.ndarray
        :type u_wave: float, list, tuple, numpy.ndarray
        :type h: float, list, tuple, numpy.ndarray
        :type peak_period: float, list, tuple, numpy.ndarray
        """
        self.uc = RESHAPE.variable2array(u_current)
        self.uw = RESHAPE.variable2array(u_wave)
        self.h = RESHAPE.variable2matrix(h, 'space')                                    # CHECK!
        self.Tp = RESHAPE.variable2array(peak_period)
        self.active = False if u_current is None and u_wave is None else True
        self.constants = constants

# =============================================================================
# Check the Porous Media Theory - Thesis, page 27! for the flow microenvironment
# =============================================================================
    @property
    def uc_matrix(self):
        """Reshaped current flow velocity."""
        return RESHAPE.variable2matrix(self.uc, 'space')  # CHECK!

    @property
    def uw_matrix(self):
        """Reshaped wave flow velocity."""
        return RESHAPE.variable2matrix(self.uw, 'space') # CHECK!

    def velocities(self, coral, in_canopy=True):
        """In-canopy flow velocities, and depth-averaged flow velocities.

        :param coral: coral animal
        :param in_canopy: determine in-canopy flow (or depth-averaged), defaults to True

        :type coral: Coral
        :type in_canopy: bool, optional
        """
        if self.active:
            alpha_w = np.ones(self.uw.shape)
            alpha_c = np.ones(self.uc.shape)
            if in_canopy:
                idx = coral.volume > 0
                for i in idx:
                    alpha_w[i] = self.wave_attenuation(self.constants,
                        coral.dc_rep[i], coral.hc[i], coral.ac[i],
                        self.uw[i], self.Tp[i], self.h[i], wac_type = 'wave'
                    )
                    alpha_c[i] = self.wave_attenuation(self.constants,
                        coral.dc_rep[i], coral.hc[i], coral.ac[i],
                        self.uc[i], 1e3, self.h[i], wac_type = 'current'
                    )
            coral.ucm = self.wave_current(alpha_w, alpha_c)
            coral.um = self.wave_current()
        else:
            coral.ucm = 9999 * np.ones(RESHAPE.space)   #CHECK!

    def wave_current(self, alpha_w=1, alpha_c=1):
        """Wave-current interaction.

        :param alpha_w: wave-attenuation coefficient, defaults to 1
        :param alpha_c: current-attenuation coefficient, defaults to 1

        :type alpha_w: float, list, tuple, numpy.ndarray, optional
        :type alpha_c: float, list, tuple, numpy.ndarray, optional

        :return: wave-current interaction
        :rtype: float, numpy.ndarray
        """
        
        # Why don't we put alpha-w and alpha_c as constants as self.constants.alpha_w ?
        # Or is it recalculated in coupling within the hydrodynamics.py file? Check!
        
        # Didn't find it anywhere else, but there is wac calculation below
        
        # Or why don't we use the WAC that is calculated below with iterations?
        
        return np.sqrt(
            (alpha_w * self.uw) ** 2 + (alpha_c * self.uc) ** 2 +
            2 * alpha_w * self.uw * alpha_c * self.uc *
            np.cos(self.constants.wcAngle)
        )
        # Alternative
        return np.sqrt(
            (self.constants.alpha_w * self.uw) ** 2 + (self.constants.alpha_c * self.uc) ** 2 +
            2 * self.constants.alpha_w * self.uw * self.constants.alpha_c * self.uc *
            np.cos(self.constants.wcAngle)
        )

    @staticmethod
    def wave_attenuation(constants, diameter, height, distance, velocity, period, depth, wac_type):
        """Wave-attenuation coefficient.

        :param diameter: representative coral diameter [m]
        :param height: coral height [m]
        :param distance: axial distance [m]
        :param velocity: flow velocity [m s-1]
        :param period: wave period [s]
        :param depth: water depth [m]
        :param wac_type: type of wave-attenuation coefficient [-]

        :type diameter: float
        :type height: float
        :type distance: float
        :type velocity: float
        :type depth: float
        :type depth: float
        :type wac_type: str
        """
        # TODO: Split this method in one solely focusing on the wave attenuation coefficient;
        #  and one implementing this method to dynamically determine the drag coefficient.
        #  Thus, reformat this method as in coral_model_v0.
        # # input check
        types = ('current', 'wave')
        if wac_type not in types:
            msg = f'WAC-type {wac_type} not in {types}.'
            raise ValueError(msg)

        # # function and derivative definitions
        def function(beta):
            """Complex-valued function to be solved, where beta is the complex representation of the wave-attenuation
            coefficient.
            """
            # components
            shear = (8. * above_motion) / (3. * np.pi * shear_length) * (abs(1. - beta) * (1. - beta))
            drag = (8. * above_motion) / (3. * np.pi * drag_length) * (abs(beta) * beta)
            inertia = 1j * beta * ((constants.Cm * lambda_planar) / (1. - lambda_planar))
            # combined
            f = 1j * (beta - 1.) - shear + drag + inertia
            # output
            return f

        def derivative(beta):
            """Complex-valued derivative to be used to solve the complex-valued function, where beta is the complex
            representation of the wave-attenuation coefficient.
            """
            # components
            shear = ((1. - beta) ** 2 / abs(1. - beta) - abs(1. - beta)) / shear_length
            drag = (beta ** 2 / abs(beta) + beta) / drag_length
            inertia = 1j * (constants.Cm * lambda_planar) / (1. - lambda_planar)
            # combined
            df = 1j + (8. * above_motion) / (3. * np.pi) * (- shear + drag) + inertia
            # output
            return df

        # # parameter definitions
        # geometric parameters
        planar_area = .25 * np.pi * diameter ** 2 # TODO: planar area = planar area
        frontal_area = diameter * height
        total_area = .5 * distance ** 2
        lambda_planar = planar_area / total_area
        lambda_frontal = frontal_area / total_area
        shear_length = height / (constants.Cs ** 2)
        # # calculations
        wac = 1.
        if depth > height:
            # initial iteration values
            above_flow = velocity
            drag_coefficient = 1.
            # iteration
            for k in range(int(constants.maxiter_k)):
                drag_length = (2 * height * (1 - lambda_planar)) / (drag_coefficient * lambda_frontal)
                above_motion = (above_flow * period) / (2 * np.pi)
                if wac_type == 'wave':
                    # noinspection PyTypeChecker
                    wac = abs(newton(
                        function, x0=complex(.1, .1), fprime=derivative,
                        maxiter=constants.maxiter_aw
                    ))
                elif wac_type == 'current':
                    x = drag_length / shear_length * (height / (depth - height) + 1)
                    wac = (x - np.sqrt(x)) / (x - 1)
                else:
                    raise ValueError(
                        f'WAC-type ({wac_type}) not in {types}.'
                    )
                porous_flow = wac * above_flow
                constricted_flow = (1 - lambda_planar) / (1 - np.sqrt(
                    (4 * lambda_planar) / (constants.psi * np.pi)
                )) * porous_flow
                reynolds = (constricted_flow * diameter) / constants.nu
                new_drag = 1 + 10 * reynolds ** (-2. / 3)
                if abs((new_drag - drag_coefficient) / new_drag) <= constants.err:
                    break
                else:
                    drag_coefficient = float(new_drag)
                    above_flow = abs(
                        (1 - constants.numericTheta) * above_flow +
                        constants.numericTheta * (
                                depth * velocity - height * porous_flow
                        ) / (depth - height)
                    )

                if k == constants.maxiter_k:
                    print(
                        f'WARNING: maximum number of iterations reached '
                        f'({constants.maxiter_k})'
                    )

        return wac

    def thermal_boundary_layer(self, coral):
        """Thermal boundary layer.

        :param coral: coral animal
        :type coral: Coral
        """
        if self.active and self.constants.tme:
            delta = self.velocity_boundary_layer(self.constants,coral)
            coral.delta_t = delta * ((self.constants.alpha / self.constants.nu) ** (1 / 3))

    @staticmethod
    def velocity_boundary_layer(constants, coral):
        """Velocity boundary layer.

        :param coral: coral animal
        :type coral: Coral
        """
        def boundary_layer(rd, nu, cf, ucm):
            """Thickness velocity boundary layer."""
            return (rd * nu) / (np.sqrt(cf) * ucm)

        return CoralOnly().in_space(
            coral=coral,
            function=boundary_layer,
            args=(constants.rd, constants.nu, constants.Cf, coral.ucm)
        )


class Temperature:
    def __init__(self, constants, temperature):
        """
        Thermal micro-environment.

        Parameters
        ----------
        temperature : numeric
            Temperature of water [K].
        """
        self.T = RESHAPE.variable2matrix(temperature, 'time')
        self.constants = constants

    def coral_temperature(self, coral):
        """Coral temperature.

        :param coral: coral animal
        :type coral: Coral
        """
        # Why is it in Kelvin? Why we have to transform Celcius to Kelvins just for this function? 
        
        if self.constants.tme:
            delta_t = RESHAPE.variable2matrix(coral.delta_t, 'space') #CHECK!
            coral.dTc = (
                    (delta_t * self.constants.ap) / (self.constants.k * self.constants.K0) *
                    coral.light
            )
            coral.temp = self.T + coral.dTc
        else:
            coral.temp = self.T


class Photosynthesis:
    """Photosynthesis."""

    def __init__(self, constants, light_in, first_year):
        """
        Photosynthetic efficiency based on photosynthetic dependencies.

        :param light_in: incoming light-intensity at the water-air interface [umol photons m-2 s-1]
        :param first_year: first year of the model simulation

        :type light_in: float, list, tuple, numpy.ndarray
        :type first_year: bool
        """
        self.I0 = RESHAPE.variable2matrix(light_in, 'time')
        self.first_year = first_year

        self.pld = 1
        self.ptd = 1
        self.pfd = 1
        # How does it correlate with initial coral input file, where pfd = False????
        
        self.constants = constants
        

    def photo_rate(self, coral, environment, year):
        """Photosynthetic efficiency.

        :param coral: coral animal
        :param environment: environmental conditions
        :param year: year of simulation

        :type coral: Coral
        :type environment: Environment
        :type year: int
        """
        # components
        self.light_dependency(coral, 'qss')
        self.thermal_dependency(coral, environment, year)
        self.flow_dependency(coral)

        # combined
        coral.photo_rate = self.pld * self.ptd * self.pfd

    def light_dependency(self, coral, output):
        """Photosynthetic light dependency.

        :param coral: coral animal
        :param output: type of output

        :type coral: Coral
        :type output: str
        """
        
        def photo_acclimation(x_old, param):
            """Photo-acclimation."""
            # input check
            params = ('Ik', 'rETRmax')
            if param not in params:
                message = f'{param} not in {params}.'
                raise ValueError(message)

            # parameter definitions
            x_max = self.constants.ik_max if param == 'Ik' else self.constants.rETR_max_qss
            beta_x = self.constants.betaI if param == 'Ik' else self.constants.beta_rETR
            
            # Why in text Table 4.1 Pm_max = 3.96, but default in the model input and environment it is 1.0?


            # calculations
            xs = x_max * (coral.light / self.I0) ** beta_x
            if output == 'qss':
                return xs
            elif output == 'new':
                return xs + (x_old - xs) * np.exp(-self.constants.iota)
            
            # but "new" is never used? - formula from Anthony, H-G(2003) for adaptation to new conditions
            
            # coral.light - colony averaged light, should be consistent with Kaniewska 2014 light!

        # # parameter definitions
        if output == 'qss':
            ik = photo_acclimation(None, 'Ik')
            rETR_max_z = photo_acclimation( None, 'rETRmax')
        else:
            msg = f'Only the quasi-steady state solution is currently implemented; use key-word \'qss\'.'
            raise NotImplementedError(msg)

# TODO: Make separate modes for Pmax and rETRmax. Note: Calcification constant should be different then as well
        # def photo_acclimation(x_old, param):
        #     """Photo-acclimation."""
        #     # input check
        #     params = ('Ik', 'Pmax')
        #     if param not in params:
        #         message = f'{param} not in {params}.'
        #         raise ValueError(message)

        #     # parameter definitions
        #     x_max = self.constants.ik_max if param == 'Ik' else self.constants.pm_max
        #     beta_x = self.constants.betaI if param == 'Ik' else self.constants.betaP
            
        #     # Why in text Table 4.1 Pm_max = 3.96, but default in the model input and environment it is 1.0?


        #     # calculations
        #     xs = x_max * (coral.light / self.I0) ** beta_x
        #     if output == 'qss':
        #         return xs
        #     elif output == 'new':
        #         return xs + (x_old - xs) * np.exp(-self.constants.iota)
            
        #     # coral.light - colony averaged light, should be consistent with Kaniewska 2014 light!

        # # # parameter definitions
        # if output == 'qss':
        #     ik = photo_acclimation(None, 'Ik')
        #     p_max = photo_acclimation( None, 'Pmax')
        # else:
        #     msg = f'Only the quasi-steady state solution is currently implemented; use key-word \'qss\'.'
        #     raise NotImplementedError(msg)
            
     


        # # calculations
        # old Platt 1976 equation without photoinhibition
        # self.pld = p_max * (np.tanh(coral.light / ik) - np.tanh(self.constants.Icomp * self.I0 / ik))
        """ Function based on Platt, 1980 to calculate photosynthetic efficiency
            P_max = set for coral group 
            Note: With transplantation P_max can change!
            
            alpha - characterises photochemical reaction, the slope of the initial part of the P-I curve
            beta - characterises photoinhibition, the slope of the last decreasing part of the P-I curve 
            
            for branching coral beta = 0 
            
            Should have 2 modes for coral here, so that: 
                massive coral alpha = self.constants.alpha (constant value, but Pmax can adjust)
                branching coral alpha = rETRmax/ik (because I decide to have them with flexible Ik, but 
                as their photoinhibition beta=0, then rETRmax observed = rETRmax(z))
    
            """
            # TODO :  make 2 modes for alpha estimation, depending on the coral group
                 
        # Massive coral
        # alpha=self.constants.alpha
        # phi = self.constants.phi - commented out, as they will be constant for both coral groups
        
        # branching coral
        #alpha = rETR_max/ik # disabled for massive now
        
        # cannot use it like this for massive, as it will have Ik at 300+, 
        # while it should be rETRmax observed/ik        
            # beta = 0 for branching but exists for massive corals. Larger beta means higher photoinhibition
            
            # a = (alpha*I)/Pmax
            # in the code I = coral.light
            
        a = (self.constants.alpha * coral.light)/rETR_max_z
            
        b = (self.constants.phi * coral.light)/rETR_max_z   
        
        self.pld = rETR_max_z  * (1.0 - np.exp(-a)) * np.exp(-b)
        
        
    

    def thermal_dependency(self, coral, env, year):
        """Photosynthetic thermal dependency.

        :param coral: coral animal
        :param env: environmental conditions
        :param year: year of simulation

        :type coral: Coral
        :type env: Environment
        :type year: int
        """

        def thermal_acc():
            """Thermal-acclimation."""
            if self.constants.tme:
                if self.first_year:
                    env.tmeMMMmin = pd.DataFrame(
                        data=pd.concat(
                            [env.temp_mmm['min']] * RESHAPE.space,
                            axis=1
                        ).values,
                        columns=[np.arange(RESHAPE.space)]
                    ) + coral.dTc
                    env.tmeMMMmax = pd.DataFrame(
                        data=pd.concat(
                            [env.temp_mmm['max']] * RESHAPE.space,
                            axis=1
                        ).values,
                        columns=[np.arange(RESHAPE.space)]
                    ) + coral.dTc
                else:
                    env.tmeMMMmin[env.tmeMMM.index == year] += coral.dTc
                    env.tmeMMMmax[env.tmeMMm.index == year] += coral.dTc

                mmm_min = env.tmeMMMmin[np.logical_and(
                    env.tmeMMM.index < year,
                    env.tmeMMM.index >= year - int(self.constants.nn / coral.Csp)
                )]
                m_min = mmm_min.mean(axis=0)
                s_min = mmm_min.std(axis=0)

                mmm_max = env.tmeMMMmax[np.logical_and(
                    env.tmeMMM.index < year,
                    env.tmeMMM.index >= year - int(self.constants.nn / coral.Csp)
                )]
                m_max = mmm_max.mean(axis=0)
                s_max = mmm_max.std(axis=0)
            else:
                mmm = env.temp_mmm[np.logical_and(
                    env.temp_mmm.index < year,
                    env.temp_mmm.index >= year - int(self.constants.nn / coral.Csp)
                )]
                m_min, m_max = mmm.mean(axis=0)
                s_min, s_max = mmm.std(axis=0)

            coral.Tlo = m_min - self.constants.k_var * s_min
            coral.Thi = m_max + self.constants.k_var * s_max

        def adapted_temp():
            """Adapted temperature response."""

            def spec():
                """Specialisation term."""
                return 4e-4 * np.exp(-.33 * (delta_temp - 10))

            response = -(coral.temp - coral.Tlo) * ((coral.temp - coral.Tlo) ** 2 - delta_temp ** 2)
            temp_cr = coral.Tlo - (1 / np.sqrt(3)) * delta_temp
            try:
                if self.constants.tme:
                    response[coral.temp <= temp_cr] = -(
                            (2 / (3 * np.sqrt(3))) * delta_temp[coral.temp <= temp_cr] ** 3
                    )
                else:
                    response[coral.temp <= temp_cr] = -(
                            (2 / (3 * np.sqrt(3))) * delta_temp ** 3
                    )
            except TypeError:
                if coral.temp <= temp_cr:
                    response = (2 / (3 * np.sqrt(3))) * delta_temp ** 3

            return response * spec()

        def thermal_env():
            """Thermal envelope."""
            return np.exp((self.constants.Ea / self.constants.R) * (1 / 300 - 1 / temp_opt))

        # # parameter definitions
        thermal_acc()
        delta_temp = coral.Thi - coral.Tlo
        temp_opt = coral.Tlo + (1 / np.sqrt(3)) * delta_temp

        # # calculations
        f1 = adapted_temp()
        f2 = thermal_env()
        self.ptd = f1 * f2

    def flow_dependency(self, coral):
        """Photosynthetic flow dependency.

        :param coral: coral animal
        :type coral: Coral
        """
        if self.constants.pfd:
            pfd = self.constants.pfd_min + (1 - self.constants.pfd_min) * np.tanh(
                2 * coral.ucm / self.constants.ucr
            )
            self.pfd = RESHAPE.variable2matrix(pfd, 'space')
        else:
            self.pfd = 1


class PopulationStates:
    """Bleaching response following the population dynamics."""  # check the competition for space here
    # TODO: Check this class; incl. writing tests

    def __init__(self, constants, dt=1):
        """Population dynamics.

        :param dt: time step [yrs], defaults to one
        :type dt: float, optional
        """
        self.dt = dt
        self.constants = constants

    def pop_states_t(self, coral):
        """Population dynamics over time.

        :param coral: coral animal
        :type coral: Coral
        """
        coral.pop_states = np.zeros((*RESHAPE.spacetime, 4))
        for n in range(RESHAPE.time):
            photosynthesis = np.zeros(RESHAPE.space)
            photosynthesis[coral.cover > 0.] = coral.photo_rate[coral.cover > 0., n] #21_09 have changed coral.cover>0 to .0. 
            coral.pop_states[:, n, :] = self.pop_states_xy(coral, photosynthesis)
            coral.p0[coral.cover > 0., :] = coral.pop_states[coral.cover > 0., n, :]

    def pop_states_xy(self, coral, ps):
        """Population dynamics over space.

        :param coral: coral animal
        :param ps: photosynthetic rate

        :type coral: Coral
        :type ps: numpy.ndarray
        """
        p = np.zeros((RESHAPE.space, 4))
        # # calculations
        # growing conditions
        # > bleached pop.      # ps>0. here represents ps>tsh that is the value of the bleaching treshold light and 1. where 1.0 is a number, not column reference
        p[ps > 0. , 3] = coral.p0[ps > 0. , 3] / (
                1 + self.dt * (8. * self.constants.r_recovery * ps[ps > 0.] / 
                               coral.Csp + self.constants.r_mortality * coral.Csp)
        )
        # > pale pop.
        p[ps > 0. , 2] = (coral.p0[ps > 0. , 2] + (
                8. * self.dt * self.constants.r_recovery * ps[ps > 0. ] / coral.Csp
        ) * p[ps > 0. , 3]) / (1. + self.dt * self.constants.r_recovery * ps[ps > 0.] * coral.Csp)
        # > recovering pop.
        p[ps > 0. , 1] = (coral.p0[ps > 0. , 1] +
                        self.dt * self.constants.r_recovery * ps[ps > 0. ] * coral.Csp * p[ps > 0. , 2]) / (
                1. + .5 * self.dt * self.constants.r_recovery * ps[ps > 0. ] * coral.Csp
        )
        # > healthy pop.
        a = self.dt * self.constants.r_growth * ps[ps > 0.] * coral.Csp / coral.cover[ps > 0.]
        b = 1. - self.dt * self.constants.r_growth * ps[ps > 0.] * coral.Csp * (
                1. - p[ps > 0., 1:].sum(axis=1) / coral.cover[ps > 0.]
        )
        c = - (coral.p0[ps > 0., 0] + .5 * self.dt * self.constants.r_recovery *
               ps[ps > 0.] * coral.Csp * p[ps > 0., 1])
        p[ps > 0., 0] = (-b + np.sqrt(b ** 2 - 4. * a * c)) / (2. * a)

        # bleaching conditions
        # > healthy pop.
        p[ps <= 0., 0] = coral.p0[ps <= 0., 0] / (1. - self.dt * self.constants.r_bleaching * ps[ps <= 0.] * coral.Csp)
        # > recovering pop.
        p[ps <= 0., 1] = coral.p0[ps <= 0., 1] / (1. - self.dt * self.constants.r_bleaching * ps[ps <= 0.] * coral.Csp)
        # > pale pop.
        p[ps <= 0., 2] = (coral.p0[ps <= 0., 2] - self.dt * self.constants.r_bleaching * ps[ps <= 0.] * coral.Csp * (
                p[ps <= 0., 0] + p[ps <= 0., 1]
        )) / (1. - .5 * self.dt * self.constants.r_bleaching * ps[ps <= 0.] * coral.Csp)
        # > bleached pop.
        p[ps <= 0., 3] = (coral.p0[ps <= 0., 3] -
                         .5 * self.dt * self.constants.r_bleaching * ps[ps <= 0] * coral.Csp * p[ps <= 0., 2]) / (
                1. - .25 * self.dt * self.constants.r_bleaching * ps[ps <= 0.] * coral.Csp
        )

        # # check on carrying capacity
        if any(p.sum(axis=1) > 1.0001 * coral.cover):
            slot_1 = np.arange(len(coral.cover))[p.sum(axis=1) > 1.0001 * coral.cover]
            slot_2 = p[p.sum(axis=1) > 1.0001 * coral.cover]
            slot_3 = coral.cover[p.sum(axis=1) > 1.0001 * coral.cover]
            print(
                f'WARNING: Total population than carrying capacity at {slot_1}. '
                f'\n\tPT = {slot_2}; K = {slot_3}'
            )

        # # output
        return p


class Calcification:
    """Calcification rate."""

    def __init__(self, constants):
        """Calcification rate."""
        self.ad = 1 # what is this? Aragonite dependency???? Why is it set to 1 then? 
        self.constants = constants


    def calcification_rate(self, coral, omega):
        """Calcification rate.

        :param coral: coral animal
        :param omega: aragonite saturation state

        :type coral: Coral
        :type omega: float, list, tuple, numpy.ndarray
        """
        
        def aragonite_dependency(calcification_object):
            """Aragonite dependency."""
            calcification_object.ad = (omega - self.constants.omega0) / (
                    self.constants.kappaA + omega - self.constants.omega0)
            calcification_object.ad = RESHAPE.variable2matrix(calcification_object.ad, 'time')
            
        zero_calc = np.zeros(coral.photo_rate.shape) # creating an array of zeros

        aragonite_dependency(self) #Why here we do not write return ad=  aragonite_dependency(self)
        calcification = self.constants.gC * coral.Csp * coral.pop_states[:, :, 0] * self.ad * coral.photo_rate 
        
        coral.calc = np.maximum(zero_calc, calcification) # this compares zero array and calculated with formula calcification an for each value should return the maxium one
        # So if calculated in calcification value is negative, it will become Zero.
        
        # or here coral.araginite_dependency instead of self.ad?
        # And where did P(T) and P(u) get lost? 
        
        # Shouldn't it be like this:
        # coral.calc = self.constants.gC * coral.Csp * coral.pop_states[:, :, 0] * self.ad * coral.photo_rate * coral.pld * coral.pfd
        # No,  because coral.photo_rate already has that dependencies

class Morphology:
    """Morphological development."""
    
    # Should be updated, so tehy are not constant, as I made before
    # But maybe I can set dependence for branching coral in a way that if no
    # space available in top canopy, then they grow up
    
    __rf_optimal = None
    __rp_optimal = None
    __rs_optimal = None

    def __init__(self, constants, calc_sum, light_in, dt_year=1):
        """
        Morphological development.

        :param calc_sum: accumulation of calcification of :param dt_year: years [kg m-2 yr-1]
        :param light_in: incoming light-intensity at water-air interface [umol photons m-2 s-1]
        :param dt_year: update interval [yr], defaults to 1

        :type calc_sum: float, int, list, tuple, numpy.ndarray
        :type light_in: float, int, list, tuple, numpy.ndarray
        :type dt_year: float, int
        """
        try:
            _ = len(calc_sum[0])
        except TypeError:
            self.calc_sum = calc_sum
        else:
            self.calc_sum = RESHAPE.matrix2array(calc_sum, 'space', 'sum')
        self.dt_year = dt_year

        self.I0 = RESHAPE.variable2matrix(light_in, 'time')
        self.vol_increase = 0 # does it mean that it is 0 at the first moment, while we initialise the model?
        
        self.constants = constants

    @staticmethod
    def __coral_object_checker(coral):
        """Check the suitability of the coral-object for the morphological development.

        :param coral: coral animal
        :type coral: Coral
        """
        # coral must be of type Coral
        if not isinstance(coral, Coral):
            msg = f'The optimal ratios are set using the Coral-object, {type(coral)} is given.'
            raise TypeError(msg)

        # coral must have light and flow condition attributes
        if not hasattr(coral, 'light') and not hasattr(coral, 'ucm'):
            msg = f'The optimal ratios are determined based on the coral\'s light and flow conditions; ' \
                  f'none are provided.'
            raise AttributeError(msg)

    @property
    def rf_optimal(self):
        """Optimal form ratio; height-to-plate diameter.

        :rtype: float, numpy.ndarray
        """
        return self.__rf_optimal

    @rf_optimal.setter
    def rf_optimal(self, coral):
        """
        :param coral: coral animal
        :type coral: Coral
        """
        self.__coral_object_checker(coral)
        
        rf = self.constants.rf
        
        # rf = self.constants.prop_form * (coral.light.mean(axis=1) / 
        #                                  self.I0.mean(axis=1)) * (self.constants.u0 / 1e-6)
        # rf[coral.ucm > 0] = self.constants.prop_form * (
        #         coral.light.mean(axis=1)[coral.ucm > 0] / self.I0.mean(axis=1)[coral.ucm > 0]
        # ) * (self.constants.u0 / coral.ucm[coral.ucm > 0])
        self.__rf_optimal = rf

    @property
    def rp_optimal(self):
        """Optimal plate ratio; base diameter-to-plate diameter.

        :rtype: float, numpy.ndarray
        """
        return self.__rp_optimal

    @rp_optimal.setter
    def rp_optimal(self, coral):
        """
        :param coral: coral animal
        :type coral: Coral
        """
        self.__coral_object_checker(coral)
        self.__rp_optimal = self.constants.rp

        # self.__rp_optimal = self.constants.prop_plate * (
        #         1. + np.tanh(self.constants.prop_plate_flow * (coral.ucm - self.constants.u0) / self.constants.u0)
        # )

    @property
    def rs_optimal(self):
        """Optimal spacing ratio; plate diameter-to-axial distance.

        :rtype: float, numpy.ndarray
        """
        return self.__rs_optimal

    @rs_optimal.setter
    def rs_optimal(self, coral):
        """
        :param coral: coral animal
        :type coral: Coral
        """
        
        # optimal spacing should be dependant on the distance petween "stems" + 95% confidence, so they don't touch each other
        
        # should then add the vertical component to it as well? check last years formulations 
        self.__coral_object_checker(coral)
        
        # self.__rs_optimal = 0.5/np.sqrt(2.0) * 0.25

        self.__rs_optimal = self.constants.prop_space * (
            1. - np.tanh(self.constants.prop_space_light * 
                              coral.light.mean(axis=1) / self.I0.mean(axis=1))
         ) * (1. + np.tanh(self.constants.prop_space_flow * 
                           (coral.ucm - self.constants.u0) / self.constants.u0))

    def delta_volume(self, coral):
        """
        :param coral: coral object
        :type coral: Coral
        """
        # self.vol_increase = .5 * coral.ac ** 2 * self.calc_sum * self.dt_year /self.constants.rho_c * coral.Bc.mean(axis=1)
        
        # why is it multiplied to 0.5? - found on page 128 thesis
        
        # BUT: How can I limit the coral volume increase here with spacing? 
        # In this case, voulume is measured in m^5? - Does not make sense
        # volume added depends on the active SA and calcification rate and skeletal density
        
        # Plus, only healthy coral calcifyes for growth, so shouldn't I also multply by P(Healthy)? P[:,0] - Already implemented into coral.calc_sum
        # The other thing is if we cannot add it to the periphery and should grow in height. 
        
        self.vol_increase = self.calc_sum * self.dt_year /self.constants.rho_c * coral.surface_area    #.mean(axis=1) 
        # I remove the mean operator, as it makes the model bug. In the old version, biomass was calculated daily, 
        # so had dimensions (313, 366) with the exact same value each day, BECAUSE tc, dc, hc, etc were given annualy, so had 
        # dimension hc (313, ).
        # So, it does not make sense to have daily calculation and receive exact same value every day
        
        # Either, it has to be daily calc multiplied by surface area of that exact day
        # Or multiplication by biomass survived by the end of the year (that is still the same, as other dc, hc are annual)
        
        # So, I remove "mean", as my surface_area is already (313, ), so of course there is nothing to take mean of. 
    
        return self.vol_increase

    def ratio_update(self, coral, ratio):
        """
        :param coral: coral object
        :param ratio: morphological ratio to update

        :type coral: Coral
        :type ratio: str
        """

        # partial differential equation - mass balance
        def mass_balance(r_old, r_opt):
            """Mass balance."""
            return (coral.volume * r_old + self.vol_increase * r_opt) /(coral.volume + self.vol_increase)

        # input check
        ratios = ('rf', 'rp', 'rs')
        if ratio not in ratios:
            msg = f'{ratio} not in {ratios}.'
            raise ValueError(msg)

        # calculations
        self.delta_volume(coral)

        # optimal ratio
        setattr(self, f'{ratio}_optimal', coral)

        # update morphological ratio
        if hasattr(self, f'{ratio}_optimal') and hasattr(coral, ratio):
            return mass_balance(getattr(coral, ratio), getattr(self, f'{ratio}_optimal'))

    def update(self, coral):
        """Update morphology.

        :param coral: coral animal
        :type coral: Coral
        """
        # # calculations
        # updated ratios
        ratios = [self.ratio_update(coral, ratio) for ratio in ('rf', 'rp', 'rs')]
        
        # updated volume
        volume = coral.volume + self.vol_increase
        
        # Adding update for planar area here? or maybe
        

        # update coral morphology
        coral.update_morphology(volume, *ratios)
        
        


class Dislodgement:
    """Dislodgement due to storm conditions."""
    # check the Porous media theory for dislodgement update

    def __init__(self, constants):
        """Dislodgement check."""
        self.dmt = None
        self.csf = None
        self.survival = None
        
        self.constants = constants


    def update(self, coral, survival_coefficient=1):
        """Update morphology due to storm damage.

        :param coral: coral animal
        :param survival_coefficient: percentage of partial survival, defualts to 1 # - So all of them survive? 

        :type coral: Coral
        :type survival_coefficient: float, optional
        """
        # # partial dislodgement
        Dislodgement.partial_dislodgement(self, coral, survival_coefficient)
        # # update
        # ulation states
        for s in range(4):
            coral.p0[:, s] *= self.survival
        # morphology
        coral.volume *= self.survival

    def partial_dislodgement(self, coral, survival_coefficient=1.):
        """Percentage surviving storm event.

        :param coral: coral animal
        :param survival_coefficient: percentage of partial survival, defualts to 1

        :type coral: Coral
        :type survival_coefficient: float, optional
        """
        # TODO: Rewrite such that the distinction between an array or a float is well build in.
        try:
            self.survival = np.ones(coral.dc.shape)
        except TypeError:
            if Dislodgement.dislodgement_criterion(self, coral):
                self.survival = survival_coefficient * self.dmt / self.csf
            else:
                self.survival = 1.
        else:
            dislodged = Dislodgement.dislodgement_criterion(self, coral)
            self.survival[dislodged] = survival_coefficient * (
                    self.dmt[dislodged] / self.csf[dislodged]
            )

    def dislodgement_criterion(self, coral):
        """Dislodgement criterion. Returns boolean (array).

        :param coral: coral animal
        :type coral: Coral
        """
        self.dislodgement_mechanical_threshold(coral)
        self.colony_shape_factor(coral)
        return self.dmt <= self.csf

    def dislodgement_mechanical_threshold(self, coral):
        """Dislodgement Mechanical Threshold.

        :param coral: coral animal
        :type coral: Coral
        """
        # # check input
        if not hasattr(coral.um, '__iter__'):
            coral.um = np.array([coral.um])
        if isinstance(coral.um, (list, tuple)):
            coral.um = np.array(coral.um)

        # # calculations
        self.dmt = 1e20 * np.ones(coral.um.shape)
        self.dmt[coral.um > 0] = self.dmt_formula(self.constants, coral.um[coral.um > 0])

    @staticmethod
    def dmt_formula(constants, flow_velocity):
        """Determine the Dislodgement Mechanical Threshold.

        :param flow_velocity: depth-averaged flow velocity
        :type flow_velocity: float, numpy.ndarray
        """
        return constants.sigma_t / (constants.rho_w * constants.Cd * flow_velocity ** 2)

    def colony_shape_factor(self, coral):
        """Colony Shape Factor.

        :param coral: coral animal
        :type coral: Coral
        """
        self.csf = CoralOnly().in_space(
            coral=coral,
            function=self.csf_formula,
            args=(coral.dc, coral.hc, coral.bc, coral.tc)
        )

    @staticmethod
    def csf_formula(dc, hc, bc, tc):
        """Determine the Colony Shape Factor.

        :param dc: diameter coral plate [m]
        :param hc: coral height [m]
        :param bc: diameter coral base [m]
        :param tc: thickness coral plate [m]

        :type dc: float, numpy.ndarray
        :type hc: float, numpy.ndarray
        :type bc: float, numpy.ndarray
        :type tc: float, numpy.ndarray

        :return: colony shape factor
        :rtype: float, numpy.ndarray
        """
        # arms of moment
        arm_top = hc - .5 * tc
        arm_bottom = .5 * (hc - tc)
        # area of moment
        area_top = dc * tc
        area_bottom = bc * (hc - tc)
        # integral
        integral = arm_top * area_top + arm_bottom * area_bottom
        # colony shape factor
        return 16. / (np.pi * bc ** 3) * integral


class Recruitment:
    """Recruitment dynamics."""
    def __init__(self, constants):
        """Recruitment initialize"""
        self.constants = constants

    def update(self, coral):
        """Update coral cover / volume after spawning event.

        :param coral: coral animal
        :type coral: Coral
        """
        coral.p0[:, 0] += Recruitment.spawning(self, coral, 'P')
        coral.volume += Recruitment.spawning(self, coral, 'V')

    def spawning(self, coral, param):
        """Contribution due to mass coral spawning.

        :param coral: coral animal
        :param param: parameter type to which the spawning is added

        :type coral: Coral
        :type param: str
        """
        # # input check
        params = ('P', 'V')
        if param not in params:
            msg = f'{param} not in {params}.'
            raise ValueError(msg)

        # # calculations
        # potential
        power = 2 if param == 'P' else 3
        potential = self.constants.prob_settle * self.constants.no_larvae * self.constants.d_larvae ** power
        # recruitment
        averaged_healthy_pop = coral.pop_states[:, -1, 0].mean()
        # living cover
        living_cover = RESHAPE.matrix2array(coral.living_cover, 'space')

        recruited = CoralOnly().in_space(
            coral=coral,
            function=self.recruited,
            args=(potential, averaged_healthy_pop, living_cover, coral.cover)
        )

        # # output
        return recruited

    @staticmethod
    def recruited(potential, averaged_healthy_pop, cover_real, cover_potential):
        """Determination of recruitment.

        :param potential: recruitment potential
        :param averaged_healthy_pop: model domain averaged healthy population
        :param cover_real: real coral cover
        :param cover_potential: potential coral cover

        :type potential: float
        :type averaged_healthy_pop: float
        :type cover_real: float
        :type cover_potential: float
        """
        return potential * averaged_healthy_pop * (1 - cover_real / cover_potential)
    