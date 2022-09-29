"""
coral_model - core

@author: Gijs G. Hendrickx
@contributor: Peter M.J. Herman
"""

import numpy as np
import pandas as pd
import distutils.util as du
from scipy.optimize import newton

from coral_model.utils import DataReshape, CoralOnly

# # data formatting -- to be reformatted in the model simulation
RESHAPE = DataReshape()
     
# coral object
class Coral:
    """Coral object, representing one coral type."""

    def __init__(self, constants, dc, hc, bc, tc, ac, pop_den, species_constant):
        """
        :param dc: diameter coral plate [m]
        :param hc: coral height [m]
        :param bc: diameter coral base [m]
        :param tc: thickness coral plate [m]
        :param ac: axial distance corals [m]
        :param pop_den: population density of corals [m-2]
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
        
        self.dc = RESHAPE.variable2array(dc)
        self.hc = RESHAPE.variable2array(hc)
        self.bc = RESHAPE.variable2array(bc)
        self.tc = RESHAPE.variable2array(tc)
        self.ac = RESHAPE.variable2array(ac)
        
        self.constants = constants #needed to return them, because they are used in photosynthesis
        
        self.pop_den = RESHAPE.variable2array(pop_den)

        self.Csp = species_constant #make into the constants list

        # initiate environmental working objects
        # > light micro-environment
        self.light = None
        # > flow micro-environment
        self.um = None
        # > thermal environment
        self.temp = None
        # > photosynthesis
        self.photo_rate = None
        self.Tlo = None
        self.Thi = None
        # > population states
        self.pop_states = None
        self.p0 = None
        # > calcification
        self.calc = None
        
        """
           
        # Key = value           ! Default   ! Definition [units]
        #--------------------------------------------------------------------------------------------------------------------
    
        # light attenuation
        coral_Kd_1 = 0.         #  0.0      ! constant within-colony light-attenuation coefficient, lit side (0.0 for massive, add value for branching)
        coral_Kd_2 = 0.         #  0.0      ! constant within-colony light-attenuation coefficient, strongly shaded side (0.0 for massive, add value for branching)
        
        # photosynthetic light dependency
        iota = .6               # .6        ! photo-acclimation rate [d-1]
        ik_max = 372.32         # 372.32    ! maximum quasi steady-state saturation light-intensity [umol photons m-2 s-1]
        pm_max = 1.             #  1.       ! maximum quasi steady-state maximum photosynthetic efficiency [-]
        rETR_max_qss = 70.0     # 70.0      ! maximum quasi steady-state maximum relative electron transport rate [-]
        beta_rETR = 0.09        # 0.09      ! exponent of the quasi steady-state maximum relative electron transport rate [-]  
        betaI = .34             # .34       ! exponent of the quasi steady-state saturation light-intensity [-]
        betaP = .09             # .09       ! exponent of the quasi steady-state maximum photosynthetic efficiency [-]
        Icomp = .01             # .01       ! fraction of I0 at compensation point [-]
        phi = 0.05              # 0.05      ! Photoinhibition coefficient, defines the slope of photosynthetic efficiency reduction [-]
        alpha = 0.2             # 0.2       ! Initial slope of photosynthetic curve until reaching the compensation point [-]
        
        # population dynamics
        r_growth = .002         # 0.002     ! growth rate [d-1]
        r_recovery = .2         # .2        ! recovery rate [d-1]
        r_mortality = .04       # .04       ! mortality rate [d-1]
        r_bleaching =  8.       # 8.        ! bleaching rate [d-1]
        
        # morphological development
        
        dc_pa_coef =            #            ! Coefficient for diameter - planar area relationship [-]
        dc_pa_exp =             #            ! Exponent for diameter - planar area relationship [-]

        pa_vol_coef=            #            ! Coefficient for planar area - volume relationship [-]
        pa_vol_exp=             #             ! Exponent for the planar area - volume relationship [-]
        
        sa_dc_coef =            #           ! Coefficient for surface area - diameter relationship
        sa_dc_exp =             #           ! Exponent for surface area - diameter relationship
        
        rf = 1.0                # 0.8       ! form ratio height to diameter [-]
        rp = 1.0                # 1.0       ! plate ratio base to top diameter [-] 

        prop_form = .1          # .1        ! overall form proportionality constant [-]
        prop_plate = .5         # .5        ! overall plate proportionality constant [-
        prop_plate_flow = .1    # .1        !  flow plate proportionality constant [-]
        prop_space = .5         # .5/np.sqrt(2.) ! overall space proportionality constant [-]
        prop_space_light = .1   # .1      ! light space proportionality constant [-]
        prop_space_flow = .1    # .1        ! flow space proportionality constant [-]
        u0 = .2                 # .2        ! base-line flow velocity [m s-1]
        rho_c = 1600.           # 1600.     ! density of coral [kg m-3]
        #
        # coral recruitment
        no_larvae = 1e6         # 1e6       ! number of larvae released during mass spawning event [-]
        prob_settle = 1e-4      # 1e-4      ! probability of settlement [-]
        d_larvae = 1e-3         # 1e-3      ! larval diameter [m]

        """        
        # light micro-environment
        self.coral_Kd_1 = None     
        self.coral_Kd_2 = None     

        # photosynthetic light dependency
        self.iota =  None
        self.ik_max =  None
        self.pm_max =  None
        self.rETR_max_qss = None
        self.beta_rETR = None
        self.betaI = None
        self.betaP =  None
        self.Icomp = None
        self.phi = None
        self.alpha = None

        # population dynamics
        self.r_growth =  None
        self.r_recovery =  None
        self.r_mortality = None
        self.r_bleaching = None

        # morphological development
        self.dc_pa_coef = None
        self.dc_pa_exp = None

        self.pa_vol_coef= None
        self.pa_vol_exp= None
        
        self.sa_dc_coef = None    
        self.sa_dc_exp = None
        
        
        self.rf = None
        self.rp = None
        self.prop_form =  None
        self.prop_plate =  None
        self.prop_plate_flow =  None
        self.prop_space =  None
        self.prop_space_light = None
        self.prop_space_flow =  None
        self.u0 = None
        self.rho_c = None

        # coral recruitment
        self.no_larvae =  None
        self.prob_settle = None
        self.d_larvae =  None
    
    def read_coral(self,inp_file):        
        self.inpfile=inp_file
            
        keyvals={}
        with open(self.inpfile) as f:
            for line in f:
                if(len(line)>1):
                    linee = line
                    if (line.count("#")>0):
                        linee,dum = line.split ("#")
                    if(len(linee)>0):
                        name, value = linee.split("=")
                        value=value.lower().strip()
                        try:
                            keyvals[name.strip()] = float(value)
                        except (ValueError):
                            keyvals[name.strip()]=bool(du.strtobool(value))

        def default(x, default_value):
            """Set default value if no custom value is provided."""
            xx = keyvals.get(x)
            if (xx is None):
                xx = default_value
            return xx

        self.coral_Kd_1 = default("coral_Kd_1", 0.)    
        self.coral_Kd_2 = default("coral_Kd_2", 0.)      

        # photosynthetic light dependency
        self.iota = default("iota", .6)
        self.ik_max = default("ik_max", 372.32)
        self.pm_max = default("pm_max", 1.)
        self.rETR_max_qss = default("rETTR_max_qss", 70.0)
        self.beta_rETR = default("beta_rETR", 0.09)  
        self.betaI = default("betaI", .34)
        self.betaP = default("beta_P", .09)
        self.Icomp = default("Icomp", .01)
        self.phi = default("phi", 0.05)
        self.alpha = default("alpha", 0.2)

        # population dynamics
        self.r_growth = default("r_growth", .002)
        self.r_recovery = default("r_recovery", .2)
        self.r_mortality = default("r_mortality", .04)
        self.r_bleaching = default("r_bleaching", 8.)

        # morphological development
        
        self.dc_pa_coef = default ("dc_pa_coef", 0.5557)
        self.dc_pa_exp = default ("dc_pa_exp", 1.9458)

        self.pa_vol_coef= default ("pa_vol_coef", 0.2589)
        self.pa_vol_exp= default ("pa_vol_exp", 1.4482)
        
        self.sa_dc_coef = default ("sa_dc_coef", 1.6785) 
        self.sa_dc_exp = default("sa_dc_exp", 1.8841)
        
        self.rf = default("rf", 1.0)
        self.rp = default("rp", 1.0)
        self.prop_form = default("prop_form", .1)
        self.prop_plate = default("prop_plate", .5)
        self.prop_plate_flow = default("prop_plate_flow", .1)
        self.prop_space = default("prop_space", .5) / np.sqrt(2.)
        self.prop_space_light = default("prop_space_light", .1)
        self.prop_space_flow = default("prop_space_flow", .1)
        self.u0 = default("u0", .2)
        self.rho_c = default("rho_c", 1600.)

        # coral recruitment
        self.no_larvae = default("no_larvae", 1e6)
        self.prob_settle = default("prob_settle", 1e-4)
        self.d_larvae = default("d_larvae", 1e-3)

    def __repr__(self):
        """Development representation."""
        return f'Morphology({self.dc}, {self.hc}, {self.bc}, {self.bc}, {self.ac})'

    def __str__(self):
        """Print representation."""
        return (f'Coral morphology with: dc = {self.dc} m; hc = {self.hc} ;'
                f'bc = {self.bc} m; tc = {self.tc} m; ac = {self.ac} m')

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
        planar_area = self.dc_pa_coef * ((self.dc)** self.dc_pa_exp) 
        return planar_area 
    
    @property
    def surface_area(self):
        """Active area of the coral. Basically, its surface area"""
        
        SA = self.sa_dc_coef  * (self.dc ** self.sa_dc_exp)
        return SA

    
    @property
    def volume(self):
        """Coral volume."""
        coral_volume = self.pa_vol_coef * ((self.planar_area)**self.pa_vol_exp)
        return coral_volume

    @volume.setter 
    def volume(self, coral_volume):
        """
        :param coral_volume: coral volume [m3]
        :type coral_volume: float, int, list, tuple, np.ndarray
        """
        self.update_morphology(coral_volume, rf=self.rf, rp=self.rp, rs=self.rs) 
        
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
        return RESHAPE.variable2matrix(self.dc, 'space') 

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
    

    def initiate_spatial_morphology(self, pop_den):
        """Initiate the morphology based on the on set of morphological dimensions and the coral cover. This method
        contains a catch that it can only be used to initiate the morphology, and cannot overwrite existing spatial
        heterogeneous morphology definitions.

        :param cover: custom coral cover definition, defaults to None
        :type cover: None, numpy.ndarray
        """

        self.p0 = np.array([
            np.ones,
            np.zeros(pop_den.shape),
            np.zeros(pop_den.shape),
            np.zeros(pop_den.shape),
            ]).transpose()
        

        self.dc = cover * self.dc #  need to distribute these variables over the whole domain somehow. But now cover is not used, so needs to be a new way
        self.hc = cover * self.hc
        self.bc = cover * self.bc
        self.tc = cover * self.tc
        self.ac = cover * self.ac 

    def Light(self, light_in, lac, depth):
        """Light micro-environment.
    
            :param light_in: incoming light-intensity at the water-air interface [u mol photons m-2 s-1]
            :param lac: light-attenuation coefficient [m-1]
            :param depth: water depth [m]
    
            :type light_in: float, list, tuple, numpy.ndarray
            :type lac: float, list, tuple, numpy.ndarray
            :type depth: float, list, tuple, numpy.ndarray
            """
        self.I0 = RESHAPE.variable2matrix(light_in, 'time')
        self.Kd = RESHAPE.variable2matrix(lac, 'time') 
        self.h = RESHAPE.variable2matrix(depth, 'space')

        def light_spreading(self):
            """Spreading of light as function of depth. """
            return self.constants.theta_max * np.exp(
                -self.Kd * (self.h - self.hc_matrix + self.tc_matrix)) #If I do not have tc, then update it here too!
          
        def rep_light_branching(self,coral):
                   
            top = self.planar_area_matrix * self.I0 * np.exp(
                -self.Kd * (self.h - self.hc_matrix))    
            
            # For the rest, according to Kaniewska, 2014 - there is within colony light attenuation
            # different for branches direction
            
            # at the same time, massive corals do not have the branches and much extra SA,
            # So within colony LAC is not needed for them 
        
            
            side1 = ((self.surface_area_matrix - self.planar_area_matrix)/2.0) * self.I0*np.exp(
                -self.Kd*(self.h - self.hc_matrix)) *np.exp(-self.coral_Kd_1*(self.tc_matrix/2.))
            
            side2 = ((self.surface_area_matrix - self.planar_area_matrix)/2.0) * self.I0*np.exp(
                -self.Kd*(self.h - self.hc_matrix)) *np.exp(-self.coral_Kd_2*(self.tc_matrix/2.))
            
            
            total_branching = top + side1 + side2
            
    
    
            def averaged_light_branching(total_light, surface_area):
                """Averaged light-intensity."""
                return total_light / surface_area
            
            # Still good and a validation point for average light intensity dispersion 
            # over the colony (Kaniewska 2014, colony averaged PAR)
            
            #TODO: check what the CoralOnly stuff is doing here
            # update the properties for the coral processes then, where calls for light? Because then will have separate for branching and for massive
            #TODO: should the self.light things be different? But then it is really cofusing
            # if it is the unified property
            
            # maybe there is the way to sspecify it just per Class Massive and Branching? 
            
            # like have common self.light
            # then have as a property specification for class Massive self.light = self.photosynthesis.massive_light() ????
                  
            coral.branching_light = CoralOnly().in_spacetime(
                coral=coral,
                function=averaged_light_branching,
                args=(total_branching, coral.surface_area),
                no_cover_value=self.I0 * np.exp(-self.Kd * self.h)) 
            
        def rep_light_massive():
            
            def averaged_light_massive():
                evnedf vdefnvkfvio
                return alvndfvl

                
    def Temperature(self, temperature):
        """
            Thermal micro-environment.

            Parameters
        ----------
        temperature : numeric
            Temperature of water [K].
        """
        self.temp = RESHAPE.variable2matrix(temperature, 'time')
        return self.temp # do I really need to return it like this, or I can do immediately: return RESHAPE.variable2matrix(temperature, 'time')


    def photosynthesis(self, light_in):
        """Photosynthesis.
           
            Photosynthetic efficiency based on photosynthetic dependencies.
    
            :param light_in: incoming light-intensity at the water-air interface [umol photons m-2 s-1]
            :param first_year: first year of the model simulation
    
            :type light_in: float, list, tuple, numpy.ndarray
            :type first_year: bool
            """
        self.I0 = RESHAPE.variable2matrix(light_in, 'time')
            
        self.pld = 1 # can I introduce these self.pld here, or they have to be all moved to the beginning where we call cof coral processes?
        self.ptd = 1 # would it use the functions below, or then I have to call for photo_rate in the end of the function? 

        def photo_rate(self, environment, year):
            """Photosynthetic efficiency.
    
            :param environment: environmental conditions
            :param year: year of simulation
    
            :type environment: Environment
            :type year: int
            """
            # components
            self.light_dependency('qss') # also light dependency is described below, so should I reshuffle the functions as well?
            
            self.thermal_dependency(environment, year)
    
            # combined
            self.photo_rate = self.pld * self.ptd
    
        def light_dependency(self, output):
            """Photosynthetic light dependency.
    
            :param output: type of output
    
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
                x_max = self.ik_max if param == 'Ik' else self.rETR_max_qss
                beta_x = self.betaI if param == 'Ik' else self.beta_rETR
                
                # Why in text Table 4.1 Pm_max = 3.96, but default in the model input and environment it is 1.0?
    
                # calculations
                xs = x_max * (self.light / self.I0) ** beta_x
                if output == 'qss':
                    return xs
                elif output == 'new':
                    return xs + (x_old - xs) * np.exp(-self.iota)
                
                # but "new" is never used? - formula from Anthony, H-G(2003) for adaptation to new conditions
                
            # # parameter definitions
            if output == 'qss':
                ik = photo_acclimation(None, 'Ik')
                rETR_max_z = photo_acclimation( None, 'rETRmax')
            else:
                msg = f'Only the quasi-steady state solution is currently implemented; use key-word \'qss\'.'
                raise NotImplementedError(msg)
    
    
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
                
            a = (self.alpha * self.light)/rETR_max_z
                
            b = (self.phi * self.light)/rETR_max_z   
            
            self.pld = rETR_max_z  * (1.0 - np.exp(-a)) * np.exp(-b)
            
    
        def thermal_dependency(self, env, year):
            """Photosynthetic thermal dependency.
    
            :param env: environmental conditions
            :param year: year of simulation
   
            :type env: Environment
            :type year: int
            """
    
            def thermal_acc():
                """Thermal-acclimation.""" # can I just remove everything from this if-loop, cause I do not have the tme anymore?
                # if self.constants.tme:
                #     if self.first_year:
                #         env.tmeMMMmin = pd.DataFrame(
                #             data=pd.concat(
                #                 [env.temp_mmm['min']] * RESHAPE.space,
                #                 axis=1
                #             ).values,
                #             columns=[np.arange(RESHAPE.space)]
                #         ) 
                #         env.tmeMMMmax = pd.DataFrame(
                #             data=pd.concat(
                #                 [env.temp_mmm['max']] * RESHAPE.space,
                #                 axis=1
                #             ).values,
                #             columns=[np.arange(RESHAPE.space)]
                #         ) 
                #     else:
                #         env.tmeMMMmin[env.tmeMMM.index == year] += coral.dTc
                #         env.tmeMMMmax[env.tmeMMm.index == year] += coral.dTc
    
                #     mmm_min = env.tmeMMMmin[np.logical_and(
                #         env.tmeMMM.index < year,
                #         env.tmeMMM.index >= year - int(self.constants.nn / self.Csp)
                #     )]
                #     m_min = mmm_min.mean(axis=0)
                #     s_min = mmm_min.std(axis=0)
    
                #     mmm_max = env.tmeMMMmax[np.logical_and(
                #         env.tmeMMM.index < year,
                #         env.tmeMMM.index >= year - int(self.constants.nn / self.Csp)
                #     )]
                #     m_max = mmm_max.mean(axis=0)
                #     s_max = mmm_max.std(axis=0)
                # else:
                mmm = env.temp_mmm[np.logical_and(
                        env.temp_mmm.index < year,
                        env.temp_mmm.index >= year - int(self.constants.nn / self.Csp)
                    )]
                m_min, m_max = mmm.mean(axis=0)
                s_min, s_max = mmm.std(axis=0)
    
                self.Tlo = m_min - self.constants.k_var * s_min
                self.Thi = m_max + self.constants.k_var * s_max # same, do I need to explicitly Return, or it knows?
    
            def adapted_temp():
                """Adapted temperature response."""
    
                def spec():
                    """Specialisation term."""
                    return 4e-4 * np.exp(-.33 * (delta_temp - 10))
    
                response = -(self.temp - self.Tlo) * ((self.temp - self.Tlo) ** 2 - delta_temp ** 2)
                temp_cr = coral.Tlo - (1 / np.sqrt(3)) * delta_temp
                try:
                    # if self.constants.tme:
                    #     response[self.temp <= temp_cr] = -(
                    #             (2 / (3 * np.sqrt(3))) * delta_temp[self.temp <= temp_cr] ** 3
                    #     )
                    # else:
                    response[coral.temp <= temp_cr] = -(
                                (2 / (3 * np.sqrt(3))) * delta_temp ** 3
                        )
                except TypeError:
                    if self.temp <= temp_cr:
                        response = (2 / (3 * np.sqrt(3))) * delta_temp ** 3
    
                return response * spec()
    
            def thermal_env():
                """Thermal envelope."""
                return np.exp((self.constants.Ea / self.constants.R) * (1 / 300 - 1 / temp_opt))
    
            # # parameter definitions
            thermal_acc()
            delta_temp = self.Thi - self.Tlo
            temp_opt = self.Tlo + (1 / np.sqrt(3)) * delta_temp
    
            # # calculations
            f1 = adapted_temp()
            f2 = thermal_env()
            self.ptd = f1 * f2

    def PopulationStates(self, dt=1):
        """Bleaching response following the population dynamics."""  # check the competition for space here
        # TODO: Check this class; incl. writing tests
    
        
            """Population dynamics.
    
            :param dt: time step [yrs], defaults to one
            :type dt: float, optional
            """
        self.dt = dt
        
        def pop_states_t(self):
            """Population dynamics over time.
    
            :param coral: coral animal
            :type coral: Coral
            """
            self.pop_states = np.zeros((*RESHAPE.spacetime, 4))
            for n in range(RESHAPE.time):
                photosynthesis = np.zeros(RESHAPE.space)
                photosynthesis[coral.cover > 0.] = self.photo_rate[coral.cover > 0., n]  # TODO": resolve cover here!
                self.pop_states[:, n, :] = self.pop_states_xy(coral, photosynthesis)
                self.p0[coral.cover > 0., :] = self.pop_states[coral.cover > 0., n, :]
    
        def pop_states_xy(self, ps):
            """Population dynamics over space.

            :param ps: photosynthetic rate

            :type ps: numpy.ndarray
            """
            
            # TODO: check where coral.cover is used and think how to remove it from the equation
            p = np.zeros((RESHAPE.space, 4))
            # # calculations
            # growing conditions
            # > bleached pop.      # ps>0. here represents ps>tsh that is the value of the bleaching treshold light and 1. where 1.0 is a number, not column reference
            p[ps > 0. , 3] = self.p0[ps > 0. , 3] / (
                    1 + self.dt * (8. * self.r_recovery * ps[ps > 0.] / 
                                   self.Csp + self.r_mortality * self.Csp)
            )
            # > pale pop.
            p[ps > 0. , 2] = (self.p0[ps > 0. , 2] + (
                    8. * self.dt * self.r_recovery * ps[ps > 0. ] / self.Csp
            ) * p[ps > 0. , 3]) / (1. + self.dt * self.r_recovery * ps[ps > 0.] * self.Csp)
            # > recovering pop.
            p[ps > 0. , 1] = (self.p0[ps > 0. , 1] +
                            self.dt * self.r_recovery * ps[ps > 0. ] * self.Csp * p[ps > 0. , 2]) / (
                    1. + .5 * self.dt * self.r_recovery * ps[ps > 0. ] * self.Csp
            )
            # > healthy pop.
            a = self.dt * self.r_growth * ps[ps > 0.] * self.Csp / coral.cover[ps > 0.]   #COVERR!
            b = 1. - self.dt * self.r_growth * ps[ps > 0.] * self.Csp * (
                    1. - p[ps > 0., 1:].sum(axis=1) / coral.cover[ps > 0.]
            )
            c = - (self.p0[ps > 0., 0] + .5 * self.dt * self.r_recovery *
                   ps[ps > 0.] * self.Csp * p[ps > 0., 1])
            p[ps > 0., 0] = (-b + np.sqrt(b ** 2 - 4. * a * c)) / (2. * a)
    
            # bleaching conditions
            # > healthy pop.
            p[ps <= 0., 0] = self.p0[ps <= 0., 0] / (1. - self.dt * self.r_bleaching * ps[ps <= 0.] * self.Csp)
            # > recovering pop.
            p[ps <= 0., 1] = self.p0[ps <= 0., 1] / (1. - self.dt * self.r_bleaching * ps[ps <= 0.] * self.Csp)
            # > pale pop.
            p[ps <= 0., 2] = (self.p0[ps <= 0., 2] - self.dt * self.r_bleaching * ps[ps <= 0.] * self.Csp * (
                    p[ps <= 0., 0] + p[ps <= 0., 1]
            )) / (1. - .5 * self.dt * self.r_bleaching * ps[ps <= 0.] * self.Csp)
            # > bleached pop.
            p[ps <= 0., 3] = (self.p0[ps <= 0., 3] -
                             .5 * self.dt * self.r_bleaching * ps[ps <= 0] * self.Csp * p[ps <= 0., 2]) / (
                    1. - .25 * self.dt * self.r_bleaching * ps[ps <= 0.] * self.Csp )
                
            return p

       
    def Calcification(self, omega):
        """Calcification rate.
        :param omega: aragonite saturation state
        :type omega: float, list, tuple, numpy.ndarray
            """
        self.ad = 1 
        def aragonite_dependency(calcification_object): # why do I need here this thing "calcification object"? 
            """Aragonite dependency."""
            calcification_object.ad = (omega - self.constants.omega0) / (
                        self.constants.kappaA + omega - self.constants.omega0)
            calcification_object.ad = RESHAPE.variable2matrix(calcification_object.ad, 'time')
                
        zero_calc = np.zeros(self.photo_rate.shape) # creating an array of zeros
    
        aragonite_dependency(self) #Why here we do not write return ad=  aragonite_dependency(self)
        calcification = self.constants.gC * self.Csp * self.pop_states[:, :, 0] * self.ad * self.photo_rate 
            
        self.calc = np.maximum(zero_calc, calcification) # this compares zero array and calculated with formula calcification an for each value should return the maxium one
            # So if calculated in calcification value is negative, it will become Zero.
            
            # or here coral.araginite_dependency instead of self.ad?
            # And where did P(T) and P(u) get lost? 
            
            # Shouldn't it be like this:
            # coral.calc = self.constants.gC * coral.Csp * coral.pop_states[:, :, 0] * self.ad * coral.photo_rate * coral.pld * coral.pfd
            # No,  because coral.photo_rate already has that dependencies
        return self.calc
    
    
    def Morphology(self, calc_sum, light_in, dt_year=1):
        """Morphological development.
   
            :param calc_sum: accumulation of calcification of :param dt_year: years [kg m-2 yr-1]
            :param light_in: incoming light-intensity at water-air interface [umol photons m-2 s-1]
            :param dt_year: update interval [yr], defaults to 1
    
            :type calc_sum: float, int, list, tuple, numpy.ndarray
            :type light_in: float, int, list, tuple, numpy.ndarray
            :type dt_year: float, int
            """
        try:
                _ = len(calc_sum[0]) # what does it do?
        except TypeError:
                self.calc_sum = calc_sum
        else:
            self.calc_sum = RESHAPE.matrix2array(calc_sum, 'space', 'sum')
            self.dt_year = dt_year
    
            self.I0 = RESHAPE.variable2matrix(light_in, 'time')
            self.vol_increase = 0 # does it mean that it is 0 at the first moment, while we initialise the model?
            
            # Why cannot I just call for self.calc and do self.calc.sum already, instead of puttion it inside here every time? 
    
        @staticmethod
        def __coral_object_checker(coral): # do I really need it? we do it for coral only, inside the Coral class. If there will be seagrass, it will have different class, for instance
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
                
        def delta_volume(self, coral):
            """
            :param coral: coral object
            :type coral: Coral
            """
            self.vol_increase = self.calc_sum * self.dt_year /self.rho_c * coral.surface_area    
        
            return self.vol_increase
    
        def ratio_update(self, coral, ratio):
            """
            :param coral: coral object
            :param ratio: morphological ratio to update
    
            :type coral: Coral
            :type ratio: str
            """
    
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
            # updated volume
            volume = coral.volume + self.vol_increase
    
            # update coral morphology
            self.update_morphology(volume, *ratios) 
            # TODO: change the function here or up in the coral to update volume and diameter calculations!   
            Can use formula for new volume - get new diameter with the constants 

    def Dislodgement(self):
        """Dislodgement due to storm conditions."""
       
            """Dislodgement check."""
        self.dmt = None
        self.csf = None
        self.survival = None
            
        # self.constants = constants # but constants are already an attribute of Coral, so I think I do not have to mention it again in the function?
    
    
        def update(self, survival_coefficient=1):
            """Update morphology due to storm damage.
    
            :param coral: coral animal
            :param survival_coefficient: percentage of partial survival, defualts to 1 # - So all of them survive? 
    
            :type coral: Coral
            :type survival_coefficient: float, optional
            """
            # # partial dislodgement
            Dislodgement.partial_dislodgement(self, coral, survival_coefficient)
            # # update
            # population states
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



class Massive(Coral):
    """ Class of Massive corals. Inherits all the initialisation properties and functions of parent super-class Coral.
    Here some class-specific functions can be added. """
    def __init__(self):
            super().__init__()
            


class Branching(Coral):
    """ Class of Branching corals. Inherits all the properties of class Coral"""
    def __init__(self):
            super().__init__()            
            



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
        self.h = RESHAPE.variable2matrix(h, 'space')
        self.Tp = RESHAPE.variable2array(peak_period)
        self.active = False if u_current is None and u_wave is None else True
        self.constants = constants

# =============================================================================
# =============================================================================
    @property
    def uc_matrix(self):
        """Reshaped current flow velocity."""
        return RESHAPE.variable2matrix(self.uc, 'space')  

    @property
    def uw_matrix(self):
        """Reshaped wave flow velocity."""
        return RESHAPE.variable2matrix(self.uw, 'space') 
           
    def velocities(self, coral):
        """Depth-averaged flow velocities.

        :param coral: coral animal
        :param in_canopy: determine in-canopy flow (or depth-averaged), defaults to True  - So is it in-canopy only, or depth averaged as well somewhow?

        :type in_canopy: bool, optional
        """
        if self.active:
            alpha_w = np.ones(self.uw.shape)
            alpha_c = np.ones(self.uc.shape)
            if in_canopy:
                idx = coral.volume > 0
                for i in idx:
                    alpha_w[i] = self.wave_attenuation(
                        coral.dc_rep[i], coral.hc[i], coral.ac[i],
                        self.uw[i], self.Tp[i], self.h[i], 'wave'
                    )
                    alpha_c[i] = self.wave_attenuation(
                        coral.dc_rep[i], coral.hc[i], coral.ac[i],
                        self.uc[i], 1e3, self.h[i], 'current'
                    )
            coral.ucm = self.wave_current(alpha_w, alpha_c)
            coral.um = self.wave_current()
        else:
            coral.ucm = 9999 * np.ones(RESHAPE.space) # what happens to um then? Because I need um now only
            
            # Should I make it always active?
            
            # then I have: 
            
         # alpha_w = np.ones(self.uw.shape)
         # alpha_c = np.ones(self.uc.shape)
            
         # coral.um = self.wave_current()
         
         # But also, if in-canopy is disabled, I do nod need the wave attenuation
         # function that depends on coral morphology and calculated drag coefficient,
         # which is strange not to have for dislodgement calculations
         
         # in case we use only the depth-averaged flow, which is um (needed for dislodgement)
         # it is strange that it is not specified in ELSE: in case in_canopy = False
         
         # So what does it calculate then?         
         

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
    def wave_attenuation(constants, diameter, height, distance, velocity, period, depth, wac_type): #why doesn't it just depend on coral? that has height, diam, etc...
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
                porous_flow = wac * above_flow # I do not have it, if there is no in-canopy flow, right?
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






class Recruitment:
    """Recruitment dynamics."""
    def __init__(self):
        """Recruitment initialize"""
        
        # addition of volume of 1 sexually-reproducable coral (5 cm diameter) 

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
        potential = self.prob_settle * self.no_larvae * self.d_larvae ** power
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
    