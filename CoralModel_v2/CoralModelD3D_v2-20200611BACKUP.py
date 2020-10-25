# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:32:23 2020

@author: hendrick
"""

import numpy as np
import pandas as pd
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import bmi.wrapper
import os
from scipy.optimize import newton
from tqdm import tqdm
import datetime
from netCDF4 import Dataset
import faulthandler
faulthandler.enable()

# =============================================================================
# # # # specify directories of ddl- and input-files
# =============================================================================
model_folder = os.path.join('MiniModel2')

# Delft3D directories
D3D_HOME = os.path.join('p:\\11202744-008-vegetation-modelling', 'code_1709',
                        'windows', 'oss_artifacts_x64_63721', 'x64')
dflow_dir = os.path.join(D3D_HOME, 'dflowfm', 'bin', 'dflowfm.dll')
dimr_path = os.path.join(D3D_HOME, 'dimr', 'bin', 'dimr_dll.dll')
# work directory
workdir = os.path.join('p:\\11202744-008-vegetation-modelling', 'students',
                       'GijsHendrickx', 'models', model_folder)
inputdir = os.path.join(workdir, 'timeseries')
# input files (Delft3D)
config_file = os.path.join(workdir, 'dimr_config.xml')
mdufile = os.path.join(workdir, 'fm', 'FlowFM.mdu')
# print directories and input-files as check
print('Model              : {0}\n'.format(workdir))
print('Delft3D home       : {0}'.format(D3D_HOME))
print('DIMR-directory     : {0}'.format(dimr_path))
print('Configuration file : {0}'.format(config_file))

# =============================================================================
# # # # prepare locations
# =============================================================================
# # print directories of input- and output-files
print('\nTime-series dir.   : {0}'.format(inputdir))

# # intermediate figures
figfolder = os.path.join(workdir, 'figures')
# check existence and create if necessary
if not os.path.exists(figfolder):
    os.mkdir(figfolder)
    print('New folder created : {0}'.format(figfolder))
print('Figure directory   : {0}'.format(figfolder))

# # output files
outputfolder = os.path.join(workdir, 'output')
# check existance and create if necessary
if not os.path.exists(outputfolder):
    os.mkdir(outputfolder)
    print('New folder created : {0}'.format(outputfolder))
print('Output directory   : {0}'.format(outputfolder))

# =============================================================================
# # # # environmental conditions - files
# =============================================================================
fLight = 'TS_PAR.txt'
fLAC = 'TS_LAC.txt'
fTemperature = 'TS_SST.txt'
fAcidity = 'TS_ARG.txt'
fStorm = 'TS_stormcat.txt'

# =============================================================================
# # # # create correct environment
# =============================================================================
os.environ['PATH'] = (os.path.join(D3D_HOME, 'share', 'bin') + ';' +
                      os.path.join(D3D_HOME, 'dflowfm', 'bin') + ';' +
                      os.path.join(D3D_HOME, 'dimr', 'bin') + ';' +
                      os.path.join(D3D_HOME, 'dwaves', 'bin') + ';' +
                      os.path.join(D3D_HOME, 'esmf', 'scripts') + ';' +
                      os.path.join(D3D_HOME, 'swan', 'scripts'))
# print created environment as check
print('\nEnvironment        : {0}\n'
      .format(os.path.join(D3D_HOME, 'share', 'bin')) +
      '                     {0}\n'
      .format(os.path.join(D3D_HOME, 'dflowfm', 'bin')) +
      '                     {0}\n'
      .format(os.path.join(D3D_HOME, 'dimr', 'bin')) +
      '                     {0}\n'
      .format(os.path.join(D3D_HOME, 'dwaves', 'bin')) +
      '                     {0}\n'
      .format(os.path.join(D3D_HOME, 'esmf', 'scripts')) +
      '                     {0}\n'
      .format(os.path.join(D3D_HOME, 'swan', 'scripts')))

# =============================================================================
# # # # define and initialize wrappers
# =============================================================================
# define DFM wrapper
modelFM = bmi.wrapper.BMIWrapper(engine=dflow_dir, configfile=mdufile)
# define DIMR wrapper
modelDIMR = bmi.wrapper.BMIWrapper(engine=dimr_path, configfile=config_file)
# initialise model
modelDIMR.initialize()

print('Model initialized.\n')

# =============================================================================
# # # # set the pointers to important model variables of FlowFM
# =============================================================================
# number of boxes, including boundary boxes
ndx = modelFM.get_var('ndx')
# number of non-boundary boxes, i.e. within-domain boxes
ndxi = modelFM.get_var('ndxi')
# x-coord. of the center of gravity of the boxes
xzw = modelFM.get_var('xzw')
# y-coord. of the center of gravity of the boxes
yzw = modelFM.get_var('yzw')
# total number of links between boxes
lnx = modelFM.get_var('lnx')
# number of links between within-domain boxes
lnxi = modelFM.get_var('lnxi')
# link martrix between adjacent boxes [ln, 2] matrix
ln = modelFM.get_var('ln')
# distance between the centers of adjacent boxes
dx = modelFM.get_var('dx')
# width of the interface between adjacent boxes
wu = modelFM.get_var('wu')
# surface area of the boxes
ba = modelFM.get_var('ba')


# =============================================================================
# # # # set time parameters for coupled model
# =============================================================================
# # time-span
# start year
Ystart = 2000
# simulation time [yrs]
Y = 100
# year range
years = np.arange(Ystart, Ystart + Y)

# # model time per vegetation step [s]
mtpervt = 43200
# storm
mtpervt_storm = 86400

# =============================================================================
# # # # define output
# =============================================================================
# # # map > full spatial extent
# # data to output file
# mean flow > { uc }
U2mfile = True
# mean coral temperature > { Tc, Tlo, Thi }
T2mfile = True
# mean photosynthesis > { PS }
PS2mfile = True
# population states > { P } > { PH, PR, PP, PB }
P2mfile = True
# calcification > { G }
G2mfile = True
# morphology > { Lc } > { dc, hc, bc, tc, ac }
M2mfile = True

# # map-file
# map-file directory
mapfile = 'CoralModel_map.nc'
mapfilef = os.path.join(outputfolder, mapfile)
# time-interval > annually

# # # history > time-series
# # data to output file
# flow > { uc }
U2hfile = True
# temperature > { Tc, Tlo, Thi }
T2hfile = True
# photosynthesis > { PS }
PS2hfile = True
# population states > { P } > { PH, PR, PP, PB }
P2hfile = True
# calcification > { G }
G2hfile = True
# morphology > { Lc } > { dc, hc, bc, tc, ac }
M2hfile = True

# # his-file
# location(s)
xynfilef = os.path.join(workdir, 'fm', 'FlowFm_obs.xyn')
xyn = pd.read_csv(xynfilef, header=None, delim_whitespace=True)
xyn.columns = ['x', 'y', 'name']
# his-file directory
hisfile = 'CoralModel_his.nc'
hisfilef = os.path.join(outputfolder, hisfile)
# time-interval > daily

# =============================================================================
# # # # vegetation boundaries
# =============================================================================
xbndmin = min(xzw[range(ndxi)])
xbndmax = max(xzw[range(ndxi)])
ybndmin = min(yzw[range(ndxi)])
ybndmax = max(yzw[range(ndxi)])

xvbndmin = xbndmin
xvbndmax = xbndmax
yvbndmin = ybndmin
yvbndmax = 700.

# =============================================================================
# # # # enabling/disabling processes
# =============================================================================
# in-canopy flow
cft = False
# thermal boundary layer
tbl = False

# # process dependencies
if not cft:
    tbl = False

# =============================================================================
# # # # model constants
# =============================================================================
# # light
# light-attenuation coefficient
Kd0 = .1
# maximum saturation intensity (umol photons m^-2 s^-1)
Ikmax = 400.

# # hydrodynamics
# Smagorinsky constant
Cs = .17
# inertia constant
Cm = 1.7
# friction coefficient
Cf = .01
# wave-current angle (degrees)
wcangle = 0.

# # temperature
# TBL coefficient
K0 = 80.
# thermal acclimation coefficient
Kvar = 2.45

# # acidity
# aragonite saturation state
omega0 = 5.

# # morphology
# overall form proportionality constant
Xf = .1
# flow plate proportionality constant
Xpu = .1
# light spacing proportionality constant
XsI = .1
# flow spacing proportionality constant
Xsu = .1

# # dislodgement
# tensile strength substratum [N m^-2]
sigmat = 2e5

# # recruitment
# probability of settlement [-]
probSett = 1e-4
# larvae due to spawning [-]
NoLarvae = 1e6
# larval diameter [m]
dLarv = 1e-3

# =============================================================================
# # # # initial conditions
# =============================================================================
# # initial morphology
dc0 = .05  # m
hc0 = .3  # m
bc0 = .5 * dc0
tc0 = .5 * hc0
ac0 = .2  # m

# # carrying capacity
K = np.zeros(ndxi)
K[np.logical_and.reduce((xzw >= xvbndmin,
                         xzw <= xvbndmax,
                         yzw >= yvbndmin,
                         yzw <= yvbndmax))] = 1.

# # coral cover
P0 = np.array([
        K,
        np.zeros(K.shape),
        np.zeros(K.shape),
        np.zeros(K.shape)
    ]).transpose()

print('Constants set.\n')

# =============================================================================
# # # # intermediate plotting function
# =============================================================================
# triangulate face coordinates for plotting
face_triang = tri.Triangulation(xzw[range(ndxi)], yzw[range(ndxi)])
# define basic plotting routine for model output


def showfld(fld, face_triang, lvls, ttl,
            show=True, save=False, unit=None):
    """
    Function to visualise the intermediate results from the model computations
    in which the grid structure used in Delft3D-FM is translated such that a
    field is presented.
    """
    f = plt.figure()
    plt.tricontourf(face_triang, fld, levels=lvls)
    plt.title(ttl)
    cbar = plt.colorbar(shrink=.5, extend='both', aspect=10)
    if unit is not None:
        cbar.set_label(unit, rotation=270, va='bottom')
    if show:
        plt.show()
    if save:
        f.savefig(os.path.join(figfolder, ('{0}.png'.format(ttl))),
                  bbox_inches='tight')

# =============================================================================
# # # # input modification
# =============================================================================

global spacetime

class DataReshape():
    
    def __init__(self, spacetime):
        """
        Reshape data to matrix with shape [space x time].

        Parameters
        ----------
        spacetime: list, np.ndarray
            Definition of with dimensions of spacetime: [space, time].

        Returns
        -------
        2D-matrix.
        """
        self.space = int(spacetime[0])
        self.time = int(spacetime[1])
    
    def param2matrix(self, x, dimension):
        """Transform parameter to 2D-matrix."""
        # # dimensions
        dims = ['space', 'time']
        if dimension not in dims:
            raise ValueError(
                'Dimension not in {}'.format(dims))
        # # transformation
        x = DataReshape.param2array(x)
        if dimension == 'space':
            M = np.tile(x, (self.time, 1)).transpose()
        elif dimension =='time':
            M = np.tile(x, (self.space, 1))
        return M
    
    def param2array(x):
        """Transform parameter to array if float or integer."""
        if isinstance(x, (float, int)):
            x = np.array([float(x)])
        return x

# =============================================================================
# # # # definition of basic functions for coral development
# =============================================================================

class Environment():
    
    def __init__(self, light=None, LAC=None, temperature=None, acidity=None,
                 stormcat=None):
        self.light = light
        self.Kd = LAC
        self.temp = temperature
        self.acid = acidity
        self.stormcat = stormcat
        
    @property
    def tempK(self):
        """Temperature in Kelvin."""
        if all(self.temp) < 100.:
            return self.temp + 273.15
        else:
            return self.temp
    
    @property
    def tempC(self):
        """Temperature in Celsius."""
        if all(self.temp) > 100.:
            return self.temp - 273.15
        else:
            return self.temp
    
    @property
    def tempMMM(self):
        MM = self.tempK.groupby([self.tempK.index.year,
                                 self.tempK.index.month]).agg(['mean'])
        MMM = MM.groupby(level=0).agg(['min', 'max'])
        MMM.columns = MMM.columns.droplevel([0, 1])
        return MMM
    
    @property
    def dates(self):
        d = self.temp.reset_index().drop('sst', axis=1)
        return pd.to_datetime(d['date'])
    
    def fromFile(self, param, file, fdir=None):
        def date2index(self):
            """Function applicable to time-series in Pandas."""
            self['date'] = pd.to_datetime(self['date'])
            self.set_index('date', inplace=True)
        if fdir is None:
            f = file
        else:
            f = os.path.join(fdir, file)
        
        if param == 'light':
            self.light = pd.read_csv(f, sep='\t')
            date2index(self.light)
        elif param == 'LAC':
            self.Kd = pd.read_csv(f, sep='\t')
            date2index(self.Kd)
        elif param == 'temperature':
            self.temp = pd.read_csv(f, sep='\t')
            date2index(self.temp)
        elif param == 'acidity':
            self.acid = pd.read_csv(f, sep='\t')
            date2index(self.acid)
        elif param == 'storm':
            self.stormcat = pd.read_csv(f, sep='\t')
            self.stormcat.set_index('year', inplace=True)
        
        
        
class Morphology():
    def __init__(self, dc, hc, bc, tc, ac):
        self.dc = dc
        self.hc = hc
        self.bc = bc
        self.tc = tc
        self.ac = ac
    
    def __repr__(self):
        return "Morphology({}, {}, {}, {}, {})".format(
            self.dc, self.hc, self.bc, self.tc, self.ac)
    
    def __str__(self):
        return ('Coral morphology with: '
                'dc = {}; hc = {}, bc = {}, tc = {}, ac = {}'
                .format(self.dc, self.hc, self.bc, self.tc, self.ac))
    
    @property
    def rf(self):
        return self.hc / self.dc
    
    @property
    def rp(self):
        return self.bc / self.dc
    
    @property
    def rs(self):
        return self.dc / self.ac
    
    @property
    def dcRep(self):
        return (self.bc * (self.hc - self.tc) + self.dc * self.tc) / self.hc
    
    @property
    def volume(self):
        return .25 * np.pi * (
            (self.hc - self.tc) * self.bc ** 2 + self.tc * self.dc ** 2)
    
    @property
    def dcMatrix(self):
        M = DataReshape(spacetime)
        return M.param2matrix(self.dc, 'space')
    
    @property
    def hcMatrix(self):
        M = DataReshape(spacetime)
        return M.param2matrix(self.hc, 'space')
    
    @property
    def bcMatrix(self):
        M = DataReshape(spacetime)
        return M.param2matrix(self.bc, 'space')
    
    @property
    def tcMatrix(self):
        M = DataReshape(spacetime)
        return M.param2matrix(self.tc, 'space')
    
    @property
    def acMatrix(self):
        M = DataReshape(spacetime)
        return M.param2matrix(self.ac, 'space')
    
    @property
    def dcRepMatrix(self):
        M = DataReshape(spacetime)
        return M.param2matrix(self.dcRep, 'space')
    
    @volume.setter
    def volume(self, volume):
        def Vc2dc(self):
            """Coral volume > diameter of the plate."""
            self.dc = ((4. * volume) / (np.pi * rf * rp * (
                1. + rp - rp ** 2))) ** (1. / 3.)
        def Vc2hc(self):
            """Coral volume > coral height."""
            self.hc = ((4. * volume * rf ** 2) / (np.pi * rp * (
                1. + rp - rp ** 2))) ** (1. / 3.)
        def Vc2bc(self):
            """Coral volume > diameter of the base."""
            self.bc = ((4. * volume * rp ** 2) / (np.pi * rf * (
                1. + rp - rp ** 2))) ** (1. / 3.)
        def Vc2tc(self):
            """Coral volume > thickness of the plate."""
            self.tc = ((4. * volume * rf ** 2 * rp ** 2) / (np.pi * (
                1. + rp - rp ** 2))) ** (1. / 3.)
        def Vc2ac(self):
            """Coral volume > axial distance."""
            self.ac = (1. / rs) * ((4. * volume) / (np.pi * rf * rp * (
                1. + rp - rp ** 2))) ** (1. / 3.)
        # # obtain previous morphological ratios
        rf = self.rf
        rp = self.rp
        rs = self.rs
        # # update morphology
        Vc2dc(self)
        Vc2hc(self)
        Vc2bc(self)
        Vc2tc(self)
        Vc2ac(self)
    
    
    def morph2vegden(self):
        """Translation from morphological dimensions to vegetation density."""
        try:
            rnveg = np.zeros(self.ac.shape)
            rnveg[self.ac > 0.] = (
                2. * self.dcRep[self.ac > 0.]) / (self.ac[self.ac > 0.] ** 2)
        except TypeError:
            if self.ac > 0.:
                rnveg = (2. * self.dcRep) / (self.ac ** 2)
            else:
                rnveg = 0.
        return rnveg
            
    
    def plot(self, axLabels=False, axTicks=False, explanation=False,
             save=False, figname=None, figdir=None):
        def outerLines(self, x0):
            x = [x0 - .5 * self.bc,
                 x0 - .5 * self.bc,
                 x0 - .5 * self.dc,
                 x0 - .5 * self.dc,
                 x0 + .5 * self.dc,
                 x0 + .5 * self.dc,
                 x0 + .5 * self.bc,
                 x0 + .5 * self.bc]
            y = [0.,
                 self.hc - self.tc,
                 self.hc - self.tc,
                 self.hc,
                 self.hc,
                 self.hc - self.tc,
                 self.hc - self.tc,
                 0.]
            return x, y
        
        def annotateText(ax, xyfrom, xyto, text=None):
            if text is None:
                text = str(np.sqrt((xyfrom[0] - xyto[0]) ** 2 +
                                   (xyfrom[1] - xyto[1]) ** 2))
            ax.annotate(
                '', xyfrom, xyto,
                arrowprops=dict(arrowstyle='<->'))
            ax.text(
                (xyto[0] + xyfrom[0]) / 2., .01 + (xyto[1] + xyfrom[1]) / 2.,
                text, fontsize=12)
        left = outerLines(self, 0.)
        right = outerLines(self, self.ac)
        
        fig, ax = plt.subplots(figsize=(8, 4))
        # zero-line
        ax.plot([-.5 * self.ac, 1.5 * self.ac], [0., 0.],
                color='gray', alpha=1.,
                linewidth=2., linestyle='solid',
                label='_nolegend_')
        # plot data
        ax.plot(left[0], left[1],
                color='black', alpha=1.,
                linewidth=1., linestyle='solid',
                label='_nolegend_')
        ax.plot(right[0], right[1],
                color='gray', alpha=.5,
                linewidth=1., linestyle='dashed',
                label='_nolegend_')
        # axes labels
        if axLabels:
            ax.set_xlabel('horizontal distance [m]')
            ax.set_ylabel('vertical distance [m]')
        # axes ticks
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if not axTicks:
            plt.tick_params(
                axis='x',
                which='both',
                bottom=False,
                top=False,
                labelbottom=False)
            plt.tick_params(
                axis='y',
                which='both',
                left=False,
                right=False,
                labelleft=False)
        # plot limits
        ax.set_xlim([-.5 * self.ac, 1.5 * self.ac])
        ax.set_ylim([0., 1.2 * self.hc])
        # explaining lines / texts / etc.
        if explanation:
            # dc
            annotateText(
                ax, (-.5 * self.dc, 1.15 * self.hc),
                (.5 * self.dc, 1.15 * self.hc),
                text=r'$d_{{c}}$')
            # hc
            annotateText(
                ax, (-.25 * self.ac, 0.), (-.25 * self.ac, self.hc),
                text=r'$h_{{c}}$')
            # bc
            annotateText(
                ax, (-.5 * self.bc, .5 * (self.hc - self.tc)),
                (.5 * self.bc, .5 * (self.hc - self.tc)),
                text=r'$b_{{c}}$')
            # tc
            annotateText(
                ax, (.25 * self.ac, self.hc - self.tc),
                (.25 * self.ac, self.hc),
                text=r'$t_{{c}}$')
            # ac
            annotateText(
                ax, (0., 1.05 * self.hc), (self.ac, 1.05 * self.hc),
                text=r'$a_{{c}}$')
        # legend / title
        
        # # save figure
        if save:
            if figname is None:
                print('WARNING: No figure name specified.')
                figname = 'Morphology'
            figfile = figname + '.png'
            if figdir is None:
                figfull = figfile
            else:
                figfull = os.path.join(figdir, figfile)
            fig.savefig(figfull, dpi=300, bbox_inches='tight')
            print('Figure saved as: "{}"'.format(figfull))

class Light():
    
    def __init__(self, I0, Kd, h, thetamax=.5*np.pi):
        """
        Light micro-environment.

        Parameters
        ----------
        I0 : float, np.ndarray
            Incoming light-intensity as measured at water-air interface
        [umol photons m-2 s-1].
        Kd : float, np.ndarray
            Light-attenuation coefficient [m-1].
        h : float, np.ndarray
            Water depth (excl. coral canopy) [m].
        thetamax : float, optional
            Maximum spreading of light as measured at water-air interface
            [rad]. The default is .5*np.pi.
        """
        M = DataReshape(spacetime)
        self.I0 = M.param2matrix(I0, 'time')
        self.Kd = M.param2matrix(Kd, 'time')
        self.h = M.param2matrix(h, 'space')
        self.thetamax = thetamax
    
    def repLight(self, coral):
        """Representative light-intensity."""
        Light.biomass(self, coral)
        # light catchment top of plate
        top = .25 * np.pi * coral.dcMatrix ** 2 * self.I0 * np.exp(
            - self.Kd * (self.h - coral.hcMatrix))
        # light catchment side of plate
        side = (np.pi * coral.dcMatrix * self.I0) / self.Kd * (
            np.exp(- self.Kd * (self.h - coral.hcMatrix)) -
            np.exp(- self.Kd * (self.h - coral.hcMatrix + coral.tcMatrix)))
        # light catchment side of base
        base = (np.pi * coral.bcMatrix * self.I0) / self.Kd * (
            np.exp(- self.Kd * (self.h - self.L)) -
            np.exp(- self.Kd * self.h))
        # total light catchment
        total = top + side + base
        # biomass-averaged
        try:
            coral.light = self.I0 * np.exp(- self.Kd * self.h)
            coral.light[coral.Bc > 0.] = (
                total[coral.Bc > 0.] / coral.Bc[coral.Bc > 0.])
        except TypeError:
            if coral.Bc > 0.:
                coral.light = total / coral.Bc
            else:
                coral.light = self.I0 * np.exp(- self.Kd * self.h)
    
    def biomass(self, coral):
        """Coral biomass."""
        Light.baseLight(self, coral)
        coral.Bc = np.pi * (.25 * coral.dcMatrix ** 2 +
                            coral.dcMatrix * coral.tcMatrix +
                            coral.bcMatrix * self.L)
    
    def baseLight(self, coral):
        def lightSpreading(self):
            theta = self.theta_max * np.exp(
                - self.light_attenuation * (self.h - coral.hcMatrix + coral.tcMatrix))
            return theta
        theta = lightSpreading(self)
        self.L = coral.hcMatrix - coral.tcMatrix - (
            coral.dcMatrix - coral.bcMatrix) / (
            2. * np.tan(.5 * theta))
        try:
            self.L[self.L < 0.] = 0.
        except TypeError:
            if self.L < 0.:
                self.L = 0.
    
class Flow():
    
    def __init__(self, Ucurr, Uwave, h=None, Tp=None, Cs=.17, Cm=1.7, Cf=.01,
                 nu=1e-6, alpha=1e-7, psi=2., wcAngle=0., rd = 500., theta=.5,
                 err=1e-6, maxiter_k=1e5, maxiter_aw=1e5):
        """
        Flow micro-environment.

        Parameters
        ----------
        Ucurr : float, np.ndarray
            Current-induced flow velocity [m s-1].
        Uwave : float, np.ndarray
            Wave-induced flow velocity [m s-1].
        h : float, np.ndarray, optional
            Water depth (excl. coral canopy) [m]. The default is None.
        Tp : float, np.ndarray, optional
            Peak wave period [s]. The default is None.
        Cs : float, optional
            Smagorinsky coefficient [-]. The default is .17.
        Cm : float, optional
            Inertia coefficient [-]. The default is 1.7.
        Cf : float, optional
            Friction coefficient [-]. The default is .01.
        nu : float, optional
            Kinematic viscosity of water [m2 s-1]. The default is 1e-6.
        alpha : float, optional
            Thermal diffusivity of water [m2 s-1]. The default is 1e-7.
        psi : float, optional
            Ratio of lateral over streamwise spacing of corals [-]. The
            default is 2..
        wcAngle : float, optional
            Angle between current- and wave-induced flows [rad]. The default
            is 0..
        rd : float, optional
            Velocity boundary layer wall-coordinate [-]. The default is 500..
        theta : float, optional
            Update ratio for above-canopy flow [-]. The default is .5.
        err : float, optional
            Maximum allowable relative error [-]. The default is 1e-6.
        maxiter_k : integer, float, optional
            Maximum number of iterations taken over the canopy layers. The
            default is 1e5.
        maxiter_aw : integer, float, optional
            Maximum number of iterations to solve the complex-valued wave-
            attenuation coefficient. The default is 1e5.
        """
        self.uc = DataReshape.param2array(Ucurr)
        self.uw = DataReshape.param2array(Uwave)
        self.h = DataReshape.param2array(h)
        self.Tp = DataReshape.param2array(Tp)
        self.Cs = Cs
        self.Cm = Cm
        self.Cf = Cf
        self.nu = nu
        self.alpha = alpha
        self.psi = psi
        self.wcAngle = wcAngle
        self.rd = rd
        self.theta = theta
        self.err = err
        self.maxiter_k = int(maxiter_k)
        self.maxiter_aw = int(maxiter_aw)
    
    def waveCurrent(self, coral, incanopy=True):
        """Wave-current interaction."""
        if incanopy:
            alphaw = Flow.waveAttenuation(self, coral, wac='wave')
            alphac = Flow.waveAttenuation(self, coral, wac='current')
            coral.ucm = np.sqrt((alphaw * self.uw) ** 2 +
                               (alphac * self.uc) ** 2 +
                               (2. * alphaw * self.uw * alphac * self.uc *
                                np.cos(self.wcAngle)))
        else:
            coral.um = np.sqrt(self.uw ** 2 + self.uc ** 2 +
                               2. * self.uw * self.uc *
                               np.cos(self.wcAngle))
    
    def waveAttenuation(self, coral, wac):
        """Wave-attenuation coefficient."""
        # # function and derivative defintions
        def func(beta):
            # plt.scatter(beta.real, beta.imag, alpha=.5)
            # components
            shear = (8. * af) / (3. * np.pi) * (
                abs(1. - beta) * (1. - beta)) / Lshear[i]
            drag = (8. * af) / (3. * np.pi) * (abs(beta) * beta) / Ldrag
            inertia = 1j * beta * ((self.Cm * lambdap[i]) / (1. - lambdap[i]))
            # combined
            f = 1j * (beta - 1.) - shear + drag + inertia
            # output
            return f
        
        def deriv(beta):
            # components
            shear = (((1. - beta) ** 2 / abs(1. - beta) - abs(1. - beta)) /
                     Lshear[i])
            drag = (beta ** 2 / abs(beta) + beta) / Ldrag
            inertia = 1j * (self.Cm * lambdap[i]) / (1. - lambdap[i])
            # combined
            df = 1j + (8. * af) / (3. * np.pi) * (- shear + drag) + inertia
            # output
            return df
        
        # # parameter definitions
        # wave- vs. current-induced
        if wac == 'wave':
            T = self.Tp
            U = self.uw
        elif wac == 'current':
            T = 1e3 * np.ones(spacetime[0])
            U = self.uc
        # geometric parameters
        Ap = .25 * np.pi * coral.dcRep ** 2
        Af = coral.dcRep * coral.hc
        AT = .5 * coral.ac ** 2
        lambdap = Ap / AT
        lambdaf = Af / AT
        Lshear = coral.hc / (self.Cs ** 2)
        
        # # calculations - single layer canopy
        if spacetime[0] == 1:
            lambdap = [lambdap]
            lambdaf = [lambdaf]
            Lshear = [Lshear]
        aw = np.ones(spacetime[0])
        for i in range(spacetime[0]):
            # no corals
            if coral.hc[i] == 0.:
                aw[i] = 1.
            # emergent
            elif self.h[i] <= coral.hc[i]:
                up = U[i] / (1. - lambdap[i])
                uc = ((1. - lambdap[i]) / (1. - np.sqrt(4. * lambdap[i]) /
                       (self.psi * np.pi))) * up
                Re = (uc * coral.dcRep[i]) / self.nu
                Cd = 1. + 10. * Re ** (- 2. / 3.)
                aw[i] = 1.
            # submerged
            else:
                # initial values for iteration
                uf = U[i]
                Cd = 1.
                # iteration
                for k in range(int(self.maxiter_k)):
                    Ldrag = (2. * coral.hc[i] * (1. - lambdap[i])) / (
                        Cd * lambdaf[i])
                    af = (uf * T[i]) / (2. * np.pi)
                    if wac == 'wave':
                        aw[i] = abs(newton(
                            func, x0=complex(.5, .5), fprime=deriv,
                            maxiter=self.maxiter_aw))
                    elif wac == 'current':
                        X = Ldrag / Lshear[i] * (
                            coral.hc[i] / (self.h[i] - coral.hc[i]) + 1.)
                        aw[i] = (X - np.sqrt(X)) / (X - 1.)
                    up = aw[i] * uf
                    uc = (1. - lambdap[i]) / (
                        1. - np.sqrt((4. - lambdap[i]) / (
                            self.psi * np.pi))) * up
                    Re = (uc * coral.dcRep[i]) / self.nu
                    Cdk = 1. + 10. * Re ** (-2. / 3.)
                    if abs((Cdk - Cd) / Cdk) <= self.err:
                        break
                    else:
                        Cd = float(Cdk)
                        uf = abs((1. - self.theta) * uf + self.theta * (
                                     self.h[i] * U[i] - coral.hc[i] * up) /
                                 (self.h[i] - coral.hc[i]))
                    if k == self.maxiter_k:
                        print(
                            'WARNING: maximum number of iterations reached: {}'
                            .format(self.maxiter_k))
        return aw
                    
    def TBL(self, coral):
        """Thermal boundary layer."""
        Flow.VBL(self, coral)
        coral.deltat = self.delta * ((self.alpha / self.nu) ** (1. / 3.))
        
    
    def VBL(self, coral):
        """Velocity boundary layer."""
        try:
            self.delta = np.zeros(coral.ucm.shape)
            self.delta[coral.ucm > 0.] = (
                (self.rd * self.nu) / (
                    np.sqrt(self.Cf) * coral.ucm[coral.ucm > 0.]))
        except TypeError:
            if coral.ucm > 0.:
                self.delta = (self.rd * self.nu) / (
                    np.sqrt(self.Cf) * coral.ucm)
            else:
                self.delta = 0.

class Temperature():
    
    def __init__(self, T, K0=80., ap=.4, k=.6089):
        """
        Thermal micro-environment.

        Parameters
        ----------
        h : float, np.ndarray
            Water depth (excl. coral canopy) [m].
        T : float, np.ndarray
            Temperature of water [K].
        K0 : float, optional
            Morphological thermal coefficient. The default is 80..
        ap : float, optional
            Absorptibity of coral [-]. The default is .4.
        k : float, optional
            Thermal conductivity [J m-1 s-1 K-1]. The default is .6089.
        """
        M = DataReshape(spacetime)
        self.T = M.param2matrix(T, 'time')
        self.K0 = K0
        self.ap = ap
        self.k = k
    
    def coralTemperature(self, coral, TME=True):
        """Coral temperature."""
        M = DataReshape(spacetime)
        deltat = M.param2matrix(coral.deltat, 'space')
        if TME:
            coral.temp = self.T + (
                (deltat * self.ap) / (self.k * self.K0)) * coral.light
        else:
            coral.temp = self.T

class Photosynthesis():
    
    def __init__(self, I0, iota=.6, Ikmax=372.32, Pmmax=1., betaI=.34,
                 betaP=.09, Ea=6e4, R=8.31446261815324, Kvar=2.45, nn=60.,
                 Pumin=.68886964, ucr=.17162374,
                 firstYear=True, TME=False):
        """
        Photosynthetic efficiency based on photosynthetic dependencies.

        Parameters
        ----------
        I0 : float, np.ndarry
            Incoming light-intensity as measured at water-air interface
            [umol photons m-2 s-1].
        iota : float, optional
            Photo-acclimation rate [d-1]. The default is .6.
        Ikmax : float, optional
            Maximum value of the quasi steady-state for the saturation light-
            intensity [umol photons m-2 s-1]. The default is 372.32.
        Pmmax : float, optional
            Maximum value of the quasi steady-state for the maximum
            photosynthetic efficiency [-]. The default is 1..
        betaI : float, optional
            Exponent of the quasi steady-state for the saturation light-
            intensity [-]. The default is .34.
        betaP : float, optional
            Exponent of the quasi steady-state for the maximum photosynthetic
            efficiency [-]. The default is .09.
        Ea : float, optional
            Activation energy [J mol-1]. The default is 6e4.
        R : float, optional
            Gas constant [J K-1 mol-1]. The default is 8.31446261815324.
        Kvar : TYPE, optional
            Thermal-acclimation coefficient [-]. The default is 2.45.
        nn : TYPE, optional
            Theraml-acclimation period [yrs]. The default is 60..
        Pumin : TYPE, optional
            Min. photosynthetic flow dependency [-]. The default is .68886964.
        ucr : TYPE, optional
            Minimum flow velocity at which photosynthesis is not limited by
            flow. The default is .17162374.
        firstYear : bool, optional
            First year of the model simulation. The default is True.
        TME : bool, optional
            Include the thermal micro-environment. The default is False.
        """
        M = DataReshape(spacetime)
        self.I0 = M.param2matrix(I0, 'time')
        self.iota = iota
        self.Ikmax = Ikmax
        self.Pmmax = Pmmax
        self.betaI = betaI
        self.betaP = betaP
        self.Ea = Ea
        self.R = R
        self.Kvar = Kvar
        self.nn = nn
        self.Pumin = Pumin
        self.ucr = ucr
        self.firstYear = firstYear
        self.TME = TME
    
    def photosyn(self, coral, environment, year):
        """Photosynthetic efficiency."""
        M = DataReshape(spacetime)
        # components
        Photosynthesis.PLD(self, coral, 'qss')
        Photosynthesis.PTD(self, coral, environment, year)
        Photosynthesis.PFD(self, coral)
        # combined
        coral.photo_rate = self.pld * self.ptd * M.param2matrix(
            self.pfd, 'space')
    
    def PLD(self, coral, output):
        """Photosynthetic light dependency."""
        def photoacc(self, Xold, param):
            """Photo-acclimation."""
            # # parameter definition
            if param == 'Ik':
                Xmax = self.ik_max
                betaX = self.betaI
            elif param == 'Pmax':
                Xmax = self.pm_max
                betaX = self.betaP
            # # calculations
            XS = Xmax * (coral.light / self.I0) ** betaX
            if output =='qss':
                return XS
            elif output == 'new':
                Xnew = XS + (Xold - XS) * np.exp(- self.iota)
                return Xnew
        # # parameter definition
        if output == 'qss':
            Ik = photoacc(self, None, 'Ik')
            Pmax = photoacc(self, None, 'Pmax')
        elif output == 'new':
            NotImplemented
        # # calculations
        self.pld = Pmax * (np.tanh(coral.light / Ik) -
                           np.tanh(.01 * self.I0 / Ik))
        
    
    def PTD(self, coral, environment, year):
        """Photosynthetic thermal dependency."""
        def thermacc(self):
            """Thermal-acclimation."""
            MMM = environment.tempMMM[np.logical_and(
                environment.tempMMM.index < year,
                environment.tempMMM.index >= year - int(self.nn / coral.Csp))]
            if self.TME:
                if self.firstYear:
                    NotImplemented
                else:
                    NotImplemented
            mmin, mmax = MMM.mean(axis=0)
            smin, smax = MMM.std(axis=0)
            self.Tlo = mmin - self.k_var * smin
            self.Thi = mmax + self.k_var * smax
        
        def adaptTemp(self):
            """Adapted temperature response."""
            def spec(self):
                """Specialisation term."""
                self.spec = 4e-4 * np.exp(-.33 * (self.DT - 10.))
            spec(self)
            self.f1 = - ((coral.temp - self.Tlo) *
                         ((coral.temp - self.Tlo) ** 2 - self.DT ** 2))
            Tcr = self.Tlo - (1. / np.sqrt(3.)) * self.DT
            try:
                if self.TME:
                    self.f1[coral.temp <= Tcr] = - (
                        (2. / (3. * np.sqrt(3.))) *
                        self.DT[coral.temp <= Tcr] ** 3)
                else:
                    self.f1[coral.temp <= Tcr] = -(
                        (2. / (3. * np.sqrt(3.))) * self.DT ** 3)
            except TypeError:
                if coral.temp <= Tcr:
                    self.f1 = (2. / (3. * np.sqrt(3.))) * self.DT ** 3
            self.f1 *= self.spec
        
        def thermEnv(self):
            """Thermal envelope."""
            self.f2 = np.exp((self.Ea / self.R) * (1. / 300. - 1. / self.Topt))
        
        # # parameter definition
        thermacc(self)
        self.DT = self.Thi - self.Tlo
        self.Topt = self.Tlo + (1. / np.sqrt(3.)) * self.DT
        # # calculations
        # components
        adaptTemp(self)
        thermEnv(self)
        # combined
        self.ptd = self.f1 * self.f2
        
    
    def PFD(self, coral):
        """Photosynthetic flow dependency."""
        self.pfd = self.Pumin + (1. - self.Pumin) * np.tanh(
            2. * coral.ucm / self.ucr)


class PopStates():
    
    def __init__(self, K, rG=.002, rR=.2, rM=.04, rB=8.):
        """
        Population dynamics.

        Parameters
        ----------
        K : float, np.ndarray
            Carrying capacity (0 <= K <= 1) [-].
        rG : float, optional
            Growth rate [d-1]. The default is .002.
        rR : float, optional
            Recovering rate [d-1]. The default is .2.
        rM : float, optional
            Mortality rate [d-1]. The default is .04.
        rB : float, optional
            Bleaching rate [d-1]. The default is 8..
        """
        self.K = DataReshape.param2array(K)
        self.rG = rG
        self.rR = rR
        self.rM = rM
        self.rB = rB
        
    def popstates_t(self, coral, dt=1., threshold=0.):
        """Population dynamics over time."""
        coral.pop_states = np.zeros((spacetime[0], spacetime[1], 4))
        for n in range(spacetime[1]):
            coral.pop_states[self.K > 0, n, :] = PopStates.popstates_xy(
                    self, coral, coral.photo_rate[self.K > 0, n],
                    dt=dt, threshold=threshold)
            coral.p0[self.K > 0, :] = coral.pop_states[self.K > 0, n, :]
    
    def popstates_xy(self, coral, PS, dt=1., threshold=0.):
        """Population dynamics over space."""
        P = np.zeros((len(self.K[self.K > 0]), 4))
        # # calculations
        # growing conditions
        # > bleached pop.
        P[PS > threshold, 3] = (
                coral.p0[PS > threshold, 3] / (1. + dt * (
                8. * self.rR * PS[PS > threshold] /
                coral.Csp + self.rM * coral.Csp)))
        # > pale pop.
        P[PS > threshold, 2] = (
            (coral.p0[PS > threshold, 2] + (
                8. * dt * self.rR * PS[PS > threshold] / coral.Csp) *
             P[PS > threshold, 3]) / (
                 1. + dt * self.rR * PS[PS > threshold] * coral.Csp))
        # > recovering pop.
        P[PS > threshold, 1] = (
            (coral.p0[PS > threshold, 1] + dt * self.rR * PS[PS > threshold] *
             coral.Csp * P[PS > threshold, 2]) /
            (1. + .5 * dt * self.rR * PS[PS > threshold] * coral.Csp))
        # > healthy pop.
        a = (dt * self.rG * PS[PS > threshold] * coral.Csp /
             self.K[PS > threshold])
        b = (1. - dt * self.rG * PS[PS > threshold] * coral.Csp * (
                1. - P[PS > threshold, 1:].sum(axis=1) /
                self.K[PS > threshold]))
        c = - (coral.p0[PS > threshold, 0] +
               .5 * dt * self.rR * PS[PS > threshold] *
               coral.Csp * P[PS > threshold, 1])
        P[PS > threshold, 0] = (
            -b + np.sqrt(b ** 2 - 4 * a * c)) / (2. * a)
        
        # bleaching conditions
        # > healthy pop.
        P[PS <= threshold, 0] = (
                coral.p0[PS <= threshold, 0] / (
                1. - dt * self.rB * PS[PS <= threshold] * coral.Csp))
        # > recovering pop.
        P[PS <= threshold, 1] = (
                coral.p0[PS <= threshold, 1] / (
                1. - dt * self.rB * PS[PS <= threshold] * coral.Csp))
        # > pale pop.
        P[PS <= threshold, 2] = (
            (coral.p0[PS <= threshold, 2] - dt * self.rB *
             PS[PS <= threshold] * coral.Csp * (
                 P[PS <= threshold, 0] + P[PS <= threshold, 1])) /
            (1. - .5 * dt * self.rB * PS[PS <= threshold] * coral.Csp))
        # > bleached pop.
        P[PS <= threshold, 3] = (
            (coral.p0[PS <= threshold, 3] -
             .5 * dt * self.rB * PS[PS <= threshold] * coral.Csp *
             P[PS <= threshold, 2]) /
            (1. - .25 * dt * self.rB * PS[PS <= threshold] * coral.Csp))
        
        # # check on carrying capacity
        if any(P.sum(axis=1) > 1.0001 * self.K):
            print('WARNING: total population than carrying capacity at {}. '
                  '(PT = {}, K = {})'
                  .format(
                      np.arange(len(self.K))[P.sum(axis=1) > 1.0001 * self.K],
                      P[P.sum(axis=1) > 1.0001 * self.K],
                      self.K[P.sum(axis=1) > 1.0001 * self.K]))
        
        # # output
        return P


class Calcification():
    
    def __init__(self, gC=.5, omega0=.14587415, kappaA=.66236107):
        """
        Calcification rate.

        Parameters
        ----------
        gC : float
            Calcification constant [kg m-2 d-1].
        omega0 : float
            Aragonite dissolution state [-].
        kappaA : float
            Modified Michaelis-Menten half-rate coefficient [-].
        """
        self.gC = gC
        self.omega0 = omega0
        self.kappaA = kappaA
    
    def calRate(self, coral, omega):
        """Calcification rate."""
        def aragoniteDependency(self):
            """Aragonite dependency."""
            self.ad = (omega - self.omega0) / (
                self.kappaA + omega - self.omega0)
            M = DataReshape(spacetime)
            self.ad = M.param2matrix(self.ad, 'time')
        
        aragoniteDependency(self)
        coral.calc = (
                self.gC * coral.Csp * coral.pop_states[:, :, 0] *
                self.ad * coral.photo_rate)


class MorDevelopment():
    
    def __init__(self, Gsum, h, Kd,
                 Xf=.1, Xp=.5, Xpu=.1, Xs=.5/np.sqrt(2.), XsI=.1, Xsu=.1,
                 ucr=.172, thetamax=.5*np.pi, rhoc=1600., dtyear=1.):
        """
        Morphological development.

        Parameters
        ----------
        Gsum : float, np.ndarray
            Accumulation of calcification of [dtyear] years [kg m-2 yr-1].
        h : float, np.ndarray
            Water depth (excl. coral canopy) [m].
        Kd : float, np.ndarray
            Light-attenuation coefficient [m-1].
        Xf : float, optional
            Overall form proportionality constant [-]. The default is .1.
        Xp : float, optional
            Overall plate proportionality constant [-]. The default is .5.
        Xpu : float, optional
            Flow plate proportionality constant [-]. The default is .1.
        Xs : float, optional
            Overall spacing proportionality constant [-].
            The default is .5/np.sqrt(2.).
        XsI : float, optional
            Light spacing proportionality constant [-]. The default is .1.
        Xsu : float, optional
            Flow spacing proportionality constant [-]. The default is .1.
        ucr : float, optional
            Critical flow velocity [m s-1]. The default is .172.
        thetamax : float, optional
            Maximum spreading of light as measured at water-air interface
            [rad]. The default is .5*np.pi.
        rhoc : float, optional
            Density of coral [kg m-3]. The default is 1600..
        dtyear : float, optional
            Update interval [yr]. The default is 1..
        """
        M = DataReshape(spacetime)
        self.Gsum = Gsum
        self.h = M.param2matrix(h, 'space')
        self.Kd = M.param2matrix(Kd, 'time')
        self.Xf = Xf
        self.Xp = Xp
        self.Xpu = Xpu
        self.Xs = Xs
        self.XsI = XsI
        self.Xsu = Xsu
        self.ucr = ucr
        self.thetamax = thetamax
        self.rhoc = rhoc
        self.dtyear = dtyear
    
    def update(self, coral, I0):
        """Update morphology."""
        M = DataReshape(spacetime)
        self.I0 = M.param2matrix(I0, 'time')
        # # optimal morphological ratios
        def rfOptimal(self):
            """Optimal form ratio."""
            self.rfOpt = self.prop_form * (
                coral.light.mean(axis=1) / self.I0.mean(axis=1)) / (
                self.ucr / coral.ucm)
        def rpOptimal(self):
            """Optimal plate ratio."""
            self.rpOpt = self.prop_plate * (1. + np.tanh(
                self.prop_plate_flow * (coral.ucm - self.ucr) / self.ucr))
        def rsOptimal(self):
            """Optimal spacing ratio."""
            self.rsOpt = self.prop_space * (
                1. - np.tanh(
                self.prop_space_light * coral.light.mean(axis=1) /
                self.I0.mean(axis=1))) * (
                    1. + np.tanh(self.prop_space_flow * (coral.ucm - self.ucr) / self.ucr))
        
        # # increased coral volume
        def deltaVolume(self):
            """Volumetric change."""
            self.dVc = (
                               .5 * coral.ac ** 2 * self.calc_sum * self.dt_year /
                               self.rho_c) * coral.Bc.mean(axis=1)
        
        # # update morphological ratio
        def ratioUpdate(self, ratio):
            """Update morphological ratios."""
            # input check
            ratios = ['rf', 'rp', 'rs']
            if ratio not in ratios:
                raise ValueError(
                    'Morphological ratio not in {}'
                    .format(ratios))
            # calculations
            deltaVolume(self)
            def PDE(self, rold, ropt):
                """Mass balance."""
                r = (coral.volume * rold + self.dVc * ropt) / (
                    coral.volume + self.dVc)
                return r
            if ratio == 'rf':
                rfOptimal(self)
                self.rf = PDE(self, coral.rf, self.rfOpt)
            elif ratio == 'rp':
                rpOptimal(self)
                self.rp = PDE(self, coral.rp, self.rpOpt)
            elif ratio == 'rs':
                rsOptimal(self)
                self.rs = PDE(self, coral.rs, self.rsOpt)
        
        # # calculations
        # update ratios
        ratios = ['rf', 'rp', 'rs']
        for r in ratios:
            ratioUpdate(self, r)
        # update morphology
        coral.volume += self.dVc


class Dislodgement():
    
    def __init__(self, sigmat, Cd=1., rhow=1025.):
        """
        Dislodgement check.

        Parameters
        ----------
        sigmat : float
            Tensile strength of substratum [N m-2].
        Cd : float, optional
            Drag coefficient [-]. The default is 1..
        rhow : float, optional
            Density of water [kg m-3]. The default is 1025..
        """
        self.sigmat = sigmat
        self.Cd = Cd
        self.rhow = rhow
    
    def update(self, coral, survivalCoefficient=1.):
        """Update morphology due to storm damage."""
        # # partial dislodgement
        Dislodgement.partialDislodgement(self, coral, survivalCoefficient)
        # # update
        # population states
        for s in range(4):
            coral.p0[:, s] *= self.pardis
        # morphology
        coral.volume *= self.pardis
    
    def partialDislodgement(self, coral, survivalCoefficient=1.):
        """Percentage surviving storm event."""
        try:
            self.pardis = np.ones(coral.dc.shape)
            dislodged = Dislodgement.dislodgementCriterion(self, coral)
            self.pardis[dislodged] = survivalCoefficient * (
                    self.dmt[dislodged] / self.csf[dislodged])
        except TypeError:
            if Dislodgement.dislodgementCriterion(self, coral):
                self.pardis = survivalCoefficient * self.dmt / self.csf
            else:
                self.pardis = 1.
    
    def dislodgementCriterion(self, coral):
        """Dislodgement criterion. Returns boolean (array)."""
        Dislodgement.DMT(self, coral)
        Dislodgement.CSF(self, coral)
        return self.dmt <= self.csf
    
    def DMT(self, coral):
        """Dislodgement Mechanical Threshold."""
        try:
            self.dmt = 1e20 * np.ones(coral.um.shape)
            self.dmt[coral.um > 0] = self.sigmat / (
                self.rhow * self.Cd * coral.um[coral.um > 0] ** 2)
        except TypeError:
            if coral.um > 0:
                self.dmt = self.sigmat / (
                    self.rhow * self.Cd * coral.um ** 2)
            else:
                self.dmt = 1e20
    
    def CSF(self, coral):
        """Colony Shape Factor."""
        # arms of moment
        at = coral.hc - .5 * coral.tc
        ab = .5 * (coral.hc - coral.tc)
        # area of moment
        At = coral.dc * coral.tc
        Ab = coral.bc * (coral.hc - coral.tc)
        # integral
        S = at * At + ab * Ab
        # colony shape factor
        try:
            self.csf = np.zeros(coral.dc.shape)
            self.csf[coral.bc > 0] = 16. / (np.pi * coral.bc ** 3) * S
        except TypeError:
            if coral.bc > 0:
                self.csf = 16. / (np.pi * coral.bc ** 3) * S
            else:
                self.csf = 0.


class Recruitment():
    
    def __init__(self, NoLarvae, probSett, dLarv):
        """
        Recruitment dynamics.

        Parameters
        ----------
        NoLarvae : float
            Number of larvae released during spawning [l m-2].
        probSett : float
            Probability of settlement [-].
        dLarv : float
            Larval diameter [m].
        """
        self.Nl = NoLarvae
        self.ps = probSett
        self.dl = dLarv
    
    def update(self, coral):
        """Update coral cover / volume after spawning event."""
        coral.p0[0] += Recruitment.spawning(self, coral, 'P')
        coral.volume += Recruitment.spawning(self, coral, 'V')
    
    def spawning(self, coral, param):
        """Contribution due to mass coral spawning."""
        # # input check
        params = ['P', 'V']
        if param not in params:
            raise ValueError(
                'Parameter definition not in {}'
                .format(params))
        
        # # calculations
        # potential
        if param == 'P':
            S = self.ps * self.Nl * self.dl ** 2
        elif param == 'V':
            S = self.ps * self.Nl * self.dl ** 3
        # recruitment
        self.PHav = coral.pop_states[:, -1, 0].mean()
        recr = np.zeros(coral.K.shape)
        recr[coral.K > 0] = S * self.PHav * (
                1. - coral.pop_states[coral.K > 0, -1, :].sum(axis=1) /
                coral.K[coral.K > 0])
        
        # # output
        return recr

print('Functions defined.\n')

# =============================================================================
# # # # environmental time-series
# =============================================================================
env = Environment()
env.fromFile('light', fLight, inputdir)
try:
    env.fromFile('LAC', fLAC, inputdir)
except FileNotFoundError:
    env.Kd = pd.DataFrame(index=env.light.index,
                                  data={'lac': Kd0})
env.fromFile('temperature', fTemperature, inputdir)
try:
    env.fromFile('acidity', fAcidity, inputdir)
except FileNotFoundError:
    env.acid = pd.DataFrame(index=env.light.index,
                                    data={'arag': omega0})
env.fromFile('storm', fStorm, inputdir)

print('Environmental conditions defined.\n')

# =============================================================================
# # # # coral definition
# =============================================================================
coral = Morphology(dc0 * np.ones(ndxi),
                   hc0 * np.ones(ndxi),
                   bc0 * np.ones(ndxi),
                   tc0 * np.ones(ndxi),
                   ac0 * np.ones(ndxi))
coral.Csp = 1.
coral.K = K
coral.P0 = P0

# =============================================================================
# # # # run the model
# =============================================================================
print('Start time : {}\n'
      .format((datetime.datetime.now() + datetime.timedelta(0, 1)).time()))

with tqdm(range(int(Y))) as pbar:
    for i in pbar:
        # # update hydrodynamic model
        pbar.set_postfix(inner_loop='update Delft3D')
        modelDIMR.update(mtpervt)
        
        # # extract variables from DFM via BMI
        pbar.set_postfix(inner_loop='extract variables')
        # flow characteristics
        is_sumvalsnd = modelFM.get_var('is_sumvalsnd')
        is_maxvalsnd = modelFM.get_var('is_maxvalsnd')
        Uwave = modelFM.get_var('Uorb')[range(ndxi)]
        Twave = modelFM.get_var('twav')[range(ndxi)]
        # param[range(ndxi), i]
        # > i = 0 : shear stress   [tau]
        # > i = 1 : flow velocity  [vel]
        # > i = 2 : water depth    [wd]
        # morphological characteristics
        diaveg = modelFM.get_var('diaveg')
        dcmodel = diaveg[range(ndxi)]
        stemheight = modelFM.get_var('stemheight')
        rnveg = modelFM.get_var('rnveg')
        
        # # calculate (mean) values from DFM data
        Ucurr = is_sumvalsnd[range(ndxi), 1] / mtpervt
        depth = is_sumvalsnd[range(ndxi), 2] / mtpervt
        
        # # dimensions
        spacetime = np.array(
            [ndxi, len(env.dates[env.dates.dt.year == years[i]])])
        
        # depth = 10. * np.ones(ndxi)
        # Ucurr = .5 * np.ones(ndxi)
        # Uwave = .5 * np.ones(ndxi)
        # Twave = 4. * np.ones(ndxi)
        
        # # environment
        pbar.set_postfix(inner_loop='coral environment')
        # light micro-environment
        lme = Light(
            env.light[env.light.index.year == years[i]].values[:, 0],
            env.Kd[env.Kd.index.year == years[i]].values[:, 0],
            depth)
        lme.repLight(coral)
        # flow micro-environment
        fme = Flow(Ucurr, Uwave, depth, Twave)
        fme.waveCurrent(coral, incanopy=True)
        fme.TBL(coral)
        # thermal micro-environment
        tme = Temperature(
            env.tempK[env.tempK.index.year == years[i]].values[:, 0])
        tme.coralTemperature(coral)
        
        # # physiology
        pbar.set_postfix(inner_loop='coral physiology')
        # photosynthetic dependencies
        if i == 0:
            phd = Photosynthesis(
                env.light[env.light.index.year == years[i]].values[:, 0],
                firstYear=True, TME=False, Ikmax=Ikmax)
        else:
            phd = Photosynthesis(
                env.light[env.light.index.year == years[i]].values[:, 0],
                firstYear=False, TME=False, Ikmax=Ikmax)
        phd.photosyn(coral, env, years[i])
        # population states
        ps = PopStates(K)
        ps.popstates_t(coral)
        # calcification
        cr = Calcification()
        cr.calRate(
            coral, env.acid[env.acid.index.year == years[i]].values[:, 0])
        
        # # morphology
        pbar.set_postfix(inner_loop='coral morphology')
        # morphological development
        mor = MorDevelopment(
            coral.calc.sum(axis=1), depth,
            env.Kd[env.Kd.index.year == years[i]].values[:, 0])
        mor.update(
            coral, env.light[env.light.index.year == years[i]].values[:, 0])
        
        # # storm damage
        if env.stormcat[env.stormcat.index == years[i]].values > 0:
            # return coral data to hydrodynamic model
            pbar.set_postfix(inner_loop='storm - update morphology in Delft3D')
            # translate model parameters
            rnveg = coral.morph2vegden()
            diaveg = coral.dcRep
            stemheight = coral.hc
            # reset counters
            is_sumvalsnd.fill(0.)
            is_maxvalsnd.fill(0.)
            # push counters and updated coral field to model
            modelFM.set_var('is_sumvalsnd', is_sumvalsnd)
            modelFM.set_var('is_maxvalsnd', is_maxvalsnd)
            modelFM.set_var('rnveg', rnveg)
            modelFM.set_var('diaveg', diaveg)
            modelFM.set_var('stemheight', stemheight)
            
            # run storm conditions
            pbar.set_postfix(inner_loop='storm - update Delft3D')
            modelDIMR.update(mtpervt_storm)
            
            # extract variables from DFM via BMI
            pbar.set_postfix(inner_loop='storm - extract variables')
            is_sumvalsnd = modelFM.get_var('is_sumvalsnd')
            is_maxvalsnd = modelFM.get_var('is_maxvalsnd')
            Uwave = modelFM.get_var('Uorb')[range(ndxi)]
            rnveg = modelFM.get_var('rnveg')
            diaveg = modelFM.get_var('diaveg')
            stemheight = modelFM.get_var('stemheight')
            # maximum flow velocity
            Ucurr = is_maxvalsnd[range(ndxi), 1]
            # Uwave = .8
            # Ucurr = 1.
            
            # storm flow environment
            pbar.set_postfix(inner_loop='storm - dislodgement')
            sfe = Flow(Ucurr, Uwave)
            sfe.waveCurrent(coral, incanopy=False)
            # storm dislodgement criterion
            sdc = Dislodgement(sigmat)
            sdc.update(coral)

        # # recruitment / spawning
        pbar.set_postfix(inner_loop='coral recruitment')
        rec = Recruitment(NoLarvae, probSett, dLarv)
        rec.update(coral)
        
        # # return coral data to hydrodynamic model
        pbar.set_postfix(inner_loop='update morphology in Delft3D')
        # translate model parameters
        rnveg = coral.morph2vegden()
        diaveg = coral.dcRep
        stemheight = coral.hc
        # reset counters
        is_sumvalsnd.fill(0.)
        is_maxvalsnd.fill(0.)
        # push counters and updated coral field to model
        modelFM.set_var('is_sumvalsnd', is_sumvalsnd)
        modelFM.set_var('is_maxvalsnd', is_maxvalsnd)
        modelFM.set_var('rnveg', rnveg)
        modelFM.set_var('diaveg', diaveg)
        modelFM.set_var('stemheight', stemheight)
        
# =============================================================================
#       # # # export model results
# =============================================================================
        # # map-file
        pbar.set_postfix(inner_loop='write file - map')
        files = [U2mfile, T2mfile, PS2mfile, P2mfile, G2mfile, M2mfile]
        if any(files):
            if i == 0:
                mset = Dataset(mapfilef, 'w', format='NETCDF4')
                mset.description = 'Mapped simulation data of the CoralModel'

                # dimensions
                mset.createDimension('time', None)
                mset.createDimension('nmesh2d_face', int(ndxi))

                # variables
                t = mset.createVariable('time', int, ('time',))
                t.long_name = 'year'
                t.units = 'years since 0 B.C.'

                x = mset.createVariable('mesh2d_x', 'f8', ('nmesh2d_face',))
                x.long_name = 'x-coordinate'
                x.units = 'm'

                y = mset.createVariable('mesh2d_y', 'f8', ('nmesh2d_face',))
                y.long_name = 'y-coordinate'
                y.units = 'm'

                if U2mfile:
                    Uset = mset.createVariable('ucm', 'f8',
                                                ('time', 'nmesh2d_face'))
                    Uset.long_name = 'in-canopy flow'
                    Uset.units = 'm s^-1'
                if T2mfile:
                    Tset = mset.createVariable('Tc', 'f8',
                                                ('time', 'nmesh2d_face'))
                    Tset.long_name = 'coral temperature'
                    Tset.units = 'K'

                    Tloset = mset.createVariable('Tlo', 'f8',
                                                  ('time', 'nmesh2d_face'))
                    Tloset.long_name = 'lower thermal limit'
                    Tloset.units = 'K'

                    Thiset = mset.createVariable('Thi', 'f8',
                                                  ('time', 'nmesh2d_face'))
                    Thiset.long_name = 'upper thermal limit'
                    Thiset.units = 'K'
                if PS2mfile:
                    PSset = mset.createVariable('PS', 'f8',
                                                ('time', 'nmesh2d_face'))
                    PSset.long_name = 'annual mean photosynthesis'
                    PSset.units = '-'
                if P2mfile:
                    PTset = mset.createVariable('PT', 'f8',
                                                ('time', 'nmesh2d_face'))
                    PTset.long_name = 'total population'
                    PTset.units = '-'

                    PHset = mset.createVariable('PH', 'f8',
                                                ('time', 'nmesh2d_face'))
                    PHset.long_name = 'healthy population'
                    PHset.units = '-'

                    PRset = mset.createVariable('PR', 'f8',
                                                ('time', 'nmesh2d_face'))
                    PRset.long_name = 'recoverying population'
                    PRset.units = '-'

                    PPset = mset.createVariable('PP', 'f8',
                                                ('time', 'nmesh2d_face'))
                    PPset.long_name = 'pale population'
                    PPset.units = '-'

                    PBset = mset.createVariable('PB', 'f8',
                                                ('time', 'nmesh2d_face'))
                    PBset.long_name = 'bleached population'
                    PBset.units = '-'
                if G2mfile:
                    Gset = mset.createVariable('G', 'f8',
                                                ('time', 'nmesh2d_face'))
                    Gset.long_name = 'calcification'
                    Gset.units = 'kg m^-2 y^-1'
                if M2mfile:
                    DCset = mset.createVariable('dc', 'f8',
                                                ('time', 'nmesh2d_face'))
                    DCset.long_name = 'plate diameter'
                    DCset.units = 'm'

                    HCset = mset.createVariable('hc', 'f8',
                                                ('time', 'nmesh2d_face'))
                    HCset.long_name = 'coral height'
                    HCset.units = 'm'

                    BCset = mset.createVariable('bc', 'f8',
                                                ('time', 'nmesh2d_face'))
                    BCset.long_name = 'base diameter'
                    BCset.units = 'm'

                    TCset = mset.createVariable('tc', 'f8',
                                                ('time', 'nmesh2d_face'))
                    TCset.long_name = 'plate thickness'
                    TCset.units = 'm'

                    ACset = mset.createVariable('ac', 'f8',
                                                ('time', 'nmesh2d_face'))
                    ACset.long_name = 'axial distance'
                    ACset.units = 'm'

                # data
                t[:] = np.array([years[i] - 1, years[i]])
                x[:] = xzw[range(ndxi)]
                y[:] = yzw[range(ndxi)]
                if U2mfile:
                    Uset[:, :] = np.array([np.zeros(ndxi), coral.ucm])
                if T2mfile:
                    Tset[:, :] = np.array([np.zeros(ndxi), coral.temp[:, -1]])
                    if tbl:
                        Tloset[:, :] = np.array([np.zeros(ndxi), phd.Tlo])
                        Thiset[:, :] = np.array([np.zeros(ndxi), phd.Thi])
                    else:
                        Tloset[:, :] = np.array(
                            [np.zeros(ndxi), phd.Tlo * np.ones(ndxi)])
                        Thiset[:, :] = np.array(
                            [np.zeros(ndxi), phd.Thi * np.ones(ndxi)])
                if PS2mfile:
                    PSset[:, :] = np.array(
                        [np.zeros(K.shape), coral.photosyn.mean(axis=1)])
                if P2mfile:
                    PTset[:, :] = np.array(
                        [K, coral.popstates[:, -1, :].sum(axis=1)])
                    PHset[:, :] = np.array(
                        [K, coral.popstates[:, -1, 0]])
                    PRset[:, :] = np.array(
                        [np.zeros(K.shape), coral.popstates[:, -1, 1]])
                    PPset[:, :] = np.array([
                        np.zeros(K.shape), coral.popstates[:, -1, 2]])
                    PBset[:, :] = np.array([
                        np.zeros(K.shape), coral.popstates[:, -1, 3]])
                if G2mfile:
                    Gset[:, :] = np.array(
                        [np.zeros(K.shape), coral.calc.sum(axis=1)])
                if M2mfile:
                    DCset[:, :] = np.array([dc0 * K, coral.dc])
                    HCset[:, :] = np.array([hc0 * K, coral.hc])
                    BCset[:, :] = np.array([bc0 * K, coral.bc])
                    TCset[:, :] = np.array([tc0 * K, coral.tc])
                    ACset[:, :] = np.array([ac0 * K, coral.ac])
            else:
                mset = Dataset(mapfilef, mode='a')
                # append data
                mset['time'][:] = np.append(mset['time'][:], years[i])
                if U2mfile:
                    mset['ucm'][-1, :] = coral.ucm
                if T2mfile:
                    mset['Tc'][-1, :] = coral.temp[:, -1]
                    if tbl:
                        mset['Tlo'][-1, :] = phd.Tlo
                        mset['Thi'][-1, :] = phd.Thi
                    else:
                        mset['Tlo'][-1, :] = phd.Tlo * np.ones(ndxi)
                        mset['Thi'][-1, :] = phd.Thi * np.ones(ndxi)
                if PS2mfile:
                    mset['PS'][-1, :] = coral.photosyn.mean(axis=1)
                if P2mfile:
                    mset['PT'][-1, :] = coral.popstates[:, -1, :].sum(axis=1)
                    mset['PH'][-1, :] = coral.popstates[:, -1, 0]
                    mset['PR'][-1, :] = coral.popstates[:, -1, 1]
                    mset['PP'][-1, :] = coral.popstates[:, -1, 2]
                    mset['PB'][-1, :] = coral.popstates[:, -1, 3]
                if G2mfile:
                    mset['G'][-1, :] = coral.calc.sum(axis=1)
                if M2mfile:
                    mset['dc'][-1, :] = coral.dc
                    mset['hc'][-1, :] = coral.hc
                    mset['bc'][-1, :] = coral.bc
                    mset['tc'][-1, :] = coral.tc
                    mset['ac'][-1, :] = coral.ac
            mset.close()

        # # his-file
        pbar.set_postfix(inner_loop='write file - his')
        files = [U2hfile, T2hfile, PS2hfile, P2hfile, G2hfile, M2hfile]
        if any(files):
            date0 = datetime.datetime(2000, 1, 1)
            if i == 0:
                hset = Dataset(hisfilef, 'w', format='NETCDF4')
                hset.description = 'Historic simulation data of the CoralModel'

                # dimensions
                hset.createDimension('time', None)
                hset.createDimension('stations', len(xyn))

                # variables
                t = hset.createVariable('time', 'f8', ('time',))
                t.long_name = 'days since January 1, 2000'
                t.units = 'days'

                x = hset.createVariable('station_x_coordinate', 'f8',
                                        ('stations',))
                y = hset.createVariable('station_y_coordinate', 'f8',
                                        ('stations',))
                if U2hfile:
                    Uset = hset.createVariable('ucm', 'f8',
                                                ('time', 'stations'))
                    Uset.long_name = 'in-canopy flow'
                    Uset.units = 'm s^-1'
                if T2hfile:
                    Tset = hset.createVariable('Tc', 'f8',
                                                ('time', 'stations'))
                    Tset.long_name = 'coral temperature'
                    Tset.units = 'K'

                    Tloset = hset.createVariable('Tlo', 'f8',
                                                  ('time', 'stations'))
                    Tloset.long_name = 'lower thermal limit'
                    Tloset.units = 'K'

                    Thiset = hset.createVariable('Thi', 'f8',
                                                  ('time', 'stations'))
                    Thiset.long_name = 'upper thermal limit'
                    Thiset.units = 'K'
                if PS2hfile:
                    PSset = hset.createVariable('PS', 'f8',
                                                ('time', 'stations'))
                    PSset.long_name = 'photosynthesis'
                    PSset.units = '-'
                if P2hfile:
                    PTset = hset.createVariable('PT', 'f8',
                                                ('time', 'stations'))
                    PTset.long_name = 'total population'
                    PTset.units = '-'

                    PHset = hset.createVariable('PH', 'f8',
                                                ('time', 'stations'))
                    PHset.long_name = 'healthy population'
                    PHset.units = '-'

                    PRset = hset.createVariable('PR', 'f8',
                                                ('time', 'stations'))
                    PRset.long_name = 'recovering population'
                    PRset.units = '-'

                    PPset = hset.createVariable('PP', 'f8',
                                                ('time', 'stations'))
                    PPset.long_name = 'pale population'
                    PPset.units = '-'

                    PBset = hset.createVariable('PB', 'f8',
                                                ('time', 'stations'))
                    PBset.long_name = 'bleached population'
                    PBset.units = '-'
                if G2hfile:
                    Gset = hset.createVariable('G', 'f8',
                                                ('time', 'stations'))
                    Gset.long_name = 'calcification'
                    Gset.units = 'kg m^-2 d^-1'
                if M2hfile:
                    DCset = hset.createVariable('dc', 'f8',
                                                ('time', 'stations'))
                    DCset.long_name = 'plate diameter'
                    DCset.units = 'm'

                    HCset = hset.createVariable('hc', 'f8',
                                                ('time', 'stations'))
                    HCset.long_name = 'coral height'
                    HCset.units = 'm'

                    BCset = hset.createVariable('bc', 'f8',
                                                ('time', 'stations'))
                    BCset.long_name = 'base diameter'
                    BCset.units = '-'

                    TCset = hset.createVariable('tc', 'f8',
                                                ('time', 'stations'))
                    TCset.long_name = 'plate thickness'
                    TCset.units = 'm'

                    ACset = hset.createVariable('ac', 'f8',
                                                ('time', 'stations'))
                    ACset.long_name = 'axial distance'
                    ACset.units = 'm'

                # data indices
                xs = xyn['x'].values
                ys = xyn['y'].values
                idx = np.zeros(len(xs))
                for s in range(len(xs)):
                    idx[s] = ((xzw - xs[s]) ** 2 + (yzw - ys[s]) ** 2).argmin()
                idx = idx.astype(int)

                # data
                idates = env.dates[
                        env.dates.dt.year == years[i]].reset_index(drop=True)
                t[:] = (idates - date0).dt.days.values
                x[:] = xs
                y[:] = ys
                if U2hfile:
                    Uset[:, :] = np.tile(coral.ucm, (len(idates), 1))[:, idx]
                if T2hfile:
                    Tset[:, :] = coral.temp[idx, :].transpose()
                    if tbl:
                        Tloset[:, :] = np.tile(
                            phd.Tlo, (len(idates), 1))[:, idx]
                        Thiset[:, :] = np.tile(
                            phd.Thi, (len(idates), 1))[:, idx]
                    else:
                        Tloset[:, :] = phd.Tlo * np.ones(
                            (len(idates), len(idx)))
                        Thiset[:, :] = phd.Thi * np.ones(
                            (len(idates), len(idx)))
                if PS2hfile:
                    PSset[:, :] = coral.photosyn[idx, :].transpose()
                if P2hfile:
                    PTset[:, :] = coral.popstates[idx, :, :].sum(
                        axis=2).transpose()
                    PHset[:, :] = coral.popstates[idx, :, 0].transpose()
                    PRset[:, :] = coral.popstates[idx, :, 1].transpose()
                    PPset[:, :] = coral.popstates[idx, :, 2].transpose()
                    PBset[:, :] = coral.popstates[idx, :, 3].transpose()
                if G2hfile:
                    Gset[:, :] = coral.calc[idx, :].transpose()
                if M2hfile:
                    DCset[:, :] = np.tile(coral.dc, (len(idates), 1))[:, idx]
                    HCset[:, :] = np.tile(coral.hc, (len(idates), 1))[:, idx]
                    BCset[:, :] = np.tile(coral.bc, (len(idates), 1))[:, idx]
                    TCset[:, :] = np.tile(coral.tc, (len(idates), 1))[:, idx]
                    ACset[:, :] = np.tile(coral.ac, (len(idates), 1))[:, idx]
            else:
                hset = Dataset(hisfilef, mode='a')
                # date conversion
                idates = env.dates[
                        env.dates.dt.year == years[i]].reset_index(drop=True)
                t = (idates - date0).dt.days.values
                # append data
                hset['time'][:] = np.append(hset['time'][:], t)
                if U2hfile:
                    hset['ucm'][t, :] = np.tile(
                        coral.ucm, (len(idates), 1))[:, idx]
                if T2hfile:
                    hset['Tc'][t, :] = coral.temp[idx, :].transpose()
                    if tbl:
                        hset['Tlo'][t, :] = np.tile(
                                phd.Tlo, (len(idates), 1))[:, idx]
                        hset['Thi'][t, :] = np.tile(
                                phd.Thi, (len(idates), 1))[:, idx]
                    else:
                        hset['Tlo'][t, :] = phd.Tlo * np.ones(
                                (len(idates), len(idx)))
                        hset['Thi'][t, :] = phd.Thi * np.ones(
                                (len(idates), len(idx)))
                if PS2hfile:
                    hset['PS'][t, :] = coral.photosyn[idx, :].transpose()
                if P2hfile:
                    hset['PT'][t, :] = coral.popstates[idx, :, :].sum(
                        axis=2).transpose()
                    hset['PH'][t, :] = coral.popstates[idx, :, 0].transpose()
                    hset['PR'][t, :] = coral.popstates[idx, :, 1].transpose()
                    hset['PP'][t, :] = coral.popstates[idx, :, 2].transpose()
                    hset['PB'][t, :] = coral.popstates[idx, :, 3].transpose()
                if G2hfile:
                    hset['G'][t, :] = coral.calc[idx, :].transpose()
                if M2hfile:
                    hset['dc'][t, :] = np.tile(
                        coral.dc, (len(idates), 1))[:, idx]
                    hset['hc'][t, :] = np.tile(
                        coral.hc, (len(idates), 1))[:, idx]
                    hset['bc'][t, :] = np.tile(
                        coral.bc, (len(idates), 1))[:, idx]
                    hset['tc'][t, :] = np.tile(
                        coral.tc, (len(idates), 1))[:, idx]
                    hset['ac'][t, :] = np.tile(
                        coral.ac, (len(idates), 1))[:, idx]
            hset.close()

# =============================================================================
# # # # finalize the model
# =============================================================================
modelDIMR.finalize()
print('\nModel finalized.')
print('\nEnd time   : {0}'.format(datetime.datetime.now().time()))