########## Import packages =======================
import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd
print(yt.__version__)

from mpl_toolkits.axes_grid1 import AxesGrid
from yt.funcs import mylog as ytlog
ytlog.setLevel(40)
import matplotlib as mpl


# from mpl_toolkits import mplot3d
# from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt 
# from matplotlib import animation
import scipy.integrate as sint
import scipy.constants as const
from scipy.ndimage import uniform_filter
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import operator
import functools
import glob
import subprocess
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import matplotlib
import math

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# Sunpy
import sunpy.visualization.colormaps.color_tables as ct
import sunpy.visualization.colormaps
import astropy.units as u

# Scipy
import scipy.io.idl as idl
from scipy import interpolate

# Others
import os
import warnings
warnings.simplefilter("ignore")
### import param file =======
import sys
plotparam_dir = '/Users/samrat/Dropbox/codes/misc/'
sys.path.append(plotparam_dir)
import params as param

##########===================================================
### ================= MAIN CODE =============================
##########===================================================

## ======= Single line ===============
normalisations = dict(length_unit=(1e8, 'cm'), temperature_unit=(1e6, 'K'), numberdensity_unit=(1e9, 'cm**-3'))

fileloc = '/Users/samrat/simulations/shear_arcade/synthesis/' ## laptop
figfileloc = '/Users/samrat/simulations/shear_arcade/PFR-synthesis/figures/' ## Laptop
# figfileloc = '/Users/samratsen/iac/my_codes/synthesis/post_flare_rain/figures/'

gamma = 5/3
normalisations = dict(length_unit=(1e8, 'cm'), temperature_unit=(1e6, 'K'), numberdensity_unit=(1e9, 'cm**-3'))

ds = yt.load(fileloc + "shear_arcade_0187.dat", units_override = normalisations, unit_system = 'cgs')


level = ds.max_level
Temp = 1e6 * ds.covering_grid(level, left_edge=ds.domain_left_edge,
                            dims = ds.domain_dimensions * ds.refine_by**level)['Te'][:,:,0].value ## Kelvin
ne = ds.covering_grid(level, left_edge=ds.domain_left_edge,
                            dims = ds.domain_dimensions * ds.refine_by**level)['number_density'][:,:,0].value ## cm^-3

cad = ds.time_unit.value / 2 # cadence in seconds

print('time = ', ds.current_time.in_units('min'))



#--------------------------------------------------------------------------------------------------------------#
# Load simulation variables (tg, nel, nx, ny, nz)
#--------------------------------------------------------------------------------------------------------------#
nx       = ne.shape[0]
ny       = ne.shape[1]
# tg and nel should be in log10 for the interpolation with the table
tg       = np.log10(Temp)      # log10(Tg [K])
nel      = np.log10(ne)  # log10(Electron number density [cm^-3])
# We reshape to 1D because interpolations are way faster like this
tg          = tg.reshape(nx*ny)
nel         = nel.reshape(nx*ny)
nh          = 0.85*10**nel

#--------------------------------------------------------------------------------------------------------------#
# ================ coronal lines ===========================
#--------------------------------------------------------------------------------------------------------------#
# myabund     = 7.85  # Abundance of Fe (it depends on which CHIANTI abund file you use)
myabund     = 7.51  # Abundance of Si (it depends on which CHIANTI abund file you use)
# tabfolder   = "/Users/samratsen/iac/my_codes/synthesis/response_tab/"
tabfolder   = "/Users/samrat/response_tab/"

tabfile     = tabfolder+"goft_table_si_4_1402.770.sav"

restore     = idl.readsav(tabfile)
density     = restore("density")
temperature = restore("temperature")
table       = restore("table")

#--------------------------------------------------------------------------------------------------------------#
# Emiss
#--------------------------------------------------------------------------------------------------------------#
goftn = interpolate.interpn((temperature, density), table,  np.array([tg, nel]).T, 
                           fill_value=np.nan, method='linear', bounds_error=False)

emiss = 10**(myabund - 12.0)*10**nel*nh*goftn    # (erg cm-3 sr-1 s-1)
emiss = emiss.reshape(nx,ny)                     # (erg cm-3 sr-1 s-1)
emiss = np.nan_to_num(emiss, nan=1e-12)          # Remove nan for really small values
emiss[emiss == 0.0] = 1e-12                      # To avoid problems with the logarithm in the plot






###================== Spectral line (single slit) ======================
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.interpolate import RegularGridInterpolator

kb = 1.3807 * 1e-16 ## cm^2 g s^-2 K^-1
# mFe = 9.2732 * 1e-23 ## g  ## Iron
mFe = 4.6495 * 1e-23 ## g  ## Silicon
c = 2.9979 * 1e10 ## cm/s

########   central wavelengths for different lines =================
# lam0 = 93.932 * 1e-8 ## cm
# lam0 = 108.355 * 1e-8 ## cm
# lam0 = 171.0730 * 1e-8 ## cm 
# lam0 = 174.5310 * 1e-8 ## cm
# lam0 = 193.509 * 1e-8 ## cm
# lam0 = 211.3170 * 1e-8 ## cm 
# lam0 = 284.163 * 1e-8 ## cm
# lam0 = 335.410 * 1e-8 ## cm
lam0 = 1402.770 * 1e-8 ## cm
#############    =========================== ####################################
width = 1.5 * 1e-8           # 1.5 * 1e-8 #  
n_lambda = 115  # 300 ### 
nxarray= np.arange(emiss.shape[0])
nyarray= np.arange(emiss.shape[1]) 

deltalam = (lam0/c) * np.sqrt(2 * kb * Temp / mFe) ## cm
lamshift = lam0 * (1 + vy / c) ## cm
cell_size = 3.26 * 1e6 ## cm

lamarray = np.linspace(lam0-width, lam0+width, n_lambda) ## cm

I_lambda = np.zeros((len(nxarray), len(nyarray), n_lambda))

for k in range(n_lambda):
    I_lambda[:, :, k] = emiss * cell_size * (1 / (deltalam * np.sqrt(np.pi))) * np.exp(
        -((lamarray[k] - lamshift) / deltalam) ** 2
    )

#### ===============Integration along y for fixed nx ===================
def integrate_along_y(nx, nyarray, lamarray):
    integrated_intensity = np.zeros(len(lamarray))
    for k in range(len(lamarray)):
        integrated_intensity[k] = np.sum(I_lambda[nx, :, k])
    return integrated_intensity
fixednx = 96
nyarray = np.arange(Temp.shape[1])
Intensity_y = integrate_along_y(fixednx, nyarray, lamarray)

#### ============== plotting =================================================
# plt.figure(figsize=(8, 4))
# plt.plot((np.array(lamarray) - lam0) * 1e8, Intensity_y / np.max(Intensity_y), label='Intensity', color='black', lw=2)
# # plt.axvline(x=1e8 * lam0, color='r', linestyle='--', lw=1, label='Fe XVI 211.3170 A')
# plt.xlabel('$\Delta \lambda$ ($\AA$)')
# plt.ylabel('Intensity')
# plt.legend()
# # plt.title('Integrated Intensity along Y direction for fixed X')
print('Spectral Resolution = ', 1e11 * 2 * width / n_lambda, 'mA')
# print('time = ', ds.current_time.in_units('min'))
# plt.show()

vdoparray = 1e-5 * (c/lam0) * (np.array(lamarray) - lam0)  ## km/s
Inorm = Intensity_y / np.max(Intensity_y)
plt.plot(-vdoparray, Inorm, color='black', lw=2)
# plt.yscale('log')
plt.xlabel('$\Delta v$ (km/s)')
plt.ylabel('Intensity')
plt.xlim(-320, 20)
print('time = ', ds.current_time.in_units('min')) 
plt.show()


###########=======================================================================
################################### Spectral Map =================================
###########=======================================================================
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.interpolate import RegularGridInterpolator

kb = 1.3807 * 1e-16 ## cm^2 g s^-2 K^-1
# mFe = 9.2732 * 1e-23 ## g  ## Iron
mFe = 4.6495 * 1e-23 ## g  ## Silicon
c = 2.9979 * 1e10 ## cm/s

lam0 = 1402.770 * 1e-8 ## cm  ########   central wavelength
width = 1.5 * 1e-8 #  3.2 * 1e-8 ## 0.65 * 1e-8 ## cm
n_lambda = 300 ## 300 # 115
nxarray= np.arange(emiss.shape[0])
nyarray= np.arange(emiss.shape[1]) 

deltalam = (lam0/c) * np.sqrt(2 * kb * Temp / mFe) ## cm
lamshift = lam0 * (1 + vy / c) ## cm
cell_size = 3.26 * 1e6 ## cm

lamarray = np.linspace(lam0-width, lam0+width, n_lambda) ## cm

I_lambda = np.zeros((len(nxarray), len(nyarray), n_lambda))

for k in range(n_lambda):
    I_lambda[:, :, k] = emiss * cell_size * (1 / (deltalam * np.sqrt(np.pi))) * np.exp(
        -((lamarray[k] - lamshift) / deltalam) ** 2
    )

def integrate_along_y(nx, nyarray, lamarray):
    integrated_intensity = np.zeros(len(lamarray))
    for k in range(len(lamarray)):
        integrated_intensity[k] = np.array(np.sum(I_lambda[nx, :, k]))
    return integrated_intensity
xpixarr = np.arange(96-20, 96+21)
nyarray = np.arange(Temp.shape[1])

Intensity_y_spmap = []
for nxx in range(0, len(xpixarr)):
    Intensity_y_spmap.append(integrate_along_y(xpixarr[nxx], nyarray, lamarray))


