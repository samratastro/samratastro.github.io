## packages
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

#### import pamater file for plotting
import sys
plotparam_dir = '/Users/samrat/Dropbox/codes/misc/'
sys.path.append(plotparam_dir)
import params as param

############################ ================================================ Synthetic IMAGING ===========================================



### SolO/HRI_EUV synthesis for 2.5D MHD
gamma = 5/3
normalisations = dict(length_unit=(1e8, 'cm'), temperature_unit=(1e6, 'K'), numberdensity_unit=(1e9, 'cm**-3'))
animationloc = figfileloc + "animation/HRI-EUV/"
### =========  AIA passband response functions ==============
# tabfolder   = "/Users/samratsen/iac/my_codes/synthesis/response_tab/"
tabfile     = tabfolder+"goft_table_eui_hri_174.sav"
restore     = idl.readsav(tabfile)
density     = restore("density")
temperature = restore("temperature")
table       = restore("table")

os.chdir(fileloc)
files = sorted(glob.glob('*.dat'))
def ds_yt(files):    
    ds = yt.load(files, units_override = normalisations, unit_system = 'cgs')
    level = ds.max_level
    Temp = 1e6 * ds.covering_grid(level, left_edge=ds.domain_left_edge,
                                dims = ds.domain_dimensions * ds.refine_by**level)['Te'][:,:,0].value ## Kelvin

    ne = ds.covering_grid(level, left_edge=ds.domain_left_edge,
                                dims = ds.domain_dimensions * ds.refine_by**level)['number_density'][:,:,0].value ## cm^-3
    return ne, Temp

start_frame = 187
end_frame = 187
cell_size = 3.272 * 1e6 ## cm

for i in range(start_frame, end_frame +1):
    ne, Temp = ds_yt(files[i])
    #### ====rearranement of ne, T ============== 
    nx       = ne.shape[0]
    ny       = ne.shape[1]
    tg       = np.log10(Temp)      # log10(Tg [K])
    nel      = np.log10(ne)  # log10(Electron number density [cm^-3])
    tg       = tg.reshape(nx*ny)
    nel      = nel.reshape(nx*ny)
    nh       = 0.85*10**nel

    myabund     = 7.85  # Abundance of Fe (it depends on which CHIANTI abund file you use)
# tabfolder   = "/Users/samratsen/iac/my_codes/synthesis/response_tab/"
    tabfolder   = "/Users/samrat/response_tab/"
    tabfile     = tabfolder+"goft_table_eui_hri_174.sav"
    restore     = idl.readsav(tabfile)
    density     = restore("density")
    temperature = restore("temperature")
    table       = restore("table")       
    Roftn = interpolate.interpn((temperature, density), table,  np.array([tg, nel]).T, 
                           fill_value=np.nan, method='linear', bounds_error=False)  
    goftn = (1/7.3) * Roftn   # ph cm^+5 s^-1 pix^-1 For Solo/HRIEUV 7.29 Dn = 1ph

    emiss = 10**(myabund - 12.0)*10**nel*nh*goftn    # (ph cm^-1 s^-1 pix^-1)
    emiss = emiss.reshape(nx,ny)                     # (ph cm^-1 s^-1 pix^-1)
    emiss = np.nan_to_num(emiss, nan=1e-12)          # Remove nan for really small values
    emiss[emiss == 0.0] = 1e-12  

    #### Resolution Degrading ====================

    def degrade_resolution(array, factor):
        # Trim the array to make dimensions divisible by the factor
        trimmed_shape = (array.shape[0] // factor * factor, array.shape[1] // factor * factor)
        trimmed_array = array[:trimmed_shape[0], :trimmed_shape[1]]
        
        # Reshape and compute the mean over blocks
        return trimmed_array.reshape(trimmed_shape[0] // factor, factor, trimmed_shape[1] // factor, factor).mean(axis=(1, 3))

    # Downsample the emissivity to SolO/HRI resolution
    factor = 4  # HRI-EUV resolution (factor=13 <==> 435 km)
    original = 1     
    emiss_simulation = degrade_resolution(emiss, original)
    emiss_downsampled = degrade_resolution(emiss, factor)
    plt.figure(figsize=(8, 6))
    
    vmin= -9
    vmax= -6
    plt.subplot(1,2,1)
    plt.imshow(np.log10(emiss_simulation[:, 0:740].T), origin='lower', \
            cmap=ct.aia_color_table(171*u.Angstrom) , vmin=vmin, vmax=vmax, extent=[0, 6.28, 0, 24], aspect=0.5)
    plt.title("simulation resolution")
    plt.text(-2, 10, f"t = {(cad * i/60):.2f} min", rotation=90, fontsize=14, color='red', fontweight='bold')
    plt.xlabel('x (Mm)')
    plt.ylabel('y (Mm)')
    plt.subplot(1,2,2)
    plt.imshow(np.log10(emiss_downsampled[:, 0:183].T), origin='lower', \
            cmap=ct.aia_color_table(171*u.Angstrom) , vmin=vmin, vmax=vmax, extent=[0, 6.28, 0, 24], aspect=0.5,\
                 interpolation='gaussian')


    ### emissivity plot ======
    # plt.figure(figsize=(5, 6))
    # plt.imshow(np.log10(emiss[:, 0:740].T), origin='lower', cmap=ct.aia_color_table(171*u.Angstrom), vmin=-9, vmax=-6,\
    #            extent=[0, 6.28, 0, 24], aspect=0.5)
    plt.colorbar(pad=0.005, shrink=0.5, label='log$(\epsilon)$ [ph cm$^{-1}$ s$^{-1}$ pix$^{-1}$]')
    plt.xlabel('x (Mm)')
    plt.title("HRI$_{EUV}$ resolution")
    
   
    plt.tight_layout()
    # plt.savefig(animationloc + f'emiss_EUI-HRI_t{i}.png', dpi=250)
    # plt.clf()
    plt.show()

#### For IRIS-SJI imaging in 2.5D

gamma = 5/3
normalisations = dict(length_unit=(1e8, 'cm'), temperature_unit=(1e6, 'K'), numberdensity_unit=(1e9, 'cm**-3'))
animationloc = figfileloc + "animation/IRIS-SJI/"
     
os.chdir(fileloc)
files = sorted(glob.glob('*.dat'))
def ds_yt(files):    
    ds = yt.load(files, units_override = normalisations, unit_system = 'cgs')
    level = ds.max_level
    Temp = 1e6 * ds.covering_grid(level, left_edge=ds.domain_left_edge,
                                dims = ds.domain_dimensions * ds.refine_by**level)['Te'][:,:,0].value ## Kelvin

    ne = ds.covering_grid(level, left_edge=ds.domain_left_edge,
                                dims = ds.domain_dimensions * ds.refine_by**level)['number_density'][:,:,0].value ## cm^-3
    return ne, Temp

tabfolder   = "/Users/samrat/response_tab/"
got_file1 = tabfolder + 'goft_table_isji140_2014-04-24_abph_extro.dat'

start_frame = 170
end_frame = 195
cell_size = 3.272 * 1e6 ## cm
    
for kk in range(start_frame, end_frame +1):
    ne, Temp = ds_yt(files[kk])
    nx       = ne.shape[0]
    ny       = ne.shape[1]
    # tg and nel should be in log10 for the interpolation with the table
    tg       = np.log10(Temp)      # log10(Tg [K])
    nel      = np.log10(ne)  # log10(Electron number density [cm^-3])
    # We reshape to 1D because interpolations are way faster like this
    tg          = tg.reshape(nx*ny)
    nel         = nel.reshape(nx*ny)
    nh          = 0.85*10**nel

    with open(got_file1) as f:
        lines1 = [line.rstrip() for line in f]

    num = lines1[7].strip()
    numn = int(num[0:3])
    numt = int(num[6:9])
    goft_mat1 = np.zeros((numn, numt))
    n_e_lg = np.zeros((numn))
    n_e_lg[0] = float(lines1[48].strip())
    logt_sji1 = np.zeros((numt))
    for i in range(0,40):
        logt_sji1[i*5:i*5+5] = [float(s) for s in lines1[8+i].split()]

    goft_mat1 = np.zeros((numn,numt))
    n_e_lg = np.zeros((numn))
    n_e_lg[0] = float(lines1[48].strip())
    for j in range(0, 6):
        n_e_lg[j] = float(lines1[40*(j+1)+8+j].strip())
        for i in range(0,40):
            goft_mat1[j,i*5:i*5+5] = [float(s) for s in lines1[40*(j+1)+8+j+i+1].split()]

    temperature = logt_sji1
    table = goft_mat1[2,:]  ## 2 corresponds to ne = 1e9 cm^-3

    f_interp = interpolate.interp1d(temperature, table, kind='linear', bounds_error=False, fill_value=np.nan)
    Roft    = f_interp(tg)   ## in DN cm^+5 s^-1 sr^-1
    goft = 12 * (6.5*1e-13) * Roft ## in ph cm^+5 s^-1 pix^-1 :: For IRIS-SJI 1 DN = 12 ph
    emiss = 10**nel*nh*goft                          
    emiss = emiss.reshape(nx,ny)                     
    emiss = np.nan_to_num(emiss, nan=1e-12)          # Remove nan for really small values
    emiss[emiss == 0.0] = 1e-12                    # To avoid problems with the logarithm in the plot

    #### ======== Resolution Degrading ====================

    def degrade_to_instrument_resolution(array, model_pixel_km=32.6, instrument_fwhm_km=240):
        fwhm_pixels = instrument_fwhm_km / model_pixel_km
        sigma = fwhm_pixels / 2.355

        blurred = gaussian_filter(array, sigma=sigma)

        # Downsample to instrument pixel scale
        factor = int(round(instrument_fwhm_km / model_pixel_km))  # ≈ 7.3 for IRIS-SJI

        return blurred[::factor, ::factor]

    # Downsample the emissivity to SolO/HRIEUV resolution with Gaussian Convolution    
    emiss_simulation = emiss.copy()
    emiss_downsampled = degrade_to_instrument_resolution(emiss)


    ### emissivity plot ======
    vmin= -10.2   # # 1.5 * (1/6) * (6.5*1e-13)
    vmax= -6.5 # 9.4     # 3.5 * (1/6) * (6.5*1e-13)

    plt.figure(figsize=(8, 6))
    plt.subplot(1,2,1)
    plt.imshow(np.log10(emiss_simulation[:, 0:740].T), vmin=vmin, vmax=vmax, \
               origin='lower', cmap='irissji1400', extent=[0, 6.28, 0, 24], aspect=0.5)
    plt.title("simulation resolution")
    plt.text(-2, 10, f"t = {(cad * kk/60):.2f} min", rotation=90, fontsize=14, color='red', fontweight='bold')
    plt.xlabel('x (Mm)')
    plt.ylabel('y (Mm)')

    plt.subplot(1,2,2)
    plt.imshow(np.log10(emiss_downsampled[:, 0:105].T), vmin=vmin, vmax=vmax,\
               origin='lower', cmap='irissji1400', extent=[0, 6.28, 0, 24], aspect=0.5, interpolation='gaussian')
    plt.xlabel('x (Mm)')
    plt.title("IRIS-SJI resolution")
    plt.colorbar(pad=0.005, shrink=0.5, label='log($\epsilon$) [ph cm$^{-1}$ s$^{-1}$ pix$^{-1}$]')
    plt.tight_layout()
    plt.savefig(animationloc + f'IRIS-SJI_t{kk}.png', dpi=250)
    plt.clf()
    # plt.show()









