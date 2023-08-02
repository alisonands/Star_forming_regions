#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 12:09:43 2023

@author: alisonandrade
"""

import astropy
import matplotlib
import aplpy
from astropy.visualization import astropy_mpl_style
import math as m
import scipy.optimize as opt
import statistics as st

from regions import Regions
from regions.core import PixCoord
from regions.shapes import RectangleSkyRegion, RectanglePixelRegion

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import Angle, SkyCoord
from astropy.utils.data import get_pkg_data_filename
from astropy import units as u
import numpy as np

matplotlib.rcParams.update({'font.size': 11})

#%%


#%% Imports/reading

regions_file = '/Users/alisonandrade/Documents/alison_17/590thesis/orion/OrionMaps/regions.reg'
orion_A_file = '/Users/alisonandrade/Documents/alison_17/590thesis/orion/OrionMaps/Planck_353GHz_2048_l212.0_b-19.0_w8_fwhm10.fits'
orion_A_cov_file = '/Users/alisonandrade/Documents/alison_17/590thesis/orion/OrionMaps/Planck_353GHz_2048_l212.0_b-19.0_w8_fwhm10_cov.fits'
#orion_A_data = get_pkg_data_filename(orion_A_file)
orion_A = fits.open(orion_A_file)
orion_A_cov = fits.open(orion_A_cov_file)
orion_A_cov_data = orion_A_cov[0].data
orion_A.info()
orion_A_cov.info()

#%%

with fits.open(orion_A_file) as orion_A:
    i_stokes = orion_A[0].data
    q_stokes = orion_A[1].data
    u_stokes = orion_A[2].data
    polint = orion_A[3].data
    pol_frac = orion_A[4].data
    bpos_ang = orion_A[5].data
    pol_disp = orion_A[6].data
    
with fits.open(orion_A_cov_file) as orion_A_cov:
    ii_stokes_cov = orion_A_cov[0].data
    iq_stokes_cov = orion_A_cov[1].data
    iu_stokes_cov = orion_A_cov[2].data
    qq_stokes_cov = orion_A_cov[3].data
    qu_stokes_cov = orion_A_cov[4].data
    uu_stokes_cov = orion_A_cov[5].data
    pi_var_cov = orion_A_cov[6].data
    pf_var_cov = orion_A_cov[7].data  ##<---  pol frac??
    bang_var_cov = orion_A_cov[8].data
    s_var_cov = orion_A_cov[9].data
    
regions = Regions.read(regions_file, format = 'ds9') ##same as parsing?

mask_cov = np.logical_and(pol_frac, polint>3*np.sqrt(pi_var_cov))
    
wcs = WCS(orion_A[4].header)
plt.imshow(pol_frac)
plt.imshow(mask_cov)
#plt.imshow(i_stokes)

#%%

orion_A = fits.open(orion_A_file)

plt.imshow(orion_A[0].data)
plt.title('I Stokes Orion A')

regions_pix = []
for x in range(len(regions)):
    region0 = (RectangleSkyRegion(regions[x].center, regions[x].width, regions[x].height, regions[x].angle))
    regions_pix.append(region0.to_pixel(wcs))
    regions_pix[x].plot()

#%%

fwhm_10 = fits.getdata(orion_A_file) ##??
plt.imshow(pol_frac)
plt.title('Polarized fraction Orion A (mask)')

ax = plt.subplot(1, 1, 1)
mask_list = []

for x in range(len(regions_pix)):
    mask = regions_pix[x].to_mask()
    mask_list.append(mask)
    ax.add_artist(mask.bbox.as_artist())

print(np.shape(pol_frac))
print(np.shape(polint))

#%%


ax = plt.subplot(1, 1, 1)
mask0 = mask_list[0].to_image(pol_frac.shape)  

mask = np.logical_and(pol_frac[0], mask0 == 1)
mask_coord = mask.nonzero()                     #loc in list where value of mask is True

pol_frac[mask_coord]

pol_frac_list = []

for x in range(len(mask_list)):
    mask0 = mask_list[x].to_image(pol_frac.shape)
    
    mask_trial = np.logical_and(pol_frac, mask0 == 1)
    mask_coord = mask_trial.nonzero() 
    pol_frac_list.append(pol_frac[mask_coord])
    pol_frac_list[x] = pol_frac_list[x][~np.isnan(pol_frac_list[x])] ##removing all nan values (i.e. from list [3])

    plt.imshow(mask_trial)
    plt.title('mask #' + str(x))
    plt.show()
    
hist, bins, patches = plt.hist(pol_frac_list[3], 100)
for x in range(len(pol_frac_list)):
    hist, bins, patches = plt.hist(pol_frac_list[x], 100)
    plt.title('pol frac reg wo/cov '+ str(x))
    plt.show()

p_max = 0.15
p_0 = (3*p_max)/(3+p_max)

for x in range(len(pol_frac_list)):
    gamma = [(1/2)*np.arccos(((2*p_obs*(1+(2/3)*p_0))/(p_0*(1+p_obs)))-1)*180/np.pi for p_obs in pol_frac_list[x]]

    hist, bins, patches = plt.hist(gamma, 100)
    plt.title('gamma reg wo/cov '+ str(x))
    plt.show()

#%%

ax = plt.subplot(1, 1, 1)
mask0 = mask_list[0].to_image(pol_frac.shape)  

mask = np.logical_and(pol_frac[0], mask0 == 1)
mask_coord = mask.nonzero()                     #loc in list where value of mask is True

pol_frac[mask_coord]

pol_frac_list = []

for x in range(len(mask_list)):
    mask0 = mask_list[x].to_image(pol_frac.shape)
    
    mask_trial = np.logical_and(mask_cov, mask0 == 1)
    mask_coord = mask_trial.nonzero() 
    pol_frac_list.append(pol_frac[mask_coord])
    pol_frac_list[x] = pol_frac_list[x][~np.isnan(pol_frac_list[x])] ##removing all nan values (i.e. from list [3])

    plt.imshow(mask_trial)
    plt.title('mask #' + str(x))
    plt.show()
    
hist, bins, patches = plt.hist(pol_frac_list[3], 100)
for x in range(len(pol_frac_list)):
    hist, bins, patches = plt.hist(pol_frac_list[x], 100)
    plt.title('pol frac reg w/cov '+ str(x))
    plt.show()

p_max = 0.15
p_0 = (3*p_max)/(3+p_max)

for x in range(len(pol_frac_list)):
    gamma = [(1/2)*np.arccos(((2*p_obs*(1+(2/3)*p_0))/(p_0*(1+p_obs)))-1)*180/np.pi for p_obs in pol_frac_list[x]]

    hist, bins, patches = plt.hist(gamma, 100)
    plt.title('gamma reg w/cov '+ str(x))
    plt.show()
    



#%%
