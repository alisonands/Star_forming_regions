#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:41:27 2023

@author: alisonandrade
"""

import matplotlib.pyplot as plt
import numpy as np

from regions import Regions
from astropy.io import fits
from astropy.wcs import WCS
import fits_file_handling as fts
from regions.shapes import RectangleSkyRegion, RectanglePixelRegion


#%%

regions_file = '/Users/alisonandrade/Documents/alison_17/590thesis/orion/OrionMaps/regions.reg'
regions = Regions.read(regions_file, format = 'ds9') ##same as parsing?

col_dens_file = '/Users/alisonandrade/Documents/alison_17/590thesis/orion/OrionMaps/Planck_353GHz_2048_l212.0_b-19.0_w8_fwhm6_dust.fits'

pol_frac_file = '/Users/alisonandrade/Documents/alison_17/590thesis/orion/OrionMaps/Orion_ja_columndensity.fits'

orion_A_file = '/Users/alisonandrade/Documents/alison_17/590thesis/orion/OrionMaps/Planck_353GHz_2048_l212.0_b-19.0_w8_fwhm10.fits'
with fits.open(orion_A_file) as orion_A:
    i_stokes = orion_A[0].data
    q_stokes = orion_A[1].data
    u_stokes = orion_A[2].data
    polint = orion_A[3].data
    pol_frac = orion_A[4].data
    bpos_ang = orion_A[5].data
    pol_disp = orion_A[6].data

with fits.open(col_dens_file) as col_dens:
    tau353 = col_dens[0].data
    
wcs = WCS(orion_A[0].header)

#tau353 = tau353_0*8.3333*10**25
#%%
plt.imshow(polint)

regions_pix = []
for x in range(len(regions)):
    region0 = (RectangleSkyRegion(regions[x].center, regions[x].width, regions[x].height, regions[x].angle))
    regions_pix.append(region0.to_pixel(wcs))
    regions_pix[x].plot()

#%%
#mask_list_pf = fts.mask_box(polint, regions_pix, wcs)
plt.imshow(pol_frac)
ax = plt.subplot(1, 1, 1)
mask_list_pf = []

for x in range(len(regions_pix)):
    mask = regions_pix[x].to_mask()
    mask_list_pf.append(mask)
    ax.add_artist(mask.bbox.as_artist())

#%%
pol_frac_list_pf = []
for x in range(len(pol_frac)):
    print(pol_frac[x])

for x in range(len(mask_list_pf)):
    mask0 = mask_list_pf[x].to_image(pol_frac.shape)
    
    mask_trial = np.logical_and(pol_frac, mask0 == 1)
    mask_coord = mask_trial.nonzero() 
    pol_frac_list_pf.append(pol_frac[mask_coord])

pol_frac_list_pf

#%%
mask_list_int = fts.mask_box(polint, regions_pix, wcs)

#%% make regions for intensity

pol_frac_list_int = []

for x in range(len(mask_list_int)):
    mask0 = mask_list_int[x].to_image(polint.shape)
    
    mask_trial = np.logical_and(polint, mask0 == 1)
    mask_coord = mask_trial.nonzero() 
    pol_frac_list_int.append(polint[mask_coord])

mask_list = []

for x in range(len(regions_pix)):
    mask = regions_pix[x].to_mask()
    mask_list.append(mask)
    #ax.add_artist(mask.bbox.as_artist())

tau353_list = []

for x in range(len(mask_list)):
    mask0 = mask_list[x].to_image(tau353.shape)
    
    mask_trial = np.logical_and(tau353, mask0 == 1)
    mask_coord = mask_trial.nonzero() 
    tau353_list.append(tau353[mask_coord])
    #tau353_list[x] = tau353_list[x][~np.isnan(tau353[x])] ##removing all nan values (i.e. from list [3])
    
    plt.imshow(mask_trial)
    plt.title('mask #' + str(x))
    plt.show()

for x in range(len(tau353_list)):
    hist, bins, patches = plt.hist(tau353_list[x], 100)
    plt.title('col density wo/cov '+ str(x))
    plt.show()

#%%

for x in range(3):
    plt.hist2d(tau353_list[x], pol_frac_list_pf[x], bins=(100, 100), cmap=plt.cm.jet)
    plt.title(x)
    plt.xlabel('pol fraction')
    plt.ylabel('pol intensity')
    plt.show()
    
for x in range(4,6):
    plt.hist2d(tau353_list[x], pol_frac_list_pf[x], bins=(100, 100), cmap=plt.cm.jet)
    plt.title(x)
    plt.xlabel('pol fraction')
    plt.ylabel('pol intensity')
    plt.show()

#plt.hist2d(pol_frac_list_pf[5], pol_frac_list_int[5], bins=(100, 100), cmap=plt.cm.jet)
#plt.xlabel('pol fraction')
#plt.ylabel('pol intensity')
#plt.title('5')
#plt.show()
    
