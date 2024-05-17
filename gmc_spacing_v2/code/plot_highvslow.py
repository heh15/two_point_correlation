### THIS SCRIPT LOOKS AT THE DIFFERENCE IN LOCATION ###
### (GALACTOCENTRIC RADIUS) FOR HIGH VS LOW MASS CLOUDS ####


home_dir = '/home/machado.35/projects/'

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
import pandas as pd
import scipy.stats
from os import listdir
import os
from astropy.io import fits
import astropy.units as u


sources = listdir(home_dir + 'gmc_spacing/data/cloud_cats_fits/homogen_150pc_low_noise/')

    
masses = []
for i in range(len(sources)):
    cat = Table.read(home_dir + 'gmc_spacing/data/cloud_cats_fits/homogen_150pc_low_noise/'+sources[i])
    masses.extend(cat['MLUM_MSUN'])
mass_div = np.nanmedian(masses)
mass_div_high = np.nanpercentile(masses, 84)

high_mass_dist = []
low_mass_dist = []
all_dist = []
for i in range(len(sources)):
    cat = Table.read(home_dir + 'gmc_spacing/data/cloud_cats_fits/homogen_150pc_low_noise/'+sources[i])
    for k in range(len(cat)):
        if cat['MLUM_MSUN'][k] > mass_div:
            high_mass_dist.append(cat['RGAL_KPC'][k])
            all_dist.append(cat['RGAL_KPC'][k])
        else:
            low_mass_dist.append(cat['RGAL_KPC'][k])
            all_dist.append(cat['RGAL_KPC'][k])
            
high_mass_dist = np.array(high_mass_dist)
low_mass_dist = np.array(low_mass_dist)
all_dist = np.array(all_dist)

high_med = np.median(high_mass_dist)
high_err = np.array([[high_med - np.nanpercentile(high_mass_dist, 16)], [np.nanpercentile(high_mass_dist, 84) - high_med]])

low_med = np.median(low_mass_dist)
low_err = np.array([[low_med - np.nanpercentile(low_mass_dist, 16)], [np.nanpercentile(low_mass_dist, 84) - low_med]])

### MAKE PLOTS ###

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.size': 12})

bins = np.linspace(0,15,41)
plt.figure(figsize=(5,5), dpi=150)
plt.hist(high_mass_dist, bins=bins, label='High Mass', alpha=1, density=True, zorder=1, histtype='stepfilled', lw=2,
         color='lightsteelblue', edgecolor='steelblue')
plt.hist(low_mass_dist, bins=bins, label='Low Mass', alpha=0.8, density=True, zorder=2, histtype='step', lw=3, color='darkred')
plt.ylabel('Probability Density')
plt.xlabel('Galactocentric Radius (kpc)')
plt.scatter(high_med, 0.02, zorder=4, color='steelblue')
plt.errorbar(high_med, 0.02, ls='', xerr=high_err, capsize=3, zorder=3, color='steelblue')
plt.scatter(low_med, 0.03, zorder=4, color='darkred')
plt.errorbar(low_med, 0.03, ls='', xerr=low_err, capsize=3, zorder=3, color='darkred')
plt.tick_params(direction='in', which='both', zorder=10)
plt.legend()
plt.grid(ls=':', zorder=0, alpha=0.6)
plt.savefig(home_dir + 'gmc_spacing/figures/mass_gal_radius.pdf', bbox_inches='tight')


### LOOKING AT THIS DIFFERENCE FOR MOST MASSIVE GALAXIES IN SAMPLE ###
tab1 = Table.read(home_dir + 'gmc_spacing/data/phangs_sample_table_v1p6.fits')
gal_masses = tab1['props_mstar']
gal_mass_cutoff = np.nanpercentile(gal_masses, 50)

high_mass_dist = []
low_mass_dist = []
all_dist = []
for i in range(len(sources)):
    cat = Table.read(home_dir + 'gmc_spacing/data/cloud_cats_fits/homogen_150pc_low_noise/'+sources[i])
    name = sources[i].split('_')[0]
    loc = np.where(tab1['name'] == name)[0][0]
    if tab1['props_mstar'][loc] > gal_mass_cutoff:
        for k in range(len(cat)):
            if cat['MLUM_MSUN'][k] > mass_div:
                high_mass_dist.append(cat['RGAL_KPC'][k])
                all_dist.append(cat['RGAL_KPC'][k])
            else:
                low_mass_dist.append(cat['RGAL_KPC'][k])
                all_dist.append(cat['RGAL_KPC'][k])
                
high_mass_dist = np.array(high_mass_dist)
low_mass_dist = np.array(low_mass_dist)
all_dist = np.array(all_dist)

high_med = np.median(high_mass_dist)
high_err = np.array([[high_med - np.nanpercentile(high_mass_dist, 16)], [np.nanpercentile(high_mass_dist, 84) - high_med]])

low_med = np.median(low_mass_dist)
low_err = np.array([[low_med - np.nanpercentile(low_mass_dist, 16)], [np.nanpercentile(low_mass_dist, 84) - low_med]])


bins = np.linspace(0,15,41)
plt.figure(figsize=(5,5), dpi=150)
plt.hist(high_mass_dist, bins=bins, label='High Mass', alpha=1, density=True, zorder=1, histtype='stepfilled', lw=2,
         color='lightsteelblue', edgecolor='steelblue')
plt.hist(low_mass_dist, bins=bins, label='Low Mass', alpha=0.8, density=True, zorder=2, histtype='step', lw=3, color='darkred')
plt.ylabel('Probability Density')
plt.xlabel('Galactocentric Radius (kpc)')
plt.scatter(high_med, 0.02, zorder=4, color='steelblue')
plt.errorbar(high_med, 0.02, ls='', xerr=high_err, capsize=3, zorder=3, color='steelblue')
plt.scatter(low_med, 0.03, zorder=4, color='darkred')
plt.errorbar(low_med, 0.03, ls='', xerr=low_err, capsize=3, zorder=3, color='darkred')
plt.tick_params(direction='in', which='both', zorder=10)
plt.legend()
plt.grid(ls=':', zorder=0, alpha=0.6)
plt.savefig(home_dir + 'gmc_spacing/figures/mass_gal_radius_hmgals.pdf', bbox_inches='tight')