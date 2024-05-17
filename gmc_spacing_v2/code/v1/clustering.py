##################################################
## This script calculates clustering amplitudes ##
########### for full sample of galaxies ###########

home_dir = '/home/machado.35/projects/'


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from astropy.table import Table
from os import listdir
import sys
import scipy.stats
sys.path.append('/home/machado.35/projects/modules')
from deprojectGalaxy import deproject


def pair_dist_calc(x,y):
    n = len(x)
    x1 = np.tile(x, (n,1)) #constant x along column
    y1 = np.tile(y, (n,1)) #constant y along column
    x2 = np.transpose(x1) #constant x along rows
    y2 = np.transpose(y1) #constant y along rows
    dist = np.sqrt((x1-x2)**2 + (y1-y2)**2)
    dist = np.tril(dist) #DISCARDS BELOW DIAGONAL TO AVOID DOUBLE COUNTING
    dist[np.where(dist == 0)] = np.nan
    return(dist)


### SETUP ###
sources = listdir(home_dir+'gmc_spacing/data/cloud_cats_csv/')
for i in range(len(sources)):
    sources[i] = sources[i].split('_')[0]

t = Table.read(home_dir + 'gmc_spacing/data/phangs_sample_table_v1p6.fits')
    

### CLUSTERING CALCULATIONS ###

bins = np.linspace(0,6500,44)
clust_amp, clust_amp_cov = np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins)))
DD, RR, DR = np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins)))
RR_cov, DR_cov = np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins)))
n_real, n_rand, n_real_rand = np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins)))
n_cov, n_real_cov = np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins)))
plt.figure(figsize=(8,5), dpi=125)
for i in range(len(sources)):
    cat = home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
    rand_cat = home_dir + 'gmc_spacing/data/sim_results/exp/'+sources[i]+'_exp_150pc.csv' # Exponential
    cov_cat = home_dir + 'gmc_spacing/data/sim_results/flat/'+sources[i]+'_flat_150pc.csv' # Flat

    name = sources[i].split('_')[0]
    row = np.where(t['name'] == sources[i].split('_')[0])[0][0]
    galRA=t[row]['orient_ra']
    galDEC=t[row]['orient_dec']
    pa=t[row]['orient_posang']
    incl=t[row]['orient_incl']
    galDist=t[row]['dist'] #Mpc
    center_coord=[galRA, galDEC]
    
    
    if os.path.isfile(cat) and os.path.isfile(rand_cat):
        cat = pd.read_csv(cat)
        if len(cat) > 5:
            rand_cat = pd.read_csv(rand_cat)
            cov_cat = pd.read_csv(cov_cat)
            
            real_x = cat['RA']
            real_y = cat['DEC']
            
            rand_x = rand_cat['rRA'][0:5000]
            rand_y = rand_cat['rDEC'][0:5000]
            
            cov_x = cov_cat['rRA'][0:5000]
            cov_y = cov_cat['rDEC'][0:5000]
            
            deproj_real = deproject(center_coord=center_coord, incl=incl, pa=pa, ra=real_x, dec=real_y, return_offset=True)
            deproj_rand = deproject(center_coord=center_coord, incl=incl, pa=pa, ra=rand_x, dec=rand_y, return_offset=True)
            deproj_cov = deproject(center_coord=center_coord, incl=incl, pa=pa, ra=cov_x, dec=cov_y, return_offset=True)

            real_rand = np.zeros((len(deproj_real[2]),len(deproj_rand[2])))
            real_cov = np.zeros((len(deproj_real[2]),len(deproj_cov[2])))
            for k in range(len(deproj_real[2])):
                real_rand[k] = np.sqrt((deproj_rand[2]-deproj_real[2][k])**2 + (deproj_rand[3]-deproj_rand[3][k])**2)
                real_cov[k] = np.sqrt((deproj_cov[2]-deproj_real[2][k])**2 + (deproj_cov[3]-deproj_cov[3][k])**2)
            real_rand = np.deg2rad(real_rand) * galDist * 1e6
            real_cov = np.deg2rad(real_cov) * galDist * 1e6
            real_sources = pair_dist_calc(deproj_real[2], deproj_real[3])
            real_sources = np.deg2rad(real_sources) * galDist * 1e6
            rand_sources = pair_dist_calc(deproj_rand[2], deproj_rand[3])
            rand_sources = np.deg2rad(rand_sources) * galDist * 1e6
            cov_sources = pair_dist_calc(deproj_cov[2], deproj_cov[3])
            cov_sources = np.deg2rad(cov_sources) * galDist * 1e6
            for ii in range(len(bins)-1):
                sep_val_low = bins[ii]
                sep_val_high =bins[ii+1]
                ## REMOVE NAN VALS ##
                real_sources = real_sources[np.isfinite(real_sources)]
                rand_sources = rand_sources[np.isfinite(rand_sources)]
                cov_sources = cov_sources[np.isfinite(cov_sources)]

                ## COUNTS WITHIN BIN ##
                n_real[i][ii] = len(np.where(np.logical_and(real_sources>=sep_val_low, real_sources<=sep_val_high))[0])
                n_rand[i][ii] = len(np.where(np.logical_and(rand_sources>=sep_val_low, rand_sources<=sep_val_high))[0])
                n_real_rand[i][ii] = len(np.where(np.logical_and(real_rand>=sep_val_low, real_rand<=sep_val_high))[0])
                
                n_cov[i][ii] = len(np.where(np.logical_and(cov_sources>=sep_val_low, cov_sources<=sep_val_high))[0])
                n_real_cov[i][ii] = len(np.where(np.logical_and(real_cov>=sep_val_low, real_cov<=sep_val_high))[0])
                
                ND = len(deproj_real[2])
                NR = len(deproj_rand[2])
                NR_cov = len(deproj_cov[2])

                DD[i][ii] = n_real[i][ii] / (ND * (ND-1) / 2)
                RR[i][ii] = n_rand[i][ii] / (NR * (NR-1)/2)
                DR[i][ii] = n_real_rand[i][ii] / (ND*NR)
                
                RR_cov[i][ii] = n_cov[i][ii] / (NR_cov * (NR_cov-1)/2)
                DR_cov[i][ii] = n_real_cov[i][ii] / (ND*NR_cov)
                if RR[i][ii] != 0:
                    clust_amp[i][ii] = (1/RR[i][ii]) * ((DD[i][ii]) - (2*DR[i][ii]) + RR[i][ii] ) + 1
                else:
                    clust_amp[i][ii] = np.nan
                
                if RR_cov[i][ii] != 0:
                    clust_amp_cov[i][ii] = (1/RR_cov[i][ii]) * ((DD[i][ii]) - (2*DR_cov[i][ii]) + RR_cov[i][ii] ) + 1
                else:
                    clust_amp_cov[i][ii] = np.nan
                    
                    

### SET PARAMETERS FOR PLOTTING UNCERTAINTIES ###


med = np.zeros(len(bins))
for i in range(len(bins)):
    med[i] = np.nanmedian(clust_amp[:,i][clust_amp[:,i]!=0])
    
low = np.zeros(len(bins))
for i in range(len(bins)):
    if len(clust_amp[:,i][clust_amp[:,i]!=0]) < 2:
        low[i] = 0
    else:
        low[i] = np.nanpercentile(clust_amp[:,i][clust_amp[:,i]!=0], 16)

        
high = np.zeros(len(bins))
for i in range(len(bins)):
    if len(clust_amp[:,i][clust_amp[:,i]!=0]) < 2:
        high[i] = 0
    else:
        high[i] = np.nanpercentile(clust_amp[:,i][clust_amp[:,i]!=0], 84)
        
        
        
med_cov = np.zeros(len(bins))
for i in range(len(bins)):
    med_cov[i] = np.nanmedian(clust_amp_cov[:,i][clust_amp_cov[:,i]!=0])
    
low_cov = np.zeros(len(bins))
for i in range(len(bins)):
    if len(clust_amp_cov[:,i][clust_amp_cov[:,i]!=0]) < 2:
        low_cov[i] = 0
    else:
        low_cov[i] = np.nanpercentile(clust_amp_cov[:,i][clust_amp_cov[:,i]!=0], 16)

        
high_cov = np.zeros(len(bins))
for i in range(len(bins)):
    if len(clust_amp_cov[:,i][clust_amp_cov[:,i]!=0]) < 2:
        high_cov[i] = 0
    else:
        high_cov[i] = np.nanpercentile(clust_amp_cov[:,i][clust_amp_cov[:,i]!=0], 84)

yerr = [abs(low-med), (high-med)]
yerr_cov = [abs(low_cov-med_cov), (high_cov-med_cov)]




### MAKE CLUSTERING PLOT ###

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.size': 12})

fig, ax = plt.subplots(1,2, figsize=(10,4), dpi=150)
ax[0].axvspan(158, 505, alpha=.3, hatch='//', color='gray', zorder=1)
ax[0].axvspan(158, 505, alpha=.3, hatch='//', color='None', edgecolor='gray', zorder=1)
ax[0].plot(bins, med_cov, label='Flat', c='navy', zorder=5, lw=3)
ax[0].plot(bins, med, label='Exponential', c='red', zorder=5, lw=3)
ax[0].fill_between(bins, low, high, color='red', edgecolor='None', alpha=0.2, zorder=3)

ax[0].vlines(305, ymin=-1, ymax=5, ls=':', color='black', alpha=0.7)

#ax[0].errorbar(bins, med, yerr=yerr, ls='none', capsize=3, c='red', zorder=4)
ax[0].hlines(y=1, xmin=-1000, xmax=7000, color='black', ls='--', lw=2, alpha=0.6, zorder=2)
ax[0].set_xlim(0,1300)
ax[0].tick_params(direction='in', which='both')
ax[0].set_xlabel('Cloud Spacing (pc)')
ax[0].set_ylabel(r'1 + $\omega (\Delta r)$')
ax[0].grid(axis='both', zorder=0, alpha=1, ls=':')
ax[0].set_ylim(0.25,4)
ax[0].legend()

ax[1].vlines(305, ymin=-3, ymax=15, ls=':', color='black', alpha=0.7)
ax[1].plot(bins, med_cov, label='Flat', c='navy', zorder=5, lw=3)
ax[1].plot(bins, med, label='Exponential', c='red', zorder=5, lw=3)
ax[1].fill_between(bins, low_cov, high_cov, color='navy', edgecolor='None', alpha=0.1, zorder=3)
ax[1].hlines(y=1, xmin=-1000, xmax=7000, color='black', ls='--', lw=2, alpha=0.6, zorder=2)
ax[1].set_xlim(0,6300)
ax[1].set_ylim(0,6)
ax[1].tick_params(direction='in', which='both')
ax[1].set_xlabel('Cloud Spacing (pc)')
ax[1].set_ylabel(r'1 + $\omega (\Delta r)$')
ax[1].grid(axis='both', zorder=0, alpha=1, ls=':')
ax[1].legend()

plt.savefig(home_dir + 'gmc_spacing/figures/clustering_amp_2panel.pdf', bbox_inches='tight')


fig = plt.figure(figsize=(5,4), dpi=150)
plt.plot(bins, med_cov, label='Flat', c='navy', zorder=5, lw=3)
plt.plot(bins, med, label='Exponential', c='red', zorder=5, lw=3)
plt.fill_between(bins, low, high, color='red', edgecolor='None', alpha=0.2, zorder=3)
plt.vlines(305, ymin=-1, ymax=5, ls=':', color='black', alpha=0.7)
#ax[0].errorbar(bins, med, yerr=yerr, ls='none', capsize=3, c='red', zorder=4)
plt.hlines(y=1, xmin=-1000, xmax=7000, color='black', ls='--', lw=2, alpha=0.6, zorder=2)
plt.xlim(0,6300)
plt.tick_params(direction='in', which='both')
plt.xlabel('Cloud Spacing (pc)')
plt.ylabel(r'1 + $\omega (\Delta r)$')
plt.grid(axis='both', zorder=0, alpha=1, ls=':')
plt.ylim(0.6,1.4)
plt.legend()
plt.savefig(home_dir + 'gmc_spacing/figures/clustering_stability.pdf', bbox_inches='tight')