####### THIS SCRIPT GENERATES HISTOGRAMS OF #######
### NUMBERS OF NEIGHBORS WITHIN 300,500,700 PC ####
####### for all catalogs (real & simulated) #######



home_dir = '/home/machado.35/projects/'

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from astropy.table import Table
from os import listdir

## Setup ##

sources = listdir(home_dir + 'gmc_spacing/data/cloud_cats_csv/')
for i in range(len(sources)):
    sources[i] = sources[i].split('_')[0]


    
radii = [300,500,700]
fig, ax = plt.subplots(1,3, figsize=(20,5), dpi=150)
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.size': 16})
for ii in range(len(radii)):
    nc_real = []
    nc_exp = []
    nc_cov = []
    for i in range(len(sources)):
        name = sources[i].split('_')[0]
        real_data = home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
        rand_2d = home_dir + 'gmc_spacing/data/sim_results/flat/'+sources[i]+'_flat_150pc.csv'
        rand_exp = home_dir + 'gmc_spacing/data/sim_results/exp/'+sources[i]+'_exp_150pc.csv'
        if os.path.isfile(real_data) and os.path.isfile(rand_exp):
            real_data = pd.read_csv(real_data)
            if len(real_data) > 5:
                rand_2d = pd.read_csv(rand_2d)
                rand_exp = pd.read_csv(rand_exp)

                col_label = 'clouds_r'+str(radii[ii])+'pc'

                real_dist = real_data['nn_dist']
                real_nc = real_data[col_label] #neighbor counts
                rand_2d_dist = rand_2d['nn_dist']
                rand_2d_nc = rand_2d[col_label]
                rand_exp_dist = rand_exp['nn_dist']
                rand_exp_nc = np.array(rand_exp[col_label])
                nc_real.extend(real_nc)
                nc_exp.extend(rand_exp_nc)
                nc_cov.extend(rand_2d_nc)
    ncreal = np.array(nc_real)
    ncexp = np.array(nc_exp)
    nccov = np.array(nc_cov)
    bins = np.linspace(0,25,26)
    ax[ii].hist(ncreal, bins=bins, density=True, lw=2, histtype='step', facecolor='grey', hatch='/', edgecolor='black', alpha=0.5, fill=True, label='Real GMCs', zorder=10)
    ax[ii].hist(nccov, bins=bins, density=True, lw=2, histtype='stepfilled', color='blue', alpha=0.3, edgecolor='None', label='Flat')
    ax[ii].hist(ncexp, bins=bins, density=True, lw=2.2, label='Exponential', histtype='stepfilled', color='red', edgecolor='None', alpha=0.3)
    ax[ii].legend()
    ax[ii].set_ylabel('Fraction of GMCs')
    ax[ii].set_xlabel('Number of Neighbors Within '+str(radii[ii])+'pc')
    ax[ii].tick_params(direction='in')
    print('MEAN NUM NEIGHBORS WITHIN '+str(radii[ii])+'pc:')
    print('REAL: '+str(np.nanmean(ncreal)))
    print('2D COVERAGE: '+str(np.nanmean(nccov)))
    print('EXPONENTIAL: '+str(np.nanmean(ncexp)))
    print()
    print('MEDIAN NUM NEIGHBORS WITHIN '+str(radii[ii])+'pc:')
    print('REAL: '+str(np.nanmedian(ncreal)))
    print('2D COVERAGE: '+str(np.nanmedian(nccov)))
    print('EXPONENTIAL: '+str(np.nanmedian(ncexp)))
    ax[ii].grid(alpha=1, ls=':', axis='both', zorder=0)
    print()

ax[0].set_xlim(0,8)
ax[1].set_xlim(0,12)
ax[2].set_xlim(0,18)
plt.savefig(home_dir + 'gmc_spacing/figures/neighbor_density_lo_noise.pdf', bbox_inches='tight')