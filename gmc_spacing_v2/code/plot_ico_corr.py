### This script measures correlations with I_CO
### and produces associated plots


home_dir = '/home/machado.35/projects/'


from astropy.table import Table
import pandas as pd
from os import listdir
import re
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib.lines import Line2D


#Location of Circumcloud Tables
circum_fp = home_dir + 'gmc_spacing/data/circumcloud/'

sources = listdir(home_dir + 'gmc_spacing/data/cloud_cats_csv/')
for i in range(len(sources)):
    sources[i] = sources[i].split('_')[0]

circum_tables = listdir(circum_fp)

masses = np.zeros(10000)
ind1 = 0
for i in range(len(sources)):
    cat_fp = home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
    cat = pd.read_csv(cat_fp)
    ind2 = ind1 + len(cat)
    mass = cat['MLUM_MSUN']
    masses[ind1:ind2] = mass
    ind1 = ind2

med_mass = np.nanmedian(masses)
mean_mass = np.nanmean(masses)
mass_16th = np.nanpercentile(masses, 16)
mass_84th = np.nanpercentile(masses, 84)


neighbor_counts = []
nn, nn3, nn5 = [],[],[]
gas = []

for i in range(len(sources)):
    name = sources[i]
    noise = '46mK'
    res = '150pc'
    
    circum_table = circum_fp + name.upper() + '_circumcloud_stats_500pc_' + noise + '_' + res +'.fits'
    if os.path.isfile(circum_table):
        source_cat = pd.read_csv(home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv')
        neighbor_counts.extend(source_cat['clouds_r500pc'])
        nn.extend(source_cat['nn_dist'])
        nn3.extend(source_cat['nn3_dist'])
        nn5.extend(source_cat['nn5_dist'])
        circum_table = Table.read(circum_table)
        gas.extend(circum_table['I_CO21'])
        

        
neighbor_counts = np.array(neighbor_counts) +1
        
#Remove NANS
gas = np.array(gas)
neighbor_counts = np.array(neighbor_counts)
nn = np.array(nn)
nn3 = np.array(nn3)
nn5 = np.array(nn5)


finite_ind = np.isfinite(gas)
gas = gas[finite_ind]
neighbor_counts = neighbor_counts[finite_ind]
nn = nn[finite_ind]
nn3 = nn3[finite_ind]
nn5 = nn5[finite_ind]

#nan neighbors
neighbor_counts = neighbor_counts[np.isfinite(nn5)]
gas = gas[np.isfinite(nn5)]
nn = nn[np.isfinite(nn5)]
nn3 = nn3[np.isfinite(nn5)]
nn5 = nn5[np.isfinite(nn5)]


#neighbor_counts = neighbor_counts[gas != 0]
#neighbor_counts = neighbor_counts[gas < 0]
nn = nn[gas != 0]
nn3 = nn3[gas != 0]
nn5 = nn5[gas != 0]
neighbor_counts = neighbor_counts[gas !=0]
gas = gas[gas != 0]


nn = nn[gas > 0]
nn3 = nn3[gas > 0]
nn5 = nn5[gas > 0]
neighbor_counts = neighbor_counts[gas > 0]
gas = gas[gas > 0]

nn = np.log10(nn)
nn3 = np.log10(nn3)
nn5 = np.log10(nn5)
gas = np.log10(gas)


ico_bins = np.linspace(-0.25,2,10)
bin_asign= np.digitize(gas, bins=ico_bins)

nc_meds = np.zeros_like(ico_bins)
nc_lows = np.zeros_like(ico_bins)
nc_highs = np.zeros_like(ico_bins)


for i in range(len(ico_bins)):
    vals = neighbor_counts[np.where(bin_asign == i)]
    nc_meds[i] = np.nanmean(vals)
    nc_lows[i] = np.nanpercentile(vals, 16)
    nc_highs[i] = np.nanpercentile(vals, 84)

nc_meds = np.log10(nc_meds)
nc_lows = np.log10(nc_lows)
nc_highs = np.log10(nc_highs)
    
nc_err = [nc_meds-nc_lows, nc_highs - nc_meds]


meds1 = np.zeros(len(ico_bins))
lows1 = np.zeros(len(ico_bins))
highs1 = np.zeros(len(ico_bins))
for i in range(len(ico_bins)):
    nn_vals = nn[np.where(bin_asign == i)]
    meds1[i] = np.nanmedian(nn_vals)
    lows1[i] = np.nanpercentile(nn_vals, 16)
    highs1[i] = np.nanpercentile(nn_vals, 84)
    
    
meds3 = np.zeros(len(ico_bins))
lows3 = np.zeros(len(ico_bins))
highs3 = np.zeros(len(ico_bins))
for i in range(len(ico_bins)):
    nn_vals = nn3[np.where(bin_asign == i)]
    meds3[i] = np.nanmedian(nn_vals)
    lows3[i] = np.nanpercentile(nn_vals, 16)
    highs3[i] = np.nanpercentile(nn_vals, 84)
    
    
meds5 = np.zeros(len(ico_bins))
lows5 = np.zeros(len(ico_bins))
highs5 = np.zeros(len(ico_bins))
for i in range(len(ico_bins)):
    nn_vals = nn5[np.where(bin_asign == i)]
    meds5[i] = np.nanmedian(nn_vals)
    lows5[i] = np.nanpercentile(nn_vals, 16)
    highs5[i] = np.nanpercentile(nn_vals, 84)
    
yerr1 = [meds1-lows1, highs1-meds1]
yerr3 = [meds3-lows3, highs3-meds3]
yerr5 = [meds5-lows5, highs5-meds5]

fitting_gas = gas[np.where(neighbor_counts != 0)]
fitting_nc = np.log10(neighbor_counts[neighbor_counts != 0])


log_ico = np.linspace(-3, 4, 15)
cloud_mass = np.array([3e5, 8e5, 3e6, 8e6])
ico_mass_conv = 6.7



fig, ax = plt.subplots(2,2, figsize=(12,14), dpi=150)
xs = np.linspace(-15,15,1000)

spr = scipy.stats.spearmanr(gas, neighbor_counts)[0]
pfit = np.polyfit(fitting_gas, np.log10(fitting_nc), 1)


#pfit = np.polyfit(gas, np.log10(neighbor_counts), 1)
ax[0,0].scatter(gas, np.log10(neighbor_counts), alpha=1, s=8, c='mintcream', edgecolor='darkgrey', label=r'$\rho =$ '+str(np.around(spr, 2)), rasterized=True)

#ax[0,0].plot(log_ico, log_ico + np.log10(500**2) - np.log10(mean_mass), label='Mean', lw=3, alpha=0.5, c='black')
ax[0,0].plot(log_ico, log_ico + np.log10(500**2) - np.log10(med_mass), label=r'Median: 3.2$\times$10$^6$ M$_\odot$', lw=2, ls='--', alpha=0.75, c='black')
ax[0,0].plot(log_ico, log_ico + np.log10(500**2) - np.log10(mass_16th), lw=1, alpha=1, c='red')
ax[0,0].plot(log_ico, log_ico + np.log10(500**2) - np.log10(mass_84th), lw=1, alpha=1, c='red')
ax[0,0].fill_between(log_ico, log_ico + np.log10(500**2) - np.log10(mass_16th), log_ico + np.log10(500**2) - np.log10(mass_84th), color='grey', alpha=0.4, zorder=1)

#ax[0,0].plot(log_ico, log_ico + np.log10(500**2) - np.log10(cloud_mass[0]), label=r'3$\times$10$^5$ M$_\odot$', alpha=0.5)
#ax[0,0].plot(log_ico, log_ico + np.log10(500**2) - np.log10(cloud_mass[1]), label=r'8$\times$10$^5$ M$_\odot$', alpha=0.5)
#ax[0,0].plot(log_ico, log_ico + np.log10(500**2) - np.log10(cloud_mass[2]), label=r'3$\times$10$^6$ M$_\odot$', alpha=0.5)
#ax[0,0].plot(log_ico, log_ico + np.log10(500**2) - np.log10(cloud_mass[3]), label=r'8$\times$10$^6$ M$_\odot$', alpha=0.5)


ax[0,0].errorbar(ico_bins, nc_meds, yerr=nc_err, xerr=None, ls='', capsize=4, lw=3, marker='X', color='navy', markerfacecolor='black', zorder=20)

#ax[0,0].errorbar(nc_meds, ncs, yerr=None, xerr=nc_err, ls='', capsize=4, zorder=24, lw=3, marker='^', color='seagreen', markerfacecolor='black')
#ax[0,0].plot(xs, pfit[0]*xs+pfit[1], color='salmon', lw=2)
ax[0,0].set_ylabel('$log$(Neighbor Counts +1) within 500pc')
ax[0,0].set_xlabel(r'$log_{10}(I_{CO})$ within 500pc')
ax[0,0].set_ylim(-0.5,2)
ax[0,0].set_xlim(-2,3.5)
#ax[0,0].set_yscale('log')
ax[0,0].tick_params(direction='in', which='both')
ax[0,0].legend()

spr1 = scipy.stats.spearmanr(gas, nn)[0]
spr3 = scipy.stats.spearmanr(gas, nn3)[0]
spr5 = scipy.stats.spearmanr(gas, nn5)[0]
pfit1 = np.polyfit(gas, nn, 1)
pfit1_med = np.polyfit(ico_bins, meds1, 1)
pfit3 = np.polyfit(gas, nn3, 1)
pfit3_med = np.polyfit(ico_bins, meds3, 1)
pfit5 = np.polyfit(gas, nn5, 1)
pfit5_med = np.polyfit(ico_bins, meds5, 1)

ax[0,1].grid(zorder=0)
ax[0,1].scatter(ico_bins, meds1, c='lightsteelblue', edgecolor='black', s=55, zorder=25, ls='None')
ax[0,1].errorbar(ico_bins, meds1, xerr=0, yerr=yerr1, ls='', c='darkblue', lw=3, capsize=4, zorder=24)
ax[0,1].scatter(gas, nn, alpha=1, zorder=10, s=8, c='mistyrose', edgecolor='rosybrown', linewidth=1, label=r'$\rho =$ '+str(np.around(spr1, 2)), rasterized=True)
#ax[0,1].plot(xs, pfit1[0]*xs+pfit1[1], color='black', lw=3, ls='dashed', zorder=20)
ax[0,1].plot(xs, pfit1_med[0]*xs+pfit1_med[1], color='black', lw=3, ls='dashed', zorder=20)
ax[0,1].set_xlim(-2,3)
ax[0,1].set_ylim(0.6,4.2)
ax[0,1].set_ylabel(r'log$_{10}$(Distance to N neighbor) (pc)')
ax[0,1].set_xlabel(r'$log_{10}(I_{CO})$ within 500pc')
ax[0,1].tick_params(direction='in', which='both')

legend_elements = [Line2D([0], [0], marker='o', markersize=4, markerfacecolor='mistyrose', color='rosybrown', lw=0, 
                          label=r'First: $\rho =$ '+str(np.around(spr1, 2))+r', $\alpha$ = '+str(np.around(pfit1[0], 3)))]


ax[0,1].legend(handles=legend_elements)







ax[1,0].grid(zorder=0)
ax[1,0].scatter(ico_bins, meds3, c='lightsteelblue', edgecolor='black', s=35, zorder=25, ls='None')
ax[1,0].errorbar(ico_bins, meds3, xerr=0, yerr=yerr3, ls='', c='darkblue', lw=3, capsize=4, zorder=24)
ax[1,0].scatter(gas, nn3, alpha=1, zorder=10, s=8, c='mistyrose', edgecolor='rosybrown', label=r'$\rho =$ '+str(np.around(spr3, 2)), rasterized=True)
#ax[1,0].plot(xs, pfit3[0]*xs+pfit3[1], color='black', lw=3, ls='dashed', zorder=20)
ax[1,0].plot(xs, pfit3_med[0]*xs+pfit3_med[1], color='black', lw=3, ls='dashed', zorder=20)
ax[1,0].set_ylabel(r'log$_{10}$(Distance to N neighbor) (pc)')
ax[1,0].set_xlabel(r'$log_{10}(I_{CO})$ within 500pc')
ax[1,0].set_xlim(-2,3)
ax[1,0].set_ylim(0.6,4.2)
ax[1,0].tick_params(direction='in', which='both')

legend_elements = [Line2D([0], [0], marker='o', markersize=4, markerfacecolor='mistyrose', color='rosybrown', lw=0, 
                          label=r'Third: $\rho =$ '+str(np.around(spr3, 2))+r', $\alpha$ = '+str(np.around(pfit3[0], 3)))]
ax[1,0].legend(handles=legend_elements)
            

ax[1,1].grid(zorder=0)
ax[1,1].scatter(ico_bins, meds5, c='lightsteelblue', edgecolor='black', s=35, zorder=25, ls='None')
ax[1,1].errorbar(ico_bins, meds5, xerr=0, yerr=yerr5, ls='', c='darkblue', lw=3, capsize=4, zorder=24)
ax[1,1].scatter(gas, nn5, alpha=1, zorder=10, s=8, c='mistyrose', edgecolor='rosybrown', label=r'$\rho =$ '+str(np.around(spr5, 2)), rasterized=True)
#ax[1,1].plot(xs, pfit5[0]*xs+pfit5[1], color='black', lw=3, ls='dashed', zorder=20)
ax[1,1].plot(xs, pfit5_med[0]*xs+pfit5_med[1], color='black', lw=3, ls='dashed', zorder=20)
ax[1,1].set_ylabel(r'log$_{10}$(Distance to N neighbor) (pc)')
ax[1,1].set_xlabel(r'$log_{10}(I_{CO})$ within 500pc')
ax[1,1].set_xlim(-2,3)
ax[1,1].set_ylim(0.6,4.2)
ax[1,1].tick_params(direction='in', which='both')

legend_elements = [Line2D([0], [0], marker='o', markersize=4, markerfacecolor='mistyrose', color='rosybrown', lw=0, 
                         label=r'Fifth: $\rho =$ '+str(np.around(spr5, 2))+r', $\alpha$ = '+str(np.around(pfit5[0], 3)))]

ax[1,1].legend(handles=legend_elements)

plt.savefig(home_dir + 'gmc_spacing/figures/multi_panel_ico_corr.pdf', bbox_inches='tight')
