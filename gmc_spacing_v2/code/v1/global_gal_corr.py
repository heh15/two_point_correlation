## THIS MEASURES AND PLOTS CORRELATIONS WITH
## GLOBAL GALAXY PROPERTIES
## IT ALSO MAKES A CORRESPONDING TABLE

home_dir = '/home/machado.35/projects/'

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import matplotlib.colors as mcolors
from astropy.table import Table
from os import listdir
import sys
import scipy.stats
sys.path.append(home_dir + 'gmc_spacing/modules')
from deprojectGalaxy import deproject
from clust_uncertainty import clustering_uncertainty
import matplotlib as mpl



galTable = Table.read(home_dir + 'gmc_spacing/data/phangs_sample_table_v1p6.fits')
galmass = galTable['props_mstar']
galsfr = galTable['props_sfr']
gallco = galTable['lco_phangs']


sources = listdir(home_dir + 'gmc_spacing/data/cloud_cats_csv/')
names = []
for i in sources:
    name = i.split('_')[0]
    names.append(name)

table_indicies = []
for i in names:
    table_indicies.append(np.where(galTable['name'] == i)[0][0])
table_indicies = np.array(table_indicies)


first, third, fifth = [], [], []
for i in sources:
    cat = pd.read_csv(home_dir + 'gmc_spacing/data/cloud_cats_csv/'+i)
    first.append(np.nanmedian(cat['nn_dist']))
    third.append(np.nanmedian(cat['nn3_dist']))
    fifth.append(np.nanmedian(cat['nn5_dist']))
    
first = np.array(first)
third = np.array(third)
fifth = np.array(fifth)
    
    
types = []
for i in galTable[table_indicies]['morph_string']:
    if i not in types:
        types.append(i)

early_type, barred, no_bar = [], [], []
for i in galTable[table_indicies]['morph_string']:
    if '0' in i:
        if i not in early_type:
            early_type.append(i)
    if 'B' in i:
        if i not in barred:
            barred.append(i)
    else:
        if '0' not in i:
            if i not in no_bar:
                no_bar.append(i)

## EARLY TYPE ##
ind = [] 
for i in early_type:
    ind.extend(np.where(galTable[table_indicies]['morph_string'] == i))
ind = np.concatenate(ind)
early_type = ind

## BARRED ##
ind = [] 
for i in barred:
    ind.extend(np.where(galTable[table_indicies]['morph_string'] == i))
ind = np.concatenate(ind)
barred = ind

## NO-BAR ##
ind = [] 
for i in no_bar:
    ind.extend(np.where(galTable[table_indicies]['morph_string'] == i))
ind = np.concatenate(ind)
no_bar = ind

clust_curves = pd.read_csv(home_dir + 'gmc_spacing/data/gal_clust_curves.csv')


clust_peaks = np.zeros(len(clust_curves))
for i in range(len(clust_curves)):
    clust_peaks[i] = np.nanmax(clust_curves.loc()[i][1:])
    
clust_upper_err = np.zeros(len(clust_curves))
clust_lower_err = np.zeros_like(clust_upper_err)

for i in range(len(clust_curves)):
    peak_loc = np.argmax(clust_curves.loc()[i][1:])
    uncs = clustering_uncertainty(clust_curves.loc()[i][0])
    clust_upper_err[i] = uncs[0][peak_loc]
    clust_lower_err[i] = uncs[1][peak_loc]
    

yerr = np.array([clust_upper_err-1, 1-clust_lower_err])

## MAKE PLOT ##
plt.rcParams['axes.linewidth'] = 2
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.size': 20})

fig, ax = plt.subplots(2,3, figsize=(20,12), dpi=125)
plt.rcParams['axes.linewidth'] = 2
c = galTable[table_indicies]['morph_t']
s= 125
cmap = 'inferno'

#Galaxy Mass
corr = scipy.stats.spearmanr(first, galmass[table_indicies])[0]
ax[0,0].scatter(galmass[table_indicies], first, c='None', alpha=0.7, label=r'$\rho =$ '+str(np.around(corr, 2)))
ax[0,0].scatter(galmass[table_indicies][early_type], first[early_type], marker='s', c=c[early_type], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s)
ax[0,0].scatter(galmass[table_indicies][barred], first[barred], marker='x', c=c[barred], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, lw=3.5)
ax[0,0].scatter(galmass[table_indicies][no_bar], first[no_bar], marker='D', c=c[no_bar], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s)
ax[0,0].set_ylabel('Characteristic GMC Spacing')
ax[0,0].set_xscale('log')
ax[0,0].set_yscale('log')
ax[0,0].grid(ls=':')
ax[0,0].legend()


corr = scipy.stats.spearmanr(clust_peaks, galmass[table_indicies])[0]
ax[1,0].scatter(galmass[table_indicies], clust_peaks, c='None', alpha=0.7, label=r'$\rho =$ '+str(np.around(corr, 2)))
ax[1,0].scatter(galmass[table_indicies][early_type], clust_peaks[early_type], marker='s', c=c[early_type], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, zorder=5)
ax[1,0].errorbar(galmass[table_indicies], clust_peaks, yerr=yerr, capsize=5, c='None', ecolor='black', zorder=0)
ax[1,0].scatter(galmass[table_indicies][barred], clust_peaks[barred], marker='x', c=c[barred], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, lw=3.5, zorder=5)
ax[1,0].scatter(galmass[table_indicies][no_bar], clust_peaks[no_bar], marker='D', c=c[no_bar], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, zorder=5)
ax[1,0].set_ylabel('Clustering Amplitude Peak')
ax[1,0].set_xlabel(r'Stellar Mass (M$_{\odot}$)')
ax[1,0].set_xscale('log')
ax[1,0].set_yscale('log')
ax[1,0].grid(ls=':')
ax[1,0].legend()

#LCO
corr = scipy.stats.spearmanr(first[np.isfinite(gallco[table_indicies])], gallco[table_indicies][np.isfinite(gallco[table_indicies])])[0]
ax[0,1].scatter(gallco[table_indicies], first, c='None', alpha=0.7, label=r'$\rho =$ '+str(np.around(corr, 2)))
ax[0,1].scatter(gallco[table_indicies][early_type], first[early_type], marker='s', c=c[early_type], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s)
ax[0,1].scatter(gallco[table_indicies][barred], first[barred], marker='x', c=c[barred], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, lw=3.5)
ax[0,1].scatter(gallco[table_indicies][no_bar], first[no_bar], marker='D', c=c[no_bar], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s)
ax[0,1].set_xscale('log')
ax[0,1].set_yscale('log')
ax[0,1].set_xlim(2000000, 2000000000)
ax[0,1].grid(ls=':')
ax[0,1].legend()


corr = scipy.stats.spearmanr(clust_peaks[np.isfinite(gallco[table_indicies])], gallco[table_indicies][np.isfinite(gallco[table_indicies])])[0]
ax[1,1].scatter(gallco[table_indicies], clust_peaks, c='None', alpha=0.7, label=r'$\rho =$ '+str(np.around(corr, 2)))
ax[1,1].errorbar(gallco[table_indicies], clust_peaks, yerr=yerr, capsize=5, c='None', ecolor='black', zorder=0)
ax[1,1].scatter(gallco[table_indicies][early_type], clust_peaks[early_type], marker='s', c=c[early_type], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, zorder=5)
ax[1,1].scatter(gallco[table_indicies][barred], clust_peaks[barred], marker='x', c=c[barred], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, lw=3.5, zorder=5)
ax[1,1].scatter(gallco[table_indicies][no_bar], clust_peaks[no_bar], marker='D', c=c[no_bar], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, zorder=5)
ax[1,1].set_xscale('log')
ax[1,1].set_yscale('log')
ax[1,1].set_xlabel(r'L$_{CO}$ (K km pc$^2$ s$^{-1}$)')
ax[1,1].set_xlim(2000000, 2000000000)
ax[1,1].grid(ls=':')
ax[1,1].legend()

#SFR
corr = scipy.stats.spearmanr(first, galsfr[table_indicies])[0]
ax[0,2].scatter(galsfr[table_indicies], first, c='None', alpha=0.7, label=r'$\rho =$ '+str(np.around(corr, 2)))
ax[0,2].scatter(galsfr[table_indicies][early_type], first[early_type], marker='s', c=c[early_type], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, zorder=5)
ax[0,2].scatter(galsfr[table_indicies][barred], first[barred], marker='x', c=c[barred], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, lw=3.5, zorder=5)
ax[0,2].scatter(galsfr[table_indicies][no_bar], first[no_bar], marker='D', c=c[no_bar], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, zorder=5)
ax[0,2].set_xscale('log')
ax[0,2].set_yscale('log')
ax[0,2].grid(ls=':')
ax[0,2].legend()


corr = scipy.stats.spearmanr(clust_peaks, galsfr[table_indicies])[0]
ax[1,2].scatter(galsfr[table_indicies], clust_peaks, c='None', alpha=0.7, label=r'$\rho =$ '+str(np.around(corr, 2)))
ax[1,2].errorbar(galsfr[table_indicies], clust_peaks, yerr=yerr, capsize=5, c='None', ecolor='black', zorder=0)
sc = ax[1,2].scatter(galsfr[table_indicies][early_type], clust_peaks[early_type], marker='s', c=c[early_type], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, zorder=5, label='Early Type')
ax[1,2].scatter(galsfr[table_indicies][barred], clust_peaks[barred], marker='x', c=c[barred], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, lw=3.5, zorder=5, label='Late Type: Barred')
ax[1,2].scatter(galsfr[table_indicies][no_bar], clust_peaks[no_bar], marker='D', c=c[no_bar], cmap=cmap, vmin=0, vmax=8, alpha=0.7, s=s, zorder=5, label='Late Type: No Bar')
ax[1,2].set_xlabel(r'SFR (M$_{\odot}$ /yr)')
ax[1,2].set_xscale('log')
ax[1,2].set_yscale('log')
ax[1,2].grid(ls=':')
ax[1,2].legend(prop={'size': 14})

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
norm = mpl.colors.Normalize(vmin=0, vmax=8)
fig.colorbar(sc, cax=cbar_ax, cmap=cmap, norm=norm, label=r'Morphological T Type')
plt.savefig(home_dir + 'gmc_spacing/figures/global_gal_corr_plot.pdf', bbox_inches='tight')



### FOR TABLE ###
items = ['First Nearest Neighbor',
         'Third Nearest Neighbor',
         'Fifth Nearest Neighbor',
         'Clustering Amplitude']

late_type = np.concatenate([no_bar, barred])

mgal_all_first = scipy.stats.spearmanr(first, galmass[table_indicies])[0]
mgal_all_third = scipy.stats.spearmanr(third, galmass[table_indicies])[0]
mgal_all_fifth = scipy.stats.spearmanr(fifth, galmass[table_indicies])[0]
mgal_all_clust = scipy.stats.spearmanr(clust_peaks, galmass[table_indicies])[0]
mgal_all = [mgal_all_first, mgal_all_third, mgal_all_fifth, mgal_all_clust]
mgal_all = np.around(mgal_all, 2)

mgal_early_first = scipy.stats.spearmanr(first[early_type], galmass[table_indicies][early_type])[0]
mgal_early_third = scipy.stats.spearmanr(third[early_type], galmass[table_indicies][early_type])[0]
mgal_early_fifth = scipy.stats.spearmanr(fifth[early_type], galmass[table_indicies][early_type])[0]
mgal_early_clust = scipy.stats.spearmanr(clust_peaks[early_type], galmass[table_indicies][early_type])[0]
mgal_early = [mgal_early_first, mgal_early_third, mgal_early_fifth, mgal_early_clust]
mgal_early = np.around(mgal_early, 2)

mgal_late_first = scipy.stats.spearmanr(first[late_type], galmass[table_indicies][late_type])[0]
mgal_late_third = scipy.stats.spearmanr(third[late_type], galmass[table_indicies][late_type])[0]
mgal_late_fifth = scipy.stats.spearmanr(fifth[late_type], galmass[table_indicies][late_type])[0]
mgal_late_clust = scipy.stats.spearmanr(clust_peaks[late_type], galmass[table_indicies][late_type])[0]
mgal_late = [mgal_late_first, mgal_late_third, mgal_late_fifth, mgal_late_clust]
mgal_late = np.around(mgal_late, 2)

lco_all_first = scipy.stats.spearmanr(first[~np.isnan(gallco[table_indicies])], gallco[table_indicies][~np.isnan(gallco[table_indicies])])[0]
lco_all_third = scipy.stats.spearmanr(third[~np.isnan(gallco[table_indicies])], gallco[table_indicies][~np.isnan(gallco[table_indicies])])[0]
lco_all_fifth = scipy.stats.spearmanr(fifth[~np.isnan(gallco[table_indicies])], gallco[table_indicies][~np.isnan(gallco[table_indicies])])[0]
lco_all_clust = scipy.stats.spearmanr(clust_peaks[~np.isnan(gallco[table_indicies])], gallco[table_indicies][~np.isnan(gallco[table_indicies])])[0]
lco_all = [lco_all_first, lco_all_third, lco_all_fifth, lco_all_clust]
lco_all = np.around(lco_all, 2)

lco_early_first = scipy.stats.spearmanr(first[early_type][~np.isnan(gallco[table_indicies][early_type])],
                                       gallco[table_indicies][early_type][~np.isnan(gallco[table_indicies][early_type])])[0]
lco_early_third = scipy.stats.spearmanr(third[early_type][~np.isnan(gallco[table_indicies][early_type])],
                                       gallco[table_indicies][early_type][~np.isnan(gallco[table_indicies][early_type])])[0]
lco_early_fifth = scipy.stats.spearmanr(fifth[early_type][~np.isnan(gallco[table_indicies][early_type])],
                                       gallco[table_indicies][early_type][~np.isnan(gallco[table_indicies][early_type])])[0]
lco_early_clust = scipy.stats.spearmanr(clust_peaks[early_type][~np.isnan(gallco[table_indicies][early_type])],
                                       gallco[table_indicies][early_type][~np.isnan(gallco[table_indicies][early_type])])[0]
lco_early = [lco_early_first, lco_early_third, lco_early_fifth, lco_early_clust]
lco_early = np.around(lco_early, 2)

lco_late_first = scipy.stats.spearmanr(first[late_type][~np.isnan(gallco[table_indicies][late_type])],
                                       gallco[table_indicies][late_type][~np.isnan(gallco[table_indicies][late_type])])[0]
lco_late_third = scipy.stats.spearmanr(third[late_type][~np.isnan(gallco[table_indicies][late_type])],
                                       gallco[table_indicies][late_type][~np.isnan(gallco[table_indicies][late_type])])[0]
lco_late_fifth = scipy.stats.spearmanr(fifth[late_type][~np.isnan(gallco[table_indicies][late_type])],
                                       gallco[table_indicies][late_type][~np.isnan(gallco[table_indicies][late_type])])[0]
lco_late_clust = scipy.stats.spearmanr(clust_peaks[late_type][~np.isnan(gallco[table_indicies][late_type])],
                                       gallco[table_indicies][late_type][~np.isnan(gallco[table_indicies][late_type])])[0]
lco_late = [lco_late_first, lco_late_third, lco_late_fifth, lco_late_clust]
lco_late = np.around(lco_late, 2)

sfr_all_first = scipy.stats.spearmanr(first, galsfr[table_indicies])[0]
sfr_all_third = scipy.stats.spearmanr(third, galsfr[table_indicies])[0]
sfr_all_fifth = scipy.stats.spearmanr(fifth, galsfr[table_indicies])[0]
sfr_all_clust = scipy.stats.spearmanr(clust_peaks, galsfr[table_indicies])[0]
sfr_all = [sfr_all_first, sfr_all_third, sfr_all_fifth, sfr_all_clust]
sfr_all = np.around(sfr_all, 2)

sfr_early_first = scipy.stats.spearmanr(first[early_type], galsfr[table_indicies][early_type])[0]
sfr_early_third = scipy.stats.spearmanr(third[early_type], galsfr[table_indicies][early_type])[0]
sfr_early_fifth = scipy.stats.spearmanr(fifth[early_type], galsfr[table_indicies][early_type])[0]
sfr_early_clust = scipy.stats.spearmanr(clust_peaks[early_type], galsfr[table_indicies][early_type])[0]
sfr_early = [sfr_early_first, sfr_early_third, sfr_early_fifth, sfr_early_clust]
sfr_early = np.around(sfr_early, 2)

sfr_late_first = scipy.stats.spearmanr(first[late_type], galsfr[table_indicies][late_type])[0]
sfr_late_third = scipy.stats.spearmanr(third[late_type], galsfr[table_indicies][late_type])[0]
sfr_late_fifth = scipy.stats.spearmanr(fifth[late_type], galsfr[table_indicies][late_type])[0]
sfr_late_clust = scipy.stats.spearmanr(clust_peaks[late_type], galsfr[table_indicies][late_type])[0]
sfr_late = [sfr_late_first, sfr_late_third, sfr_late_fifth, sfr_late_clust]
sfr_late = np.around(sfr_late, 2)


columns = [items, mgal_all, mgal_early, mgal_late, 
           lco_all, lco_early, lco_late, 
           sfr_all, sfr_early, sfr_late]
colnames = ('','All', 'Early Type', 'Late Type', 
            'All ', 'Early Type ', 'Late Type ', 
            ' All', ' Early Type', ' Late Type')
units = [None, None, None, None, None, None, None, None, None]

gal_corr_tab = Table(columns, names=colnames)
gal_corr_tab.write(home_dir + 'gmc_spacing/tables/gal_corr_tab.tex', format='latex', overwrite=True)


tar -zcvf archive-name.tar.gz source-directory-name