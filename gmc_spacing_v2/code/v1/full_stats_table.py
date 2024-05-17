### This script makes the full stats table
### Includes galaxy name, N_gmcs, Nth nearest neighbor distances,
### peak clustering amplitude, and spatial scale of peak amplitude
### currently Table 5 in draft (as of 6/2/22)



home_dir = '/home/machado.35/projects/'

from astropy.table import Table
import astropy.units as u
import numpy as np
from os import listdir
import os
import pandas as pd
import math
from astropy.io import fits


def search_curve(curve, threshold):
    peak = np.argmax(curve)
    sub_curve = curve - threshold
    upper = []
    lower = []
    i = 1
    while i < len(curve):
        if sub_curve[peak+i] > 0:
            upper.append(peak+i)
            i+=1
        else:
            if sub_curve[peak+i] < 0:
                i=len(curve)
    i = 1
    while i < len(curve):
        if sub_curve[peak-i] > 0:
            lower.append(peak-i)
            i+=1
        else:
            if sub_curve[peak-i]:
                i = len(curve)
    
    return(upper, lower)



## GALAXY NAMES ##
galname = listdir(home_dir + 'gmc_spacing/data/cloud_cats_csv/')
for i in range(len(galname)):
    galname[i] = galname[i].split('_')[0]
galname = np.array(galname, dtype=str)


## GALAXY DISTANCES ##
t = Table.read(home_dir + 'gmc_spacing/data/phangs_sample_table_v1p6.fits')
galdist = np.zeros(len(galname))
for i in range(len(galname)):
    for k in range(len(t)):
        if t['name'][k] == galname[i]:
            galdist[i] = t['dist'][k]
galdist = np.around(galdist, 2)


## GALAXY MASSES ##
galmass = np.zeros(len(galname))
for i in range(len(galname)):
    for k in range(len(t)):
        if t['name'][k] == galname[i]:
            galmass[i] = t['props_mstar'][k]
galmass = np.around(np.log10(galmass), 2)

## NUM CLOUDS ##
numclouds = np.zeros(len(galname), dtype=int)
for i in range(len(galname)):
    cat = home_dir + 'gmc_spacing/data/cloud_cats_csv/' +str(galname[i])+'_v4cprops_distcalc_46mK_150pc.csv'
    if os.path.isfile(cat):
        cat = pd.read_csv(cat)
        numclouds[i] = len(cat['nn_dist'])

## NEIGHBOR DISTANCES ##
first = np.zeros(len(galname))
first_upper = np.zeros(len(galname))
first_lower = np.zeros(len(galname))
third = np.zeros(len(galname))
third_upper = np.zeros(len(galname))
third_lower = np.zeros(len(galname))
fifth = np.zeros(len(galname))
fifth_upper = np.zeros(len(galname))
fifth_lower = np.zeros(len(galname))

for i in range(len(galname)):
    cat = '/data/rubin/machado.35/phangs/products/max_cat/homogen_150pc_low_noise/'+str(galname[i])+'_v4cprops_distcalc_46mK_150pc.csv'
    if os.path.isfile(cat):
        cat = pd.read_csv(cat)
        first[i] = np.nanmedian(cat['nn_dist'])
        first_upper[i] = np.nanpercentile(cat['nn_dist'], 84)
        first_lower[i] = np.nanpercentile(cat['nn_dist'], 16)
        
        third[i] = np.nanmedian(cat['nn3_dist'])
        third_upper[i] = np.nanpercentile(cat['nn3_dist'], 84)
        third_lower[i] = np.nanpercentile(cat['nn3_dist'], 16)
        
        fifth[i] = np.nanmedian(cat['nn5_dist'])
        fifth_upper[i] = np.nanpercentile(cat['nn5_dist'], 84)
        fifth_lower[i] = np.nanpercentile(cat['nn5_dist'], 16)

first = np.around(first.astype(int), -1)
first_upper = np.around(first_upper.astype(int), -1)
first_lower = np.around(first_lower.astype(int), -1)

third = np.around(third.astype(int), -1)
third_upper = np.around(third_upper.astype(int), -1)
third_lower = np.around(third_lower.astype(int), -1)

fifth = np.around(fifth.astype(int), -1)
fifth_upper = np.around(fifth_upper.astype(int), -1)
fifth_lower = np.around(fifth_lower.astype(int), -1)



first_col = []
for i in range(len(first)):
    first_col.append(str(first[i])+r'$^{+'+str(first_upper[i]-first[i])+r'}_{-'+str(first[i]-first_lower[i])+r'}$')
    
third_col = []
for i in range(len(third)):
    third_col.append(str(third[i])+r'$^{+'+str(third_upper[i]-third[i])+r'}_{-'+str(third[i]-third_lower[i])+r'}$')
    
fifth_col = []
for i in range(len(fifth)):
    fifth_col.append(str(fifth[i])+r'$^{+'+str(fifth_upper[i]-fifth[i])+r'}_{-'+str(fifth[i]-fifth_lower[i])+r'}$')


### PEAK CLUSTERING AMPLITUDE ###
clust_peak = np.zeros(len(galname))

### SPACING AT CLUSTERING PEAK ###
peak_scale = np.zeros(len(galname), dtype=int)

### PEAK CLUSTERING RANGE > 1.125 ###
peak_scale_upper = np.zeros(len(galname), dtype=int)
peak_scale_lower = np.zeros(len(galname), dtype=int)

clust_curves = pd.read_csv(home_dir + 'gmc_spacing/data/gal_clust_curves.csv')
for i in range(len(galname)):
    curve = clust_curves.loc[i][1:]
    clust_peak[i] = np.around(np.nanmax(curve), 2)
    peak_scale[i] = int((np.argmax(curve)+1) * 150)
    upper, lower = search_curve(curve, 1.125)
    if len(upper) == 0:
        upper = np.argmax(curve)
    if len(lower) == 0:
        lower = np.argmax(curve)
    peak_scale_upper[i] = 150*np.max(upper)+1
    peak_scale_lower[i] = 150*np.min(lower)+1

peak_scale_col = []
for i in range(len(peak_scale)):
    peak_scale_col.append(str(peak_scale[i])+r'$^{+'+str(peak_scale_upper[i]-peak_scale[i])+r'}_{-'+str(peak_scale[i]-peak_scale_lower[i])+r'}$')
    
for i in range(len(galname)):
    galname[i] = galname[i].upper()

galtab = Table([galname, numclouds, first_col, third_col, fifth_col, clust_peak, peak_scale],
     names = ('Galaxy Name', 'Number of GMCs', '1st Nearest', '3rd Nearest', '5th Nearest', 'Clustering Amplitude Peak',
             'Clustering Peak Scale'),
      units= ([None,None,u.pc,u.pc,u.pc, r'1+$\omega$', u.pc])
     )



galtab.write(home_dir + 'gmc_spacing/tables/galaxy_clustering_table.tex', format='latex')