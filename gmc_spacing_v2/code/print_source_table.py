### SOURCE TABLE ###
# Galaxy Name
# Distance
# Galaxy Mass
# Completeness
# N_GMCs

home_dir = '/home/machado.35/projects/'

from astropy.table import Table
import astropy.units as u
import numpy as np
from os import listdir
import os
import pandas as pd
import math
from astropy.io import fits

## GALAXY NAMES ##
galname = sorted(listdir(home_dir + 'gmc_spacing/data/cloud_cats_csv/'))
for i in range(len(galname)):
    galname[i] = galname[i].split('_')[0].upper()
galname = np.array(galname, dtype=str)


## GALAXY DISTANCES ##
t = Table.read(home_dir + 'gmc_spacing/data/phangs_sample_table_v1p6.fits')
galdist = np.zeros(len(galname))
for i in range(len(galname)):
    for k in range(len(t)):
        if t['name'][k] == galname[i].lower():
            galdist[i] = t['dist'][k]
galdist = np.around(galdist, 1)


## GALAXY MASSES ##
galmass = np.zeros(len(galname))
for i in range(len(galname)):
    for k in range(len(t)):
        if t['name'][k] == galname[i].lower():
            galmass[i] = t['props_mstar'][k]
galmass = np.around(np.log10(galmass), 1)


## COMPLETENESS ##
sources = listdir(home_dir + 'gmc_spacing/data/cloud_cats_csv/')
for i in range(len(sources)):
    sources[i] = sources[i].split('_')[0]

comp_vals = np.zeros(len(sources))
for i in range(len(sources)):
    
    row = np.where(t['name'] == sources[i])[0][0]
    t_val = t['lco_phangs'][row]
    dist = t['dist'][row]
    
    
    mask_fp = home_dir + 'gmc_spacing/data/mom0/'+sources[i]+'_12m+7m+tp_co21_150pc_broad_mom0.fits'
    cat_fp = home_dir + 'gmc_spacing/data/cloud_cats_fits/homogen_150pc_low_noise/'+sources[i]+'_12m+7m+tp_co21_150pc_46mK_props.fits'
    if os.path.isfile(mask_fp) and os.path.isfile(cat_fp):
        mask = fits.open(mask_fp)
        cat = Table.read(cat_fp)
    else:
        mask_fp = home_dir + 'gmc_spacing/data/mom0/'+sources[i]+'_7m+tp_co21_150pc_broad_mom0.fits'
        cat_fp = home_dir + 'gmc_spacing/data/cloud_cats_fits/homogen_150pc_low_noise/'+sources[i]+'_7m+tp_co21_150pc_46mK_props.fits'
        mask = fits.open(mask_fp)
        cat = Table.read(cat_fp)

    ### MASK ###
    total_map = np.sum(mask[0].data[np.isfinite(mask[0].data)])
    #total_map = total_map *u.K *u.km / u.s

    pix_s1 = (abs(np.deg2rad(mask[0].header['CDELT1']))*dist*u.Mpc).to(u.pc)
    pix_s2 = (abs(np.deg2rad(mask[0].header['CDELT2']))*dist*u.Mpc).to(u.pc)
    pix_size = pix_s1 * pix_s2
    total_map = total_map * (pix_size.value)
    
    ### CAT ###
    cat_max = np.sum(cat['FLUX_NOEX']) #* u.K * u.km / u.s *u.pc**2
    
    
    
    completeness = np.around(cat_max/total_map, 2)
    comp_vals[i] = completeness

    
## NUMCLOUDS ##
numclouds = np.zeros(len(sources))
for i in range(len(sources)):
    cat = pd.read_csv(home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv')
    numclouds[i] = int(len(cat))

    
## BUILD TABLE ##
columns = [galname, galdist, galmass, comp_vals, numclouds.astype(int)]
colnames = ('Galaxy Name', 'Distance (Mpc)', r'Mass ($log_{10}$$M_{\odot}$)', 'Completeness', 'Number of Clouds')
units = ([None, u.Mpc, u.Msun, None, None])
tab = Table(columns, names=colnames)

tab.write(home_dir + 'gmc_spacing/tables/tab1.tex', format='latex', overwrite=True)