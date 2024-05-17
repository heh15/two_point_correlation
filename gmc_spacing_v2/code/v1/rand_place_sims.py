#########################################################################################################################
#########################################################################################################################
### THIS SCRIPT RUNS MONTE CARLO SIMULATIONS TO GENERATE CLOUD POPULATIONS BASED ON A GIVEN HYPOTHETICAL DISTRIBUTION ###
######################## THIS SCRIPT ALSO CALCULATES ALL SPACING STATISTIS FOR SIMULATED CATALOGS #######################
#########################################################################################################################
#########################################################################################################################


home_dir = '/home/machado.35/projects/' #Set this as the location of where the gmc_spacing folder exists


#### IMPORTS ####
import numpy as np
from astropy.table import Table
import pandas as pd
import os.path
import sys
from os import listdir
sys.path.append(home_dir+'gmc_spacing/modules')
from deprojectGalaxy import deproject
import spacing_calc


### SETUP ###
sources = listdir(home_dir+'gmc_spacing/data/cloud_cats_csv/')
for i in range(len(sources)):
    sources[i] = sources[i].split('_')[0]
res = [150]

tab_loc = home_dir+'gmc_spacing/data/phangs_sample_table_v1p6.fits'
t = Table.read(tab_loc)

## Determine number of realizations
#####################################
## Chosing a fixed number of realizations that is standardized across ALL galaxies.
#####################################
## Number of realizations is set such that each galaxy will have a total of 
## AT LEAST 1e6 random placements across all realizations.
#####################################

res = [150]
nclouds = []
for k in range(len(res)):
    for i in range(len(sources)):
        image = home_dir+'gmc_spacing/data/masks/coverage/'+sources[i]+'_12m+7m+tp_co21_'+str(res[k])+'pc_coverage2d.fits'
        cat_fp = home_dir+'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
        if os.path.isfile(image):
            if os.path.isfile(cat_fp):
                real_cat = pd.read_csv(cat_fp)
                nclouds.append(len(real_cat))

### This is set by catalog with LARGEST number of clouds
n_iters = int(1e6/max(nclouds))



#### RUN RANDOM PLACEMENTS ####

#Flat
output_path = home_dir+'gmc_spacing/data/sim_results/flat/'

for k in range(len(res)):
    for i in range(len(sources)):
        image = home_dir+'gmc_spacing/data/masks/coverage/'+sources[i]+'_12m+7m+tp_co21_'+str(res[k])+'pc_coverage2d.fits'
        cat_fp = home_dir+'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
        if os.path.isfile(image):
            pass
        else:
            image = home_dir + 'gmc_spacing/data/masks/coverage/'+sources[i]+'_7m+tp_co21_'+str(res[k])+'pc_coverage2d.fits' # ONLY FOR NGC0253
            if not os.path.isfile(image):
                image = home_dir + 'gmc_spacing/data/masks/coverage/'+sources[i]+'_7m_co21_'+str(res[k])+'pc_coverage2d.fits' # ONLY FOR NGC0300
            if os.path.isfile(cat_fp):
                row = int(np.where(t['name'] == sources[i])[0])
                galRA=t[row]['orient_ra']
                galDEC=t[row]['orient_dec']
                pa=t[row]['orient_posang']
                incl=t[row]['orient_incl']
                galDist=t[row]['dist']
                outfile_name = output_path+sources[i]+'_flat'+'_'+str(res[k])+'pc.csv'
                real_cat = pd.read_csv(cat_fp)
                n_real = len(real_cat)
                n_clouds = n_iters * n_real
                spacing_calc.place_nclouds_in_map(image, n_clouds, outfile_name, galRA, galDEC, pa, incl, galDist)

#Exp
output_path = home_dir+'gmc_spacing/data/sim_results/exp/'

n_clouds = 1e6
for k in range(len(res)):
    for i in range(len(sources)):
        image = home_dir+'gmc_spacing/data/masks/exp/'+sources[i]+'_'+str(res[k])+'pc_exp_mask.fits'
        cat_fp = home_dir+'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
        if os.path.isfile(image):
            row = int(np.where(t['name'] == sources[i])[0])
            galRA=t[row]['orient_ra']
            galDEC=t[row]['orient_dec']
            pa=t[row]['orient_posang']
            incl=t[row]['orient_incl']
            galDist=t[row]['dist']
            real_cat = pd.read_csv(cat_fp)
            n_real = len(real_cat)
            n_clouds = n_iters * n_real
            outfile_name = output_path+sources[i]+'_exp'+'_'+str(res[k])+'pc.csv'
            spacing_calc.place_nclouds_in_map(image, n_clouds, outfile_name, galRA, galDEC, pa, incl, galDist)
            
            

### CALCULATE SPACING STATISTICS FOR SIMULATED CATALOGS ###

#FLAT
radii = [300, 500, 700]
for k in range(len(res)):
    for i in range(len(sources)):
        image = home_dir+'gmc_spacing/data/masks/coverage/'+sources[i]+'_12m+7m+tp_co21_'+str(res[k])+'pc_coverage2d.fits'
        cat_fp = home_dir+'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
        sim_fp = home_dir+'gmc_spacing/data/sim_results/flat/'+sources[i]+'_flat_'+str(res[k])+'pc.csv'
        if os.path.isfile(image):
            pass
        else:
            image = home_dir + 'gmc_spacing/data/masks/coverage/'+sources[i]+'_7m+tp_co21_'+str(res[k])+'pc_coverage2d.fits' # ONLY FOR NGC0253
            if not os.path.isfile(image):
                image = home_dir + 'gmc_spacing/data/masks/coverage/'+sources[i]+'_7m_co21_'+str(res[k])+'pc_coverage2d.fits' # ONLY FOR NGC0300
            if os.path.isfile(cat_fp):
                if os.path.isfile(sim_fp):
                    row = int(np.where(t['name'] == sources[i])[0])
                    dist = t['dist'][row]*1e6 #from Mpc to pc
                    df = pd.read_csv(sim_fp) #Simulated catalog
                    cat = pd.read_csv(cat_fp) #Real catalog
                    x = df['rRA']
                    y = df['rDEC']
                    N_tot = len(x)
                    n_it = len(cat)
                    n_loops = int(np.floor(N_tot / n_it))
                    mindist = np.zeros_like(x)*np.nan
                    mindist3 = np.zeros_like(x)*np.nan
                    mindist5 = np.zeros_like(x)*np.nan
                    radial_density_300 = np.zeros_like(x)*np.nan
                    radial_density_500 = np.zeros_like(x)*np.nan
                    radial_density_700 = np.zeros_like(x)*np.nan
                    for ii in np.arange(n_loops):
                        tx = x[ii*n_it:(ii+1)*n_it]
                        ty = y[ii*n_it:(ii+1)*n_it]
                        mindist[ii*n_it:(ii+1)*n_it] = spacing_calc.pair_dist_calc(x=tx,y=ty, neighbor_num=1)[0] #1st nearest
                        mindist3[ii*n_it:(ii+1)*n_it] = spacing_calc.pair_dist_calc(x=tx,y=ty, neighbor_num=3)[0] #3rd nearest
                        mindist5[ii*n_it:(ii+1)*n_it] = spacing_calc.pair_dist_calc(x=tx,y=ty, neighbor_num=5)[0] #5th nearest
                        radial_density_300[ii*n_it:(ii+1)*n_it] = spacing_calc.rad_density(x=tx,y=ty, rad=radii[0], galdist=dist)
                        radial_density_500[ii*n_it:(ii+1)*n_it] = spacing_calc.rad_density(x=tx,y=ty, rad=radii[1], galdist=dist)
                        radial_density_700[ii*n_it:(ii+1)*n_it] = spacing_calc.rad_density(x=tx,y=ty, rad=radii[2], galdist=dist)
                    df['nn_dist'] = np.deg2rad(mindist) * dist
                    df['nn3_dist'] = np.deg2rad(mindist3) * dist
                    df['nn5_dist'] = np.deg2rad(mindist5) * dist
                    df['clouds_r'+str(radii[0])+'pc'] = radial_density_300
                    df['clouds_r'+str(radii[1])+'pc'] = radial_density_500
                    df['clouds_r'+str(radii[2])+'pc'] = radial_density_700
                    df.to_csv(sim_fp)

                    
                    
#EXP DISK
radii = [300,500,700]
for k in range(len(res)):
    for i in range(len(sources)):
        image = home_dir+'gmc_spacing/data/masks/exp/'+sources[i]+'_'+str(res[k])+'pc_exp_mask.fits'
        cat_fp = home_dir+'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
        sim_fp = home_dir+'gmc_spacing/data/sim_results/exp/'+sources[i]+'_exp_'+str(res[k])+'pc.csv'
        if os.path.isfile(image):
            if os.path.isfile(cat_fp):
                if os.path.isfile(sim_fp):
                    row = int(np.where(t['name'] == sources[i])[0])
                    dist = t['dist'][row]*1e6 #from Mpc to pc
                    df = pd.read_csv(sim_fp) #Simulated catalog
                    cat = pd.read_csv(cat_fp) #Real catalog
                    x = df['rRA']
                    y = df['rDEC']
                    N_tot = len(x)
                    n_it = len(cat)
                    n_loops = int(np.floor(N_tot / n_it))
                    mindist = np.zeros_like(x)*np.nan
                    mindist3 = np.zeros_like(x)*np.nan
                    mindist5 = np.zeros_like(x)*np.nan
                    radial_density_300 = np.zeros_like(x)*np.nan
                    radial_density_500 = np.zeros_like(x)*np.nan
                    radial_density_700 = np.zeros_like(x)*np.nan
                    for ii in np.arange(n_loops):
                        tx = x[ii*n_it:(ii+1)*n_it]
                        ty = y[ii*n_it:(ii+1)*n_it]
                        mindist[ii*n_it:(ii+1)*n_it] = spacing_calc.pair_dist_calc(x=tx,y=ty, neighbor_num=1)[0] #1st nearest
                        mindist3[ii*n_it:(ii+1)*n_it] = spacing_calc.pair_dist_calc(x=tx,y=ty, neighbor_num=3)[0] #3rd nearest
                        mindist5[ii*n_it:(ii+1)*n_it] = spacing_calc.pair_dist_calc(x=tx,y=ty, neighbor_num=5)[0] #5th nearest
                        radial_density_300[ii*n_it:(ii+1)*n_it] = spacing_calc.rad_density(x=tx,y=ty, rad=radii[0], galdist=dist)
                        radial_density_500[ii*n_it:(ii+1)*n_it] = spacing_calc.rad_density(x=tx,y=ty, rad=radii[1], galdist=dist)
                        radial_density_700[ii*n_it:(ii+1)*n_it] = spacing_calc.rad_density(x=tx,y=ty, rad=radii[2], galdist=dist)
                    df['nn_dist'] = np.deg2rad(mindist) * dist
                    df['nn3_dist'] = np.deg2rad(mindist3) * dist
                    df['nn5_dist'] = np.deg2rad(mindist5) * dist
                    df['clouds_r'+str(radii[0])+'pc'] = radial_density_300
                    df['clouds_r'+str(radii[1])+'pc'] = radial_density_500
                    df['clouds_r'+str(radii[2])+'pc'] = radial_density_700
                    df.to_csv(sim_fp)
