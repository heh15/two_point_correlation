#####
## This script generates an individual clustering curve for each galaxy
## and saves it in a .csv file. This script also produces clustering curve plots
## with uncertainties for each galaxy.


## THIS TAKES A WHILE TO RUN ##
## Of order several hours - maybe even longer than a day##
## This finds uncertainties on clustering amplitude by measuring clustering 
## for exponential vs exponential cases 1000x.

home_dir = '/home/machado.35/projects/'

import numpy as np
import pandas as pd
from astropy.table import Table
import matplotlib.pyplot as plt
from os import listdir
import os
import sys
sys.path.append(home_dir + 'gmc_spacing/modules')
from clust_uncertainty import clustering_uncertainty
from deprojectGalaxy import deproject
import scipy.stats


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


sources = listdir(home_dir + 'gmc_spacing/data/cloud_cats_csv')
for i in range(len(sources)):
    sources[i] = sources[i].split('_')[0]

galTable = Table.read(home_dir + 'gmc_spacing/data/phangs_sample_table_v1p6.fits')

### CALCULATING INDIVIDUAL GALAXY CURVES 
bins = np.linspace(0,6500,44)
clust_amp, clust_amp_cov = np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins)))
DD, RR, DR = np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins)))
RR_cov, DR_cov = np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins)))
n_real, n_rand, n_real_rand = np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins)))
n_cov, n_real_cov = np.zeros((len(sources), len(bins))), np.zeros((len(sources), len(bins)))
for i in range(len(sources)):
    cat = home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
    rand_cat = home_dir + 'gmc_spacing/data/sim_results/exp/'+sources[i]+'_exp_150pc.csv'
    cov_cat = home_dir + 'gmc_spacing/data/sim_results/flat/'+sources[i]+'_flat_150pc.csv'

    name = sources[i].split('_')[0]
    row = np.where(galTable['name'] == sources[i].split('_')[0])[0][0]
    galRA=galTable[row]['orient_ra']
    galDEC=galTable[row]['orient_dec']
    pa=galTable[row]['orient_posang']
    incl=galTable[row]['orient_incl']
    galDist=galTable[row]['dist'] #Mpc
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

                    
names = []
for i in sources:
    name = i.split('_')[0]
    names.append(name)

df = pd.DataFrame()
df['galname'] = np.array(names)
for i in range(len(bins)):
    df['bin'+str(i)] = clust_amp[:,i]
df.to_csv(home_dir + 'gmc_spacing/data/gal_clust_curves.csv', index = False)
                    

    

### CALCULATING UNCERTAINTIES FOR CURVES
t = Table.read(home_dir + 'gmc_spacing/data/phangs_sample_table_v1p6.fits')
bins = np.linspace(0,6500,44)
colnames = []
for i in range(len(bins)):
    colnames.append('bin '+ str(i))
DD, RR, DR = np.zeros_like(bins), np.zeros_like(bins), np.zeros_like(bins)
clust_amp = np.zeros_like(bins)
galaxies = {}
n_loops = 1000

for i in range(len(sources)):
    
    ##Open Catalogs##
    
    # Catalog of REAL clouds #
    real = home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
    # Simulated Catalog #
    cat = home_dir + 'gmc_spacing/data/sim_results/exp/'+sources[i]+'_exp_150pc.csv'
    
    
    # Galaxy Coords from Sample Table -- Needed for Deprojection #
    name = sources[i].split('_')[0]
    row = np.where(t['name'] == sources[i].split('_')[0])[0][0]
    galRA=t[row]['orient_ra']
    galDEC=t[row]['orient_dec']
    pa=t[row]['orient_posang']
    incl=t[row]['orient_incl']
    galDist=t[row]['dist'] #Mpc
    center_coord=[galRA, galDEC]
    
    
    if os.path.isfile(cat):
        # Check to make sure files exist #
        real = pd.read_csv(real)
        cat = pd.read_csv(cat)
        if len(real) > 5:
            df = pd.DataFrame(columns=colnames)
            # Make sure the real cloud catalog has more than 5 objects #
            for loop_num in range(n_loops):
                
                # Set number of control clouds equal to number of real clouds #
                num_real_clouds = len(real)
                start = loop_num * num_real_clouds
                end = start + num_real_clouds
                
                # Subset of control to measure clustering against - n_subset = n_realclouds #
                data_x = cat['rRA'][start:end]
                data_y = cat['rDEC'][start:end]
                
                # Control clouds -- capped a 5000 for compute times #
                rand_x = cat['rRA'][0:5000]
                rand_y = cat['rDEC'][0:5000]
                
                # Deproject control subset #
                deproj_data = deproject(center_coord=center_coord, incl=incl, pa=pa, ra=data_x, dec=data_y, return_offset=True)
                
                # Deproject full control sample #
                deproj_rand = deproject(center_coord=center_coord, incl=incl, pa=pa, ra=rand_x, dec=rand_y, return_offset=True)

                # Number of subset objects and number of control objects #
                ND = len(deproj_data[2])
                NR = len(deproj_rand[2])
                
                # Initialize cross-comparison arrays -- used to store DR calculations #
                data_rand = np.zeros((len(deproj_data[2]),len(deproj_rand[2])))
                
                for k in range(len(deproj_data[2])):
                    # Calculate distance between every subset object and every control object #
                    data_rand[k] = np.sqrt((deproj_rand[2]-deproj_data[2][k])**2 + (deproj_rand[3]-deproj_data[3][k])**2)
                
                
                # Put distances in physical units #
                data_rand = np.deg2rad(data_rand) * galDist * 1e6
                
                # Calculate real-real, subset-subset, and control-control object distances with physical units #
                data_sources = np.deg2rad(pair_dist_calc(deproj_data[2], deproj_data[3])) * galDist * 1e6
                rand_sources = np.deg2rad(pair_dist_calc(deproj_rand[2], deproj_rand[3])) * galDist * 1e6

                
                # Calculate clustering amplitude within each bin #
                for ii in range(len(bins)-1):
                    
                    # Setting bounds of bin #
                    sep_val_low = bins[ii]
                    sep_val_high =bins[ii+1]
                    ## REMOVE NAN VALS ##
                    data_sources = data_sources[np.isfinite(data_sources)]
                    rand_sources = rand_sources[np.isfinite(rand_sources)]

                    ## COUNTS WITHIN BIN ##
                    n_data = len(np.where(np.logical_and(data_sources>=sep_val_low, data_sources<=sep_val_high))[0])
                    n_rand = len(np.where(np.logical_and(rand_sources>=sep_val_low, rand_sources<=sep_val_high))[0])
                    n_data_rand = len(np.where(np.logical_and(data_rand>=sep_val_low, data_rand<=sep_val_high))[0])

                    # Calculate clustering amplitude terms for subset #
                    DD[ii] = n_data / (ND * (ND-1) / 2)
                    RR[ii] = n_rand / (NR * (NR-1)/2)
                    DR[ii] = n_data_rand / (ND*NR)

                    if RR[ii] != 0:
                        clust_amp[ii] = (1/RR[ii]) * ((DD[ii]) - (2*DR[ii]) + RR[ii]) + 1
                    else:
                        clust_amp[ii] = np.nan
                        
                df.loc[loop_num] = clust_amp
            df.to_csv(home_dir + 'gmc_spacing/data/clust_curves/'+name+'.csv', index=False)
            print(name+ ' done.')

        else:
            galaxies[name] = 'Less than 5 clouds'
            print(name)
            
            
            
clust_curves = home_dir + 'gmc_spacing/data/clust_curves/'
clust_curves = listdir(clust_curves)
real_curves = home_dir + 'gmc_spacing/data/gal_clust_curves.csv'
real_curves = pd.read_csv(real_curves)

bins = np.linspace(0,6500,44)
ax_list = [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2), (3,0), (3,1), (3,2)]
ii = 0
fignum = 1
fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(12,9), dpi=125)
fig.tight_layout(pad=1.0)
for i in range(len(real_curves.name)):
    real_curve = real_curves.iloc[i][1:]
    rand_curves = pd.read_csv(home_dir + 'gmc_spacing/data/clust_curves/'+clust_curves[i])
    upper = []
    lower = []
    for j in rand_curves.columns:
        upper.append(np.nanpercentile(rand_curves[j], 84))
        lower.append(np.nanpercentile(rand_curves[j], 16))
    upper = np.array(upper)
    lower = np.array(lower)
    bins_trim = bins[0:43]
    real_trim = real_curve[0:43]
    upper_trim = (upper+real_curve-1)[0:43].astype(float)
    lower_trim = (lower+real_curve-1)[0:43].astype(float)
    if ii >= len(ax_list):
        plt.savefig(home_dir + 'gmc_spacing/figures/gal_clust_curves_'+str(fignum)+'.pdf', bbox_inches='tight')
        fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(12,9), dpi=125)
        fig.tight_layout(pad=1.0)
        ii = 0
        ax[ax_list[ii]].plot(bins_trim,real_trim, c='darkred', zorder=3, label=real_curves.name[i].upper())
        ax[ax_list[ii]].fill_between(bins_trim, upper_trim, lower_trim, alpha=0.6, color='slategrey')
        if ax_list[ii][0]==3:
            ax[ax_list[ii]].set_xlabel('Spacing (pc)')
        if ax_list[ii][1]==0:
            ax[ax_list[ii]].set_ylabel(r'Clustering Amplitude (1+$\omega$)')
        ax[ax_list[ii]].hlines(1,1,6500, ls='--', color='black', zorder=4, lw=2)
        ax[ax_list[ii]].set_xlim(0,6300)
        ax[ax_list[ii]].tick_params(direction='in', which='both')
        ax[ax_list[ii]].legend()
        ii += 1
        fignum+=1
    else:
        ax[ax_list[ii]].plot(bins_trim,real_trim, c='darkred', zorder=3, label=real_curves.name[i].upper())
        ax[ax_list[ii]].fill_between(bins_trim, upper_trim, lower_trim, alpha=0.6, color='slategrey')
        if ax_list[ii][0]==3:
            ax[ax_list[ii]].set_xlabel('Spacing (pc)')
        if ax_list[ii][1]==0:
            ax[ax_list[ii]].set_ylabel(r'Clustering Amplitude (1+$\omega$)')
        ax[ax_list[ii]].hlines(1,1,6500, ls='--', color='black', zorder=4, lw=2)
        ax[ax_list[ii]].set_xlim(0,6300)
        ax[ax_list[ii]].tick_params(direction='in', which='both')
        ax[ax_list[ii]].legend()
        ii += 1
        
plt.savefig(gmc_spacing + 'gmc_spacing/figures/gal_clust_curves_'+str(fignum+1)+'.pdf', bbox_inches='tight')
    