### THIS SCRIPT GENERATES NEAREST NEIGHBOR ###
#### DISTANCE HISTOGRAMS FOR ALL CATALOGS ####
### IT ALSO GENERATES A SUMMARY TABLE OF ###
### SPACING STATS FOR ALL CATALOGS ###




home_dir = '/home/machado.35/projects/'


from os import listdir
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy import units as u
from astropy.table import Table


sources = listdir(home_dir + 'gmc_spacing/data/cloud_cats_csv/')
for i in range(len(sources)):
    sources[i] = sources[i].split('_')[0]
    

real_first = []
cov_first = []
exp_first = []

real_third = []
cov_third = []
exp_third = []

real_fifth = []
cov_fifth = []
exp_fifth = []

for i in range(len(sources)):
    name = sources[i]
    res = '150pc'
    real_data = home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
    rand_2d = home_dir + 'gmc_spacing/data/sim_results/flat/'+sources[i]+'_flat_150pc.csv'
    rand_exp = home_dir + 'gmc_spacing/data/sim_results/exp/'+sources[i]+'_exp_150pc.csv'
    if os.path.isfile(real_data) and os.path.isfile(rand_exp):
        real_data = pd.read_csv(real_data)
        if len(real_data) > 5:
            rand_2d = pd.read_csv(rand_2d)
            rand_exp = pd.read_csv(rand_exp)


            real_dist = real_data['nn_dist']
            real_dist3 = real_data['nn3_dist']
            real_dist5 = real_data['nn5_dist']

            rand_2d_dist = rand_2d['nn_dist']
            rand_2d_dist3 = rand_2d['nn3_dist']
            rand_2d_dist5 = rand_2d['nn5_dist']

            rand_exp_dist = rand_exp['nn_dist']
            rand_exp_dist3 = rand_exp['nn3_dist']
            rand_exp_dist5 = rand_exp['nn5_dist']

            real_first.extend(real_dist)
            real_third.extend(real_dist3)
            real_fifth.extend(real_dist5)

            exp_first.extend(rand_exp_dist)
            exp_third.extend(rand_exp_dist3)
            exp_fifth.extend(rand_exp_dist5)

            cov_first.extend(rand_2d_dist)
            cov_third.extend(rand_2d_dist3)
            cov_fifth.extend(rand_2d_dist5)

        
real_first = np.array(real_first)
real_third = np.array(real_third)
real_fifth = np.array(real_fifth)

exp_first = np.array(exp_first)
exp_third = np.array(exp_third)
exp_fifth = np.array(exp_fifth)

cov_first = np.array(cov_first)
cov_third = np.array(cov_third)
cov_fifth = np.array(cov_fifth)

#MAKE TABLE
#REAL
real_mean, real_mean3, real_mean5 = np.nanmean(real_first), np.nanmean(real_third), np.nanmean(real_fifth)
real_16, real_50, real_84 = np.nanpercentile(real_first, 16), np.nanpercentile(real_first, 50), np.nanpercentile(real_first, 84)
real3_16, real3_50, real3_84 = np.nanpercentile(real_third, 16), np.nanpercentile(real_third, 50), np.nanpercentile(real_third, 84)
real5_16, real5_50, real5_84 = np.nanpercentile(real_fifth, 16), np.nanpercentile(real_fifth, 50), np.nanpercentile(real_fifth, 84)

#EXP
exp_mean, exp_mean3, exp_mean5 = np.nanmean(exp_first), np.nanmean(exp_third), np.nanmean(exp_fifth)
exp_16, exp_50, exp_84 = np.nanpercentile(exp_first, 16), np.nanpercentile(exp_first, 50), np.nanpercentile(exp_first, 84)
exp3_16, exp3_50, exp3_84 = np.nanpercentile(exp_third, 16), np.nanpercentile(exp_third, 50), np.nanpercentile(exp_third, 84)
exp5_16, exp5_50, exp5_84 = np.nanpercentile(exp_fifth, 16), np.nanpercentile(exp_fifth, 50), np.nanpercentile(exp_fifth, 84)

#COV
cov_mean, cov_mean3, cov_mean5 = np.nanmean(cov_first), np.nanmean(cov_third), np.nanmean(cov_fifth)
cov_16, cov_50, cov_84 = np.nanpercentile(cov_first, 16), np.nanpercentile(cov_first, 50), np.nanpercentile(cov_first, 84)
cov3_16, cov3_50, cov3_84 = np.nanpercentile(cov_third, 16), np.nanpercentile(cov_third, 50), np.nanpercentile(cov_third, 84)
cov5_16, cov5_50, cov5_84 = np.nanpercentile(cov_fifth, 16), np.nanpercentile(cov_fifth, 50), np.nanpercentile(cov_fifth, 84)


first_mean = np.array([real_mean, exp_mean, cov_mean]).astype(int)
first_16 = np.array([real_16, exp_16, cov_16]).astype(int)
first_50 = np.array([real_50, exp_50, cov_50]).astype(int)
first_84 = np.array([real_84, exp_84, cov_84]).astype(int)

third_mean = np.array([real_mean3, exp_mean3, cov_mean]).astype(int)
third_16 = np.array([real3_16, exp3_16, cov3_16]).astype(int)
third_50 = np.array([real3_50, exp3_50, cov3_50]).astype(int)
third_84 = np.array([real3_84, exp3_84, cov3_84]).astype(int)

fifth_mean = np.array([real_mean5, exp_mean5, cov_mean5]).astype(int)
fifth_16 = np.array([real5_16, exp5_16, cov5_16]).astype(int)
fifth_50 = np.array([real5_50, exp5_50, cov5_50]).astype(int)
fifth_84 = np.array([real5_84, exp5_84, cov5_84]).astype(int)

catalogs = np.array(['Real', 'Exponential', 'Flat'])

columns = [catalogs, first_mean, first_16, first_50, first_84, third_mean, third_16, third_50, third_84,
                        fifth_mean, fifth_16, fifth_50, fifth_84]
colnames = ('Catalog','1st: Mean', '1st: 16th', '1st: 50th', '1st: 84th', '3rd: Mean', '3rd: 16th', '3rd: 50th', '3rd: 84th', 'Fifth: Mean', 'Fifth: 16th', 'Fifth: 50th', 'Fifth: 84th')
units = [u.pc, u.pc, u.pc, u.pc, u.pc, u.pc, u.pc, u.pc, u.pc, u.pc, u.pc, u.pc]

spacing_tab = Table(columns, names=colnames)
spacing_tab.write(home_dir + 'gmc_spacing/tables/spacing_tab.tex', format='latex', overwrite=True)


### MAKE PLOTS ###
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.size': 12})
fig, ax = plt.subplots(3,1, figsize=(5,15), dpi=155)

bins = np.linspace(0,2000,41) #50pc bins
ax[0].hist(real_first, bins=bins, density=True, lw=2, histtype='step',
          facecolor='grey', hatch='/', edgecolor='black', alpha=0.7,
           fill=True, label='Real GMCs', zorder=10)
ax[0].hist(exp_first, bins=bins, density=True, lw=2, histtype='step',
          facecolor='None', edgecolor='maroon', alpha=1,
           fill=False, label='Random (Exponential)', zorder=12, ls='--')
ax[0].hist(cov_first, bins=bins, density=True, lw=2, histtype='step',
          facecolor='None', edgecolor='navy', alpha=1,
           fill=False, label='Random (Flat)', zorder=11)
ax[0].legend()
ax[0].set_ylabel('Distance to Nearest Neighbor (pc)')
ax[0].set_xlabel('Probability Density')
ax[0].tick_params(direction='in')
ax[0].grid(alpha=1, ls=':', axis='both', zorder=0)

ax[1].hist(real_third, bins=bins, density=True, lw=2, histtype='step',
          facecolor='grey', hatch='/', edgecolor='black', alpha=0.7,
           fill=True, label='Real GMCs', zorder=10)
ax[1].hist(exp_third, bins=bins, density=True, lw=2, histtype='step',
          facecolor='None', edgecolor='maroon', alpha=1,
           fill=False, label='Random (Exponential)', zorder=12, ls='--')
ax[1].hist(cov_third, bins=bins, density=True, lw=2, histtype='step',
          facecolor='None', edgecolor='navy', alpha=1,
           fill=False, label='Random (Flat)', zorder=11)
ax[1].set_ylabel('Distance to 3rd Nearest Neighbor (pc)')
ax[1].set_xlabel('Probability Density')
ax[1].tick_params(direction='in')
ax[1].grid(alpha=1, ls=':', axis='both', zorder=0)

ax[2].hist(real_fifth, bins=bins, density=True, lw=2, histtype='step',
          facecolor='grey', hatch='/', edgecolor='black', alpha=0.7,
           fill=True, label='Real GMCs', zorder=10)
ax[2].hist(exp_fifth, bins=bins, density=True, lw=2, histtype='step',
          facecolor='None', edgecolor='maroon', alpha=1,
           fill=False, label='Random (Exponential)', zorder=12, ls='--')
ax[2].hist(cov_fifth, bins=bins, density=True, lw=2, histtype='step',
          facecolor='None', edgecolor='navy', alpha=1,
           fill=False, label='Random (Flat)', zorder=11)
ax[2].set_ylabel('Distance to 5th Nearest Neighbor (pc)')
ax[2].set_xlabel('Probability Density')
ax[2].tick_params(direction='in')
ax[2].grid(alpha=1, ls=':', axis='both', zorder=0)


plt.savefig(home_dir + 'gmc_spacing/figures/spacing_histogram_3panel.pdf', bbox_inches='tight')