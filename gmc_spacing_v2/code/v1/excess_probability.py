## CALCULATES EXCESS PROBABILITY AGAINST CONTROL DISTRIBUTIONS ##

home_dir = '/home/machado.35/projects/' #Set this as the location of where the gmc_spacing folder exists

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import pandas as pd
from os import listdir
import os
import matplotlib.patches as mpatches
import matplotlib.lines as mlines



## Functions necessary for excess probability calcuations -- mainly data cleaning ##

def pull_med(x):
    z = np.sort(x)
    z = np.trim_zeros(z)
    if len(z) == 0:
        return 0
    else:
        return np.nanmedian(z)

def pull_low(x):
    z = np.sort(x)
    z = np.trim_zeros(z)
    if len(z) == 0:
        return 0
    else:
        return np.nanpercentile(z, 16)
    
def pull_high(x):
    z = np.sort(x)
    z = np.trim_zeros(z)
    if len(z) == 0:
        return 0
    else:
        return np.nanpercentile(z, 84)
    

def trim(x):
    z = np.sort(x)
    if np.sum(z == 0) == len(z):
        return(z)
    else:
        z = np.trim_zeros(z)
        return z
    
    

## Plotting pre processing ##
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.size': 12})


sources = listdir(home_dir+'gmc_spacing/data/cloud_cats_csv/')
for i in range(len(sources)):
    sources[i] = sources[i].split('_')[0]

#Grab all 1st, 3rd, 5th real GMCs
first_real, third_real, fifth_real = [], [], []

for i in range(len(sources)):
    
    cat = pd.read_csv(home_dir+'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv')
    first_real.extend(cat['nn_dist'])
    third_real.extend(cat['nn3_dist'])
    fifth_real.extend(cat['nn5_dist'])
    
first_real = np.array(first_real)
third_real = np.array(third_real)
fifth_real = np.array(fifth_real)

med_real_first = pull_med(first_real)
med_real_third = pull_med(third_real)
med_real_fifth = pull_med(fifth_real)


#Grab all 1st, 3rd, 5th FLAT RANDOM GMCs
first_flat, third_flat, fifth_flat = [], [], []

for i in range(len(sources)):
    name = sources[i].split('_')[0]
    data = home_dir+'gmc_spacing/data/sim_results/flat/'+sources[i]+'_flat_150pc.csv'
    if os.path.isfile(data):
        cat = pd.read_csv(data)
        first_flat.extend(cat['nn_dist'])
        third_flat.extend(cat['nn3_dist'])
        fifth_flat.extend(cat['nn5_dist'])

first_flat = np.array(first_flat)
third_flat = np.array(third_flat)
fifth_flat = np.array(fifth_flat)


#Grab all 1st, 3rd, 5th EXPONENTIAL RANDOM GMCs
first_exp, third_exp, fifth_exp = [], [], []

for i in range(len(sources)):
    name = sources[i].split('_')[0]
    data = home_dir+'gmc_spacing/data/sim_results/exp/'+sources[i]+'_exp_150pc.csv'
    if os.path.isfile(data):
        cat = pd.read_csv(data)
        first_exp.extend(cat['nn_dist'])
        third_exp.extend(cat['nn3_dist'])
        fifth_exp.extend(cat['nn5_dist'])

first_exp = np.array(first_exp)
third_exp = np.array(third_exp)
fifth_exp = np.array(fifth_exp)


#Calculate excess probability curves for all galaxies
bins = np.linspace(0,2100,15)
binsize = bins[1]-bins[0]
cov_medians = np.zeros(shape=(3,len(bins),len(sources)))
exp_medians = np.zeros(shape=(3,len(bins),len(sources)))


for i in range(len(sources)):
    name = sources[i].split('_')[0]
    data = home_dir+'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
    if os.path.isfile(data):
        flat_cat = home_dir+'gmc_spacing/data/sim_results/flat/'+sources[i]+'_flat_150pc.csv'
        if os.path.isfile(flat_cat):
            exp_cat = home_dir+'gmc_spacing/data/sim_results/exp/'+sources[i]+'_exp_150pc.csv'
            flat_cat = pd.read_csv(flat_cat)
            exp_cat = pd.read_csv(exp_cat)
            real_cat = pd.read_csv(data)
            
            ####################################
            ###Retrieve Nth nearest distances###
            ####################################
            
            real_first = real_cat['nn_dist']
            real_third = real_cat['nn3_dist']
            real_fifth = real_cat['nn5_dist']
            
            flat_first = flat_cat['nn_dist']
            flat_third = flat_cat['nn3_dist']
            flat_fifth = flat_cat['nn5_dist']
            
            exp_first = exp_cat['nn_dist']
            exp_third = exp_cat['nn3_dist']
            exp_fifth = exp_cat['nn5_dist']
            
            
            
            ##############################
            ###Normalize counts per bin###
            ##############################
            true_first = np.histogram(real_first, bins=bins, density=True)
            true_third = np.histogram(real_third, bins=bins, density=True)
            true_fifth = np.histogram(real_fifth, bins=bins, density=True)
            
            rand_flat_first = np.histogram(flat_first, bins=bins, density=True)
            rand_flat_third = np.histogram(flat_third, bins=bins, density=True)
            rand_flat_fifth = np.histogram(flat_fifth, bins=bins, density=True)
            
            rand_exp_first = np.histogram(exp_first, bins=bins, density=True)
            rand_exp_third = np.histogram(exp_third, bins=bins, density=True)
            rand_exp_fifth = np.histogram(exp_fifth, bins=bins, density=True)
            
            
            ############################
            ###Calculate excess probs###
            ############################
            
            #flat catalog
            excess_flat_first = true_first[0]/rand_flat_first[0]
            excess_flat_third = true_third[0]/rand_flat_third[0]
            excess_flat_fifth = true_fifth[0]/rand_flat_fifth[0]
            
            #exponential catalog
            excess_exp_first = true_first[0]/rand_exp_first[0]
            excess_exp_third = true_third[0]/rand_exp_third[0]
            excess_exp_fifth = true_fifth[0]/rand_exp_fifth[0]
            
            
            #############################
            ###Remove Nan and Inf Vals###
            #############################
            
            bad_vals_flat_first = np.where(np.isfinite(excess_flat_first) == False)
            bad_vals_flat_third = np.where(np.isfinite(excess_flat_third) == False)
            bad_vals_flat_fifth = np.where(np.isfinite(excess_flat_fifth) == False)
            
            excess_flat_first[bad_vals_flat_first] = 0
            excess_flat_third[bad_vals_flat_third] = 0
            excess_flat_fifth[bad_vals_flat_fifth] = 0
            
            bad_vals_exp_first = np.where(np.isfinite(excess_exp_first) == False)
            bad_vals_exp_third = np.where(np.isfinite(excess_exp_third) == False)
            bad_vals_exp_fifth = np.where(np.isfinite(excess_exp_fifth) == False)
            
            excess_exp_first[bad_vals_exp_first] = 0
            excess_exp_third[bad_vals_exp_third] = 0
            excess_exp_fifth[bad_vals_exp_third] = 0
            
            
            ###########################################################
            ###Find median excess prob (across galaxies) in each bin###
            ###########################################################
            for k in range(len(bins[1:])):
                
                cov_medians[0][k][i] = excess_flat_first[k]
                cov_medians[1][k][i] = excess_flat_third[k]
                cov_medians[2][k][i] = excess_flat_fifth[k]
                
                exp_medians[0][k][i] = excess_exp_first[k]
                exp_medians[1][k][i] = excess_exp_third[k]
                exp_medians[2][k][i] = excess_exp_fifth[k]
                
                
                
flat_meds_first, flat_meds_third, flat_meds_fifth = [], [], []
flat_lows_first, flat_lows_third, flat_lows_fifth = [], [], []
flat_highs_first, flat_highs_third, flat_highs_fifth = [], [], []

exp_meds_first, exp_meds_third, exp_meds_fifth = [], [], []
exp_lows_first, exp_lows_third, exp_lows_fifth = [], [], []
exp_highs_first, exp_highs_third, exp_highs_fifth = [], [], []


### First Nearest Neighbor ###
for i in range(len(cov_medians[0][:-1])):
    flat_meds_first.append(pull_med(cov_medians[0][i]))
    flat_lows_first.append(pull_low(cov_medians[0][i]))
    flat_highs_first.append(pull_high(cov_medians[0][i]))

for i in range(len(exp_medians[0][:-1])):
    exp_meds_first.append(pull_med(exp_medians[0][i]))
    exp_lows_first.append(pull_low(exp_medians[0][i]))
    exp_highs_first.append(pull_high(exp_medians[0][i]))



### Third Nearest Neighbor ###
for i in range(len(cov_medians[0][:-1])):
    flat_meds_third.append(pull_med(cov_medians[1][i]))
    flat_lows_third.append(pull_low(cov_medians[1][i]))
    flat_highs_third.append(pull_high(cov_medians[1][i]))

for i in range(len(exp_medians[0][:-1])):
    exp_meds_third.append(pull_med(exp_medians[1][i]))
    exp_lows_third.append(pull_low(exp_medians[1][i]))
    exp_highs_third.append(pull_high(exp_medians[1][i]))

    
    
### Fifth Nearest Neighbor ###
for i in range(len(cov_medians[0][:-1])):
    flat_meds_fifth.append(pull_med(cov_medians[2][i]))
    flat_lows_fifth.append(pull_low(cov_medians[2][i]))
    flat_highs_fifth.append(pull_high(cov_medians[2][i]))

for i in range(len(exp_medians[0][:-1])):
    exp_meds_fifth.append(pull_med(exp_medians[2][i]))
    exp_lows_fifth.append(pull_low(exp_medians[2][i]))
    exp_highs_fifth.append(pull_high(exp_medians[2][i]))
    
    
    
#Real CDFs
real_sort = np.sort(first_real)
real_cdf = np.arange(len(first_real)) / float(len(first_real))

real_sort3 = np.sort(third_real)
real_cdf3 = np.arange(len(third_real)) / float(len(third_real))

real_sort5 = np.sort(fifth_real)
real_cdf5 = np.arange(len(fifth_real)) / float(len(fifth_real))


#Exponential CDF
flat_sort = np.sort(first_flat)
flat_cdf = np.arange(len(first_flat)) / float(len(first_flat))

flat_sort3 = np.sort(third_flat)
flat_cdf3 = np.arange(len(third_flat)) / float(len(third_flat))

flat_sort5 = np.sort(fifth_flat)
flat_cdf5 = np.arange(len(fifth_flat)) / float(len(fifth_flat))



## GRIDDING ##
xvals = np.linspace(0,2000,501)
real_gridded = np.interp(xvals, real_sort, real_cdf)
real_gridded3 = np.interp(xvals, real_sort3, real_cdf3)
real_gridded5 = np.interp(xvals, real_sort5, real_cdf5)

flat_gridded = np.interp(xvals, flat_sort, flat_cdf)
flat_gridded3 = np.interp(xvals, flat_sort3, flat_cdf3)
flat_gridded5 = np.interp(xvals, flat_sort5, flat_cdf5)


# Free up memory #
del first_flat, third_flat, fifth_flat






##################################################
### EXCESS PROBABILITY CURVES FOR FLAT CATALOG ###
##################################################
fig, ax = plt.subplots(2,2, figsize=(12,10), dpi=150)



#TOP LEFT
pan_colors = ['mediumseagreen', 'teal', 'darkslategrey']

ax[0,0].axvspan(np.nanpercentile(first_real, 16), np.nanpercentile(first_real, 84), alpha=.3, color='gray', hatch='//', zorder=1)
ax[0,0].axvspan(np.nanpercentile(first_real, 16), np.nanpercentile(first_real, 84), alpha=.3, color='None', edgecolor='gray', hatch='//', zorder=1)
ax[0,0].axvspan(np.nanpercentile(first_real, 5), np.nanpercentile(first_real, 95), alpha=0.15, color='grey')
ax[0,0].vlines(med_real_first, 0.01, 1000, color='black', ls=':', lw=2, alpha=0.7)
ax[0,0].plot(bins[1:], flat_meds_first, label='Real / Flat', c=pan_colors[0], lw=3, ls='dashdot')
ax[0,0].fill_between(bins[1:], flat_lows_first, flat_highs_first, color=pan_colors[0], edgecolor='None', alpha=0.4)

#ax[0,0].plot(xvals, real_gridded/flat_gridded, c='maroon', ls='--', lw=3, label="CDF Ratio (Real / Flat)")

ax[0,0].set_yscale('log')
ax[0,0].set_ylim(0.1, 10)
ax[0,0].set_xlim(100,1700)
#plt.xlim(np.nanpercentile(real_gmcs, 5), np.nanpercentile(real_gmcs, 95))
ax[0,0].legend()
ax[0,0].tick_params(direction='in', which='both')
ax[0,0].axhline(y=1, color='black', linestyle='--', linewidth=3, label='Equal Prob.', zorder=5)
ax[0,0].set_ylabel('Probability Ratio')
ax[0,0].set_xlabel('Distance to 1st Nearest Neighbor (pc)')


#TOP RIGHT
ax[0,1].axvspan(np.nanpercentile(third_real, 16), np.nanpercentile(third_real, 84), alpha=.3, hatch='//', color='gray', zorder=1)
ax[0,1].axvspan(np.nanpercentile(third_real, 16), np.nanpercentile(third_real, 84), alpha=.3, hatch='//', color='None', edgecolor='gray', zorder=1)
ax[0,1].axvspan(np.nanpercentile(third_real, 5), np.nanpercentile(third_real, 95), alpha=0.15, color='grey')
ax[0,1].vlines(med_real_third, 0.01, 1000, color='black', ls=':', lw=2, alpha=0.7)
ax[0,1].plot(bins[1:], flat_meds_third, label='Real / Flat', c=pan_colors[1], lw=3, ls='--')
ax[0,1].fill_between(bins[1:], flat_lows_third, flat_highs_third, color=pan_colors[1], edgecolor='None', alpha=0.2)

#ax[0,1].plot(xvals, real_gridded3/flat_gridded3, c='maroon', ls='--', lw=3, label="CDF Ratio (Real / Flat)")

ax[0,1].set_yscale('log')
ax[0,1].set_ylim(0.1, 10)
#plt.xlim(np.nanpercentile(real_gmcs, 5), np.nanpercentile(real_gmcs, 95))
ax[0,1].legend()
ax[0,1].tick_params(direction='in', which='both')
ax[0,1].set_xlim(100,1700)
ax[0,1].axhline(y=1, color='black', linestyle='--', linewidth=3, label='Equal Prob.', zorder=5)
ax[0,1].set_ylabel('Probability Ratio')
ax[0,1].set_xlabel('Distance to 3rd Nearest Neighbor (pc)')



#BOTTOM LEFT
ax[1,0].axvspan(np.nanpercentile(fifth_real, 16), np.nanpercentile(fifth_real, 84), alpha=.3, hatch='//', color='gray', zorder=1)
ax[1,0].axvspan(np.nanpercentile(fifth_real, 16), np.nanpercentile(fifth_real, 84), alpha=.3, hatch='//', color='None', edgecolor='gray', zorder=1)
ax[1,0].axvspan(np.nanpercentile(fifth_real, 5), np.nanpercentile(fifth_real, 95), alpha=0.15, color='grey', zorder=0)
ax[1,0].vlines(med_real_fifth, 0.01, 1000, color='black', ls=':', lw=2, alpha=0.7)
ax[1,0].plot(bins[1:], flat_meds_fifth, label='Real / Flat', c=pan_colors[2], lw=3)
ax[1,0].fill_between(bins[1:], flat_lows_fifth, flat_highs_fifth, color=pan_colors[2], edgecolor='None', alpha=0.4, zorder=2)

#ax[1,0].plot(xvals, real_gridded5/flat_gridded5, c='maroon', ls='--', lw=3, label="CDF Ratio (Real / Flat)")

ax[1,0].set_yscale('log')
ax[1,0].set_ylim(0.1, 100)
ax[1,0].set_xlim(100,1700)
#plt.xlim(np.nanpercentile(real_gmcs, 5), np.nanpercentile(real_gmcs, 95))
ax[1,0].legend()
ax[1,0].tick_params(direction='in', which='both')
ax[1,0].axhline(y=1, color='black', linestyle='--', linewidth=3, label='Equal Prob.', zorder=5)
ax[1,0].set_ylabel('Probability Ratio')
ax[1,0].set_xlabel('Distance to 5th Nearest Neighbor (pc)')



#BOTTOM RIGHT
colors = ['mediumseagreen', 'teal', 'darkslategrey']

ax[1,1].plot(bins[1:], flat_meds_first, label='1st', c=colors[0], ls='dashdot', lw=3)
ax[1,1].vlines(med_real_first, 0.01,1e4, color=colors[0], lw=2, alpha=0.5, zorder=0, ls=':')
ax[1,1].plot(bins[1:], flat_meds_third, label='3rd', c=colors[1], ls='--', lw=3)
ax[1,1].vlines(med_real_third, 0.01,1e4, color=colors[1], lw=2, alpha=0.5, zorder=0, ls=':')
ax[1,1].plot(bins[1:], flat_meds_fifth, label='5th', c=colors[2], lw=3)
ax[1,1].vlines(med_real_fifth, 0.01,1e4, color=colors[2], lw=2, alpha=0.5, zorder=0, ls=':')

ax[1,1].scatter(med_real_first, 0.2, color=colors[0])
xerr = [[med_real_first - np.nanpercentile(first_real, 16)], [np.nanpercentile(first_real, 84) - med_real_first]]
ax[1,1].errorbar(med_real_first, 0.2, yerr=None, xerr=xerr, color=colors[0], capsize=3, zorder=5)

ax[1,1].scatter(med_real_third, 0.2333, color=colors[1])
xerr = [[med_real_third - np.nanpercentile(third_real, 16)], [np.nanpercentile(third_real, 84)] - med_real_third]
ax[1,1].errorbar(med_real_third, 0.2333, yerr=None, xerr=xerr, color=colors[1], capsize=3, zorder=5)

ax[1,1].scatter(med_real_fifth, 0.3, color=colors[2])
xerr = [[med_real_fifth - np.nanpercentile(fifth_real, 16)], [np.nanpercentile(fifth_real, 84)] - med_real_fifth]
ax[1,1].errorbar(med_real_fifth, 0.3, yerr=None, xerr=xerr, color=colors[2], capsize=3, zorder=5)

ax[1,1].set_yscale('log')
ax[1,1].set_ylim(0.1, 100)
ax[1,1].set_xlim(100,1700)
#plt.xlim(np.nanpercentile(real_gmcs, 5), np.nanpercentile(real_gmcs, 95))
ax[1,1].legend(loc='lower right')
ax[1,1].tick_params(direction='in', which='both')
ax[1,1].axhline(y=1, color='black', linestyle='--', linewidth=3, label='Equal Prob.', zorder=5)
ax[1,1].set_ylabel('Probability Ratio')
ax[1,1].set_xlabel('Distance to Nth Nearest Neighbor (pc)')


plt.savefig(home_dir+'gmc_spacing/figures/four_panel_excess_prob_flat_150pc_bins.pdf', bbox_inches='tight')





### CDF RATIOS FOR EXPONENTIAL COMPARISON ###

#Real CDFs
real_sort = np.sort(first_real)
real_cdf = np.arange(len(first_real)) / float(len(first_real))

real_sort3 = np.sort(third_real)
real_cdf3 = np.arange(len(third_real)) / float(len(third_real))

real_sort5 = np.sort(fifth_real)
real_cdf5 = np.arange(len(fifth_real)) / float(len(fifth_real))


#Exponential CDF
exp_sort = np.sort(first_exp)
exp_cdf = np.arange(len(first_exp)) / float(len(first_exp))

exp_sort3 = np.sort(third_exp)
exp_cdf3 = np.arange(len(third_exp)) / float(len(third_exp))

exp_sort5 = np.sort(fifth_exp)
exp_cdf5 = np.arange(len(fifth_exp)) / float(len(fifth_exp))

## GRIDDING ##
xvals = np.linspace(0,2000,501)
real_gridded = np.interp(xvals, real_sort, real_cdf)
real_gridded3 = np.interp(xvals, real_sort3, real_cdf3)
real_gridded5 = np.interp(xvals, real_sort5, real_cdf5)

exp_gridded = np.interp(xvals, exp_sort, exp_cdf)
exp_gridded3 = np.interp(xvals, exp_sort3, exp_cdf3)
exp_gridded5 = np.interp(xvals, exp_sort5, exp_cdf5)

# Free up memory #
del first_exp, third_exp, fifth_exp





#########################################################
### EXCESS PROBABILITY CURVES FOR EXPONENTIAL CATALOG ###
#########################################################
fig, ax = plt.subplots(2,2, figsize=(12,10), dpi=150)



#TOP LEFT
pan_colors = ['rosybrown', 'indianred', 'darkred']

ax[0,0].axvspan(np.nanpercentile(first_real, 16), np.nanpercentile(first_real, 84), alpha=.3, color='gray', hatch='//', zorder=1)
ax[0,0].axvspan(np.nanpercentile(first_real, 16), np.nanpercentile(first_real, 84), alpha=.3, color='None', edgecolor='gray', hatch='//', zorder=1)
ax[0,0].axvspan(np.nanpercentile(first_real, 5), np.nanpercentile(first_real, 95), alpha=0.15, color='grey')
ax[0,0].vlines(med_real_first, 0.01, 1000, color='black', ls=':', lw=2, alpha=0.7)

ax[0,0].plot(bins[1:], exp_meds_first, label='Real / Exponential', c=pan_colors[0], lw=3)
ax[0,0].fill_between(bins[1:], exp_lows_first, exp_highs_first, color=pan_colors[0], edgecolor='None', alpha=0.4)

#ax[0,0].plot(xvals, real_gridded/exp_gridded, c='royalblue', ls='--', lw=3, label="CDF Ratio (Real / Exp)")

ax[0,0].set_yscale('log')
ax[0,0].set_ylim(0.1, 10)
ax[0,0].set_xlim(100,1700)
#plt.xlim(np.nanpercentile(real_gmcs, 5), np.nanpercentile(real_gmcs, 95))
ax[0,0].legend()
ax[0,0].tick_params(direction='in', which='both')
ax[0,0].axhline(y=1, color='black', linestyle='--', linewidth=3, label='Equal Prob.', zorder=5)
ax[0,0].set_ylabel('Probability Ratio')
ax[0,0].set_xlabel('Distance to 1st Nearest Neighbor (pc)')


#TOP RIGHT
ax[0,1].axvspan(np.nanpercentile(third_real, 16), np.nanpercentile(third_real, 84), alpha=.3, hatch='//', color='gray', zorder=1)
ax[0,1].axvspan(np.nanpercentile(third_real, 16), np.nanpercentile(third_real, 84), alpha=.3, hatch='//', color='None', edgecolor='gray', zorder=1)
ax[0,1].axvspan(np.nanpercentile(third_real, 5), np.nanpercentile(third_real, 95), alpha=0.15, color='grey')
ax[0,1].vlines(med_real_third, 0.01, 1000, color='black', ls=':', lw=2, alpha=0.7)

ax[0,1].plot(bins[1:], exp_meds_third, label='Real / Exponential', c=pan_colors[1], lw=3)
ax[0,1].fill_between(bins[1:], exp_lows_third, exp_highs_third, color=pan_colors[1], edgecolor='None', alpha=0.4)

#ax[0,1].plot(xvals, real_gridded3/exp_gridded3, c='royalblue', ls='--', lw=3, label="CDF Ratio (Real / Exp)")

ax[0,1].set_yscale('log')
ax[0,1].set_ylim(0.1, 10)
#plt.xlim(np.nanpercentile(real_gmcs, 5), np.nanpercentile(real_gmcs, 95))
ax[0,1].legend()
ax[0,1].tick_params(direction='in', which='both')
ax[0,1].set_xlim(100,1700)
ax[0,1].axhline(y=1, color='black', linestyle='--', linewidth=3, label='Equal Prob.', zorder=5)
ax[0,1].set_ylabel('Excess Probability')
ax[0,1].set_xlabel('Distance to 3rd Nearest Neighbor (pc)')



#BOTTOM LEFT
ax[1,0].axvspan(np.nanpercentile(fifth_real, 16), np.nanpercentile(fifth_real, 84), alpha=.3, hatch='//', color='gray', zorder=1)
ax[1,0].axvspan(np.nanpercentile(fifth_real, 16), np.nanpercentile(fifth_real, 84), alpha=.3, hatch='//', color='None', edgecolor='gray', zorder=1)
ax[1,0].axvspan(np.nanpercentile(fifth_real, 5), np.nanpercentile(fifth_real, 95), alpha=0.15, color='grey', zorder=0)
ax[1,0].vlines(med_real_fifth, 0.01, 1000, color='black', ls=':', lw=2, alpha=0.7)

ax[1,0].plot(bins[1:], exp_meds_fifth, label='Real / Exponential', c=pan_colors[2], lw=3)
ax[1,0].fill_between(bins[1:], exp_lows_fifth, exp_highs_fifth, color=pan_colors[2], edgecolor='None', alpha=0.4, zorder=2)

#ax[1,0].plot(xvals, real_gridded5/exp_gridded5, c='royalblue', ls='--', lw=3, label="CDF Ratio (Real / Exp)")

ax[1,0].set_yscale('log')
ax[1,0].set_ylim(0.1, 10)
ax[1,0].set_xlim(100,1700)
#plt.xlim(np.nanpercentile(real_gmcs, 5), np.nanpercentile(real_gmcs, 95))
ax[1,0].legend()
ax[1,0].tick_params(direction='in', which='both')
ax[1,0].axhline(y=1, color='black', linestyle='--', linewidth=3, label='Equal Prob.', zorder=5)
ax[1,0].set_ylabel('Probability Ratio')
ax[1,0].set_xlabel('Distance to 5th Nearest Neighbor (pc)')



#BOTTOM RIGHT
colors = ['rosybrown', 'indianred', 'darkred']

ax[1,1].plot(bins[1:], exp_meds_first, label='1st', c=colors[0], ls='dashdot', lw=3)
ax[1,1].vlines(med_real_first, 0.01,1e4, color=colors[0], lw=2, alpha=0.5, zorder=0, ls=':')
ax[1,1].plot(bins[1:], exp_meds_third, label='3rd', c=colors[1], ls='--', lw=3)
ax[1,1].vlines(med_real_third, 0.01,1e4, color=colors[1], lw=2, alpha=0.5, zorder=0, ls=':')
ax[1,1].plot(bins[1:], exp_meds_fifth, label='5th', c=colors[2], lw=3)
ax[1,1].vlines(med_real_fifth, 0.01,1e4, color=colors[2], lw=2, alpha=0.5, zorder=0, ls=':')

ax[1,1].scatter(med_real_first, 0.2, color=colors[0])
xerr = [[med_real_first - np.nanpercentile(first_real, 16)], [np.nanpercentile(first_real, 84) - med_real_first]]
ax[1,1].errorbar(med_real_first, 0.2, yerr=None, xerr=xerr, color=colors[0], capsize=3, zorder=5)

ax[1,1].scatter(med_real_third, 0.2333, color=colors[1])
xerr = [[med_real_third - np.nanpercentile(third_real, 16)], [np.nanpercentile(third_real, 84)] - med_real_third]
ax[1,1].errorbar(med_real_third, 0.2333, yerr=None, xerr=xerr, color=colors[1], capsize=3, zorder=5)

ax[1,1].scatter(med_real_fifth, 0.3, color=colors[2])
xerr = [[med_real_fifth - np.nanpercentile(fifth_real, 16)], [np.nanpercentile(fifth_real, 84)] - med_real_fifth]
ax[1,1].errorbar(med_real_fifth, 0.3, yerr=None, xerr=xerr, color=colors[2], capsize=3, zorder=5)

ax[1,1].set_yscale('log')
ax[1,1].set_ylim(0.1, 10)
ax[1,1].set_xlim(100,1700)
#plt.xlim(np.nanpercentile(real_gmcs, 5), np.nanpercentile(real_gmcs, 95))
ax[1,1].legend(loc='lower right')
ax[1,1].tick_params(direction='in', which='both')
ax[1,1].axhline(y=1, color='black', linestyle='--', linewidth=3, label='Equal Prob.', zorder=5)
ax[1,1].set_ylabel('Probability Ratio')
ax[1,1].set_xlabel('Distance to Nth Nearest Neighbor (pc)')


plt.savefig(home_dir+'gmc_spacing/figures/four_panel_excess_prob_exp_150pc_bins.pdf', bbox_inches='tight')
