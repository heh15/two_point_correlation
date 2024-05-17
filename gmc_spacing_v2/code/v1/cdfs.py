### THIS SCRIPT GENERATES CDFS FOR ALL CATALOGS ###


home_dir = '/home/machado.35/projects/'


from os import listdir
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


sources = listdir(home_dir + 'gmc_spacing/data/cloud_cats_csv/')
for i in range(len(sources)):
    sources[i] = sources[i].split('_')[0]
    

    
#1st Nearest Neighbor    
real_clouds = []
exp_clouds = []
cov_clouds = []
for i in range(len(sources)):
    name = sources[i].split('_')[0]
    clouds = home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
    rand_cov = home_dir + 'gmc_spacing/data/sim_results/flat/'+sources[i]+'_flat_150pc.csv'
    rand_exp = home_dir + 'gmc_spacing/data/sim_results/exp/'+sources[i]+'_exp_150pc.csv'
    if os.path.isfile(clouds) and os.path.isfile(rand_cov) and os.path.isfile(rand_exp):
        real_clouds.extend(pd.read_csv(clouds)['nn_dist'])
        cov_clouds.extend(pd.read_csv(rand_cov)['nn_dist'])
        exp_clouds.extend(pd.read_csv(rand_exp)['nn_dist'])
        
        
## HIGH MASS LOW MASS ##
masses = []
high_mass = []
low_mass = []
for i in range(len(sources)):
    cat = pd.read_csv(home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv')
    if len(cat) > 5:
        masses.extend(cat['MLUM_MSUN'])
mass_div = np.nanmedian(masses)
for i in range(len(sources)):
    cat = pd.read_csv(home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv')
    for k in range(len(cat)):
        if len(cat) > 5:
            if cat['MLUM_MSUN'][k] > mass_div:
                high_mass.extend([cat['nn_dist'][k]])
            else:
                low_mass.extend([cat['nn_dist'][k]])

#Exponential CDF
exp_sort = np.sort(exp_clouds)
exp_cdf = np.arange(len(exp_clouds))/ float(len(exp_clouds))

#2D Coverage CDF
cov_sort = np.sort(cov_clouds)
cov_cdf = np.arange(len(cov_clouds))/float(len(cov_clouds))

#Real Data CDF
real_sort = np.sort(real_clouds)
real_cdf = np.arange(len(real_clouds))/float(len(real_clouds))

#High Mass CDF
hm_sort = np.sort(high_mass)
hm_cdf = np.arange(len(high_mass))/float(len(high_mass))

#Low Mass CDF
lm_sort = np.sort(low_mass)
lm_cdf = np.arange(len(low_mass))/float(len(low_mass))

## GRIDDING ##
xvals = np.linspace(0,2000,501)
real_gridded = np.interp(xvals, real_sort, real_cdf)
cov_gridded = np.interp(xvals, cov_sort, cov_cdf)
exp_gridded = np.interp(xvals, exp_sort, exp_cdf)
hm_gridded = np.interp(xvals, hm_sort, hm_cdf)
lm_gridded = np.interp(xvals, lm_sort, lm_cdf)

# Free up memory #
del exp_clouds, exp_sort, exp_cdf
del real_clouds, real_sort, real_cdf
del cov_clouds, cov_sort, cov_cdf
del high_mass, hm_sort, hm_cdf
del low_mass, lm_sort, lm_cdf




## 3rd Nearest Neighbor
real_clouds3 = []
exp_clouds3 = []
cov_clouds3 = []
for i in range(len(sources)):
    name = sources[i].split('_')[0]
    clouds = home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
    rand_cov = home_dir + 'gmc_spacing/data/sim_results/flat/'+sources[i]+'_flat_150pc.csv'
    rand_exp = home_dir + 'gmc_spacing/data/sim_results/exp/'+sources[i]+'_exp_150pc.csv'
    if os.path.isfile(clouds) and os.path.isfile(rand_cov) and os.path.isfile(rand_exp):
        real_clouds3.extend(pd.read_csv(clouds)['nn3_dist'])
        cov_clouds3.extend(pd.read_csv(rand_cov)['nn3_dist'])
        exp_clouds3.extend(pd.read_csv(rand_exp)['nn3_dist'])

## HIGH MASS LOW MASS ##
masses3 = []
high_mass3 = []
low_mass3 = []
for i in range(len(sources)):
    cat = pd.read_csv(home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv')
    masses3.extend(cat['MLUM_MSUN'])
mass_div3 = np.nanmedian(masses3)
for i in range(len(sources)):
    cat = pd.read_csv(home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv')
    for k in range(len(cat)):
        if cat['MLUM_MSUN'][k] > mass_div3:
            high_mass3.extend([cat['nn3_dist'][k]])
        else:
            low_mass3.extend([cat['nn3_dist'][k]])
            
#Exponential CDF
exp_sort3 = np.sort(exp_clouds3)
exp_cdf3 = np.arange(len(exp_clouds3))/ float(len(exp_clouds3))

#2D Coverage CDF
cov_sort3 = np.sort(cov_clouds3)
cov_cdf3 = np.arange(len(cov_clouds3))/float(len(cov_clouds3))

#Real Data CDF
real_sort3 = np.sort(real_clouds3)
real_cdf3 = np.arange(len(real_clouds3))/float(len(real_clouds3))

#High Mass CDF
hm_sort3 = np.sort(high_mass3)
hm_cdf3 = np.arange(len(high_mass3))/float(len(high_mass3))

#Low Mass CDF
lm_sort3 = np.sort(low_mass3)
lm_cdf3 = np.arange(len(low_mass3))/float(len(low_mass3))


## GRIDDING ##
xvals = np.linspace(0,2000,501)
cov_gridded3 = np.interp(xvals, cov_sort3, cov_cdf3)
real_gridded3 = np.interp(xvals, real_sort3, real_cdf3)
exp_gridded3 = np.interp(xvals, exp_sort3, exp_cdf3)
hm_gridded3 = np.interp(xvals, hm_sort3, hm_cdf3)
lm_gridded3 = np.interp(xvals, lm_sort3, lm_cdf3)

# Free up memory #
del exp_clouds3, exp_sort3, exp_cdf3
del real_clouds3, real_sort3, real_cdf3
del cov_clouds3, cov_sort3, cov_cdf3
del high_mass3, hm_sort3, hm_cdf3
del low_mass3, lm_sort3, lm_cdf3



### 5th Nearest Neighbor

real_clouds5 = []
exp_clouds5 = []
cov_clouds5 = []
for i in range(len(sources)):
    name = sources[i].split('_')[0]
    clouds = home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
    rand_cov = home_dir + 'gmc_spacing/data/sim_results/flat/'+sources[i]+'_flat_150pc.csv'
    rand_exp = home_dir + 'gmc_spacing/data/sim_results/exp/'+sources[i]+'_exp_150pc.csv'
    if os.path.isfile(clouds) and os.path.isfile(rand_cov) and os.path.isfile(rand_exp):
        real_clouds5.extend(pd.read_csv(clouds)['nn5_dist'])
        cov_clouds5.extend(pd.read_csv(rand_cov)['nn5_dist'])
        exp_clouds5.extend(pd.read_csv(rand_exp)['nn5_dist'])
        
        
## HIGH MASS LOW MASS ##
masses5 = []
high_mass5 = []
low_mass5 = []
for i in range(len(sources)):
    cat = pd.read_csv(home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv')
    masses5.extend(cat['MLUM_MSUN'])
mass_div5 = np.nanmedian(masses5)
for i in range(len(sources)):
    cat = pd.read_csv(home_dir + 'gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv')
    for k in range(len(cat)):
        if cat['MLUM_MSUN'][k] > mass_div5:
            high_mass5.extend([cat['nn5_dist'][k]])
        else:
            low_mass5.extend([cat['nn5_dist'][k]])
            
#Exponential CDF
exp_sort5 = np.sort(exp_clouds5)
exp_cdf5 = np.arange(len(exp_clouds5))/ float(len(exp_clouds5))

#2D Coverage CDF
cov_sort5 = np.sort(cov_clouds5)
cov_cdf5 = np.arange(len(cov_clouds5))/float(len(cov_clouds5))

#Real Data CDF
real_sort5 = np.sort(real_clouds5)
real_cdf5 = np.arange(len(real_clouds5))/float(len(real_clouds5))

#High Mass CDF
hm_sort5 = np.sort(high_mass5)
hm_cdf5 = np.arange(len(high_mass5))/float(len(high_mass5))

#Low Mass CDF
lm_sort5 = np.sort(low_mass5)
lm_cdf5 = np.arange(len(low_mass5))/float(len(low_mass5))


## GRIDDING ##
xvals = np.linspace(0,2000,501)
real_gridded5 = np.interp(xvals, real_sort5, real_cdf5)
cov_gridded5 = np.interp(xvals, cov_sort5, cov_cdf5)
exp_gridded5 = np.interp(xvals, exp_sort5, exp_cdf5)
hm_gridded5 = np.interp(xvals, hm_sort5, hm_cdf5)
lm_gridded5 = np.interp(xvals, lm_sort5, lm_cdf5)

# Free up memory #
del exp_clouds5, exp_sort5, exp_cdf5
del real_clouds5, real_sort5, real_cdf5
del cov_clouds5, cov_sort5, cov_cdf5
del high_mass5, hm_sort5, hm_cdf5
del low_mass5, lm_sort5, lm_cdf5


## MAKE PLOT ## 
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.size': 12})
fig, ax = plt.subplots(3,1, figsize=(5,15), dpi=155)



ax[0].plot(xvals, real_gridded, c='black', lw=2.5, label='Real GMCs - ALL')
ax[0].plot(xvals, hm_gridded, c='grey', ls='dashdot', alpha=0.6, lw=2, label='Real - High Mass')
ax[0].plot(xvals, lm_gridded, c='grey', ls='dotted', alpha=0.6, lw=2, label='Real - Low Mass')
ax[0].plot(xvals, exp_gridded, c='red', ls='dashed', lw=3, label='Exponential')
ax[0].plot(xvals, cov_gridded, c='navy', ls='dashed', lw=3, label='Flat')
ax[0].legend()
ax[0].tick_params(direction='in')
ax[0].grid(alpha=0.5)
ax[0].set_xlim(0,2000)
ax[0].set_ylabel(r'$f_{GMCs}$')
#ax[0].set_xlabel('Distance (pc)')
ax[0].set_title('First Nearest Neighbor')



ax[1].plot(xvals, real_gridded3, c='black', lw=2.5, label='Real GMCs - ALL')
ax[1].plot(xvals, hm_gridded3, c='grey', ls='dashdot', alpha=0.6, lw=2, label='Real GMCs - High Mass')
ax[1].plot(xvals, lm_gridded3, c='grey', ls='dotted', alpha=0.6, lw=2, label='Real GMCs - Low Mass')
ax[1].plot(xvals, exp_gridded3, c='red', ls='dashed', lw=3, label='Random - Exponential')
ax[1].plot(xvals, cov_gridded3, c='navy', ls='dashed', lw=3, label='Random - Flat')
ax[1].set_ylabel(r'$f_{GMCs}$')
#ax[1].legend()
ax[1].tick_params(direction='in')
ax[1].grid(alpha=0.5)
ax[1].set_xlim(0,2000)
#ax[1].set_xlabel('Distance (pc)')
ax[1].set_title('Third Nearest Neighbor')


ax[2].plot(xvals, real_gridded5, c='black', lw=2.5, label='Real GMCs - ALL')
ax[2].plot(xvals, hm_gridded5, c='grey', ls='dashdot', alpha=0.6, lw=2, label='Real GMCs - High Mass')
ax[2].plot(xvals, lm_gridded5, c='grey', ls='dotted', alpha=0.6, lw=2, label='Real GMCs - Low Mass')
ax[2].plot(xvals, exp_gridded5, c='red', ls='dashed', lw=3, label='Random - Exponential')
ax[2].plot(xvals, cov_gridded5, c='navy', ls='dashed', lw=3, label='Random - Flat')
#ax[2].legend()
ax[2].set_ylabel(r'$f_{GMCs}$')
ax[2].tick_params(direction='in')
ax[2].grid(alpha=0.5)
ax[2].set_xlim(0,2000)
ax[2].set_xlabel('Distance to Nth Nearest Neighbor (pc)')
ax[2].set_title('Fifth Nearest Neighbor')
plt.savefig(home_dir + 'gmc_spacing/figures/3panel_cdf_mass_split.pdf', bbox_inches='tight')