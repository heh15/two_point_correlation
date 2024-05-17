######################################################
#### This script calculates Jeans & Toomre Length ####
####### and generates relevant plot for paper  #######


home_dir = '/home/machado.35/projects/'

from astropy.table import Table
from os import listdir
import os
import numpy as np
import sys
sys.path.append(home_dir + 'gmc_spacing/modules/')
from toomre_length import kappa2, lambda_T
from jeans_length import jeans_len
import astropy.units as u
import astropy.constants as c
import pandas as pd
from matplotlib.lines import Line2D
import scipy.stats
import matplotlib.pyplot as plt


def clouds_in_annulus(r_inner, r_outer, cloud_cat):
    
    inds = np.where(np.logical_and(cloud_cat['RGAL_KPC']>=r_inner, cloud_cat['RGAL_KPC']<=r_outer))
    
    return(inds)


source_list = listdir(home_dir + 'gmc_spacing/data/cloud_cats_csv/')
sources = {}
for i in range(len(source_list)):
    source_list[i] = source_list[i].split('_')[0]
    

megatables = listdir(home_dir + 'gmc_spacing/data/megatables/v3p0/annulus/')
for i in range(len(source_list)):
    sources[source_list[i]] = megatables[i]

galtable = Table.read(home_dir + 'gmc_spacing/data/phangs_sample_table_v1p6.fits')

fp = home_dir + 'gmc_spacing/data/megatables/v3p0/annulus/'
cloud_fp = home_dir + 'gmc_spacing/data/cloud_cats_fits/homogen_150pc_low_noise/'
df = pd.DataFrame(columns=list(sources.keys()))
df_jeans = pd.DataFrame(columns=list(sources.keys()))


## CALCULATE JEANS & TOOMRE LENGTH ##

for i in range(len(source_list)):
    name = source_list[i]
    megatable = Table.read(fp+sources[name])
    cloud_cat = cloud_fp+name+'_12m+7m+tp_co21_150pc_46mK_props.fits'
    if os.path.isfile(cloud_cat) == False:
        cloud_cat = cloud_fp+name+'_7m+tp_co21_150pc_46mK_props.fits'
    cloud_cat = Table.read(cloud_cat)
    toomre_length = np.zeros(65)
    jeans_length = np.zeros(65)
    r_inner = 0
    for k in range(len(megatable)):
        v_circ = (megatable['V_circ_CO21_URC'][k]*u.km/u.s).to(u.pc/u.s)
        beta = megatable['beta_CO21_URC'][k]
        r_gal = (megatable['r_gal'][k] * u.kpc).to(u.pc)
        k2 = kappa2(v_circ=v_circ, beta=beta, r_gal=r_gal)
        sigma_atom = megatable['Sigma_atom'][k]
        if np.isfinite(sigma_atom) == False:
            sigma_atom = 0
        sigma_mol = megatable['Sigma_mol'][k]
        sigma_gal = (sigma_atom + sigma_mol) * u.Msun / u.pc**2
        l_T = lambda_T(kappa2=k2, sigma_gal=sigma_gal)
        toomre_length[k] = l_T.value
        rho = sigma_gal / 100*u.pc
        r_outer = megatable['r_gal'][k]
        clouds = clouds_in_annulus(r_inner, r_outer, cloud_cat)
        vel = np.nanmean(cloud_cat['SIGV_KMS'])*u.km/u.s
        l_J = jeans_len(vel=vel, rho=rho)
        jeans_length[k] = l_J.value
        r_inner = r_outer
    df[name] = toomre_length
    df_jeans[name] = jeans_length

df.to_csv(home_dir + 'gmc_spacing/data/toomre_length_measurements.csv', index=False)
df_jeans.to_csv(home_dir + 'gmc_spacing/data/jeans_length_mean_vel.csv', index=False)


## MAKE PLOT ## 
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8,5), dpi=125)
med_spacings = []
toomres = []
jeans_5 = []
jeans_10 = []
for i in range(len(source_list)):
    name = source_list[i]
    megatable = Table.read(fp+sources[name])
    cloud_cat = cloud_fp+name+'_12m+7m+tp_co21_150pc_46mK_props.fits'
    cloud_cat_dist = pd.read_csv(home_dir + 'gmc_spacing/data/cloud_cats_csv/'+name+'_v4cprops_distcalc_46mK_150pc.csv')
    if os.path.isfile(cloud_cat) == False:
        cloud_cat = cloud_fp+name+'_7m+tp_co21_150pc_46mK_props.fits'
    cloud_cat = Table.read(cloud_cat)
    toomre_length = np.zeros(65)
    jeans_length = np.zeros(65)
    r_inner = 0
    for k in range(len(megatable)):
        v_circ = (megatable['V_circ_CO21_URC'][k]*u.km/u.s).to(u.pc/u.s)
        beta = megatable['beta_CO21_URC'][k]
        r_gal = (megatable['r_gal'][k] * u.kpc).to(u.pc)
        k2 = kappa2(v_circ=v_circ, beta=beta, r_gal=r_gal)
        sigma_atom = megatable['Sigma_atom'][k]
        if np.isfinite(sigma_atom) == False:
            sigma_atom = 6
        sigma_mol = megatable['Sigma_mol'][k]
        sigma_gal = (sigma_atom + sigma_mol) * u.Msun / u.pc**2
        l_T = lambda_T(kappa2=k2, sigma_gal=sigma_gal)
        toomre_length[k] = l_T.value
        rho = sigma_gal / 100*u.pc
        r_outer = megatable['r_gal'][k]
        clouds = clouds_in_annulus(r_inner, r_outer, cloud_cat)
        vel = 5*u.km/u.s
        l_J = jeans_len(vel=vel, rho=rho)
        jeans_length[k] = l_J.value
        
        cloud_dist = np.array(cloud_cat_dist['nn_dist'])[clouds]
        med_spacing = np.nanmedian(cloud_dist)
        
        med_spacings.append(med_spacing)
        toomres.append(l_T.value)
        jeans_5.append(l_J.value)
        
        ax[0].scatter(med_spacing, l_T.value, c='firebrick', edgecolor='maroon', alpha=0.5, s=18, zorder=5)
        ax[0].scatter(med_spacing, l_J.value, c='lightblue', edgecolor='teal', alpha=0.7, s=18, zorder=5)
        
        vel = 10*u.km/u.s
        l_J = jeans_len(vel=vel, rho=rho)
        jeans_10.append(l_J.value)
        jeans_length[k] = l_J.value
        
        ax[1].scatter(med_spacing, l_T.value, c='firebrick', edgecolor='maroon', alpha=0.5, s=18, zorder=5)
        ax[1].scatter(med_spacing, l_J.value, c='lightblue', edgecolor='teal', alpha=0.7, s=18, zorder=5)
        
        r_inner = r_outer
    
ax[0].set_xlabel('Median Cloud Spacing (pc)')
ax[0].set_ylabel('Scale Length (pc)')

ax[0].set_xscale('log')
ax[0].set_yscale('log')

ax[1].set_xlabel('Median Cloud Spacing (pc)')
ax[1].set_ylabel('Scale Length (pc)')

ax[1].set_xscale('log')
ax[1].set_yscale('log')

med_spacings = np.array(med_spacings)
toomres = np.array(toomres)
jeans_5 = np.array(jeans_5)
jeans_10 = np.array(jeans_10)

toomres = toomres[np.isfinite(med_spacings)]
jeans_5 = jeans_5[np.isfinite(med_spacings)]
jeans_10 = jeans_10[np.isfinite(med_spacings)]
med_spacings = med_spacings[np.isfinite(med_spacings)]

toomre_corr = scipy.stats.spearmanr(med_spacings[np.isfinite(toomres)], toomres[np.isfinite(toomres)])
jeans_corr_5 = scipy.stats.spearmanr(med_spacings[np.isfinite(jeans_5)], jeans_5[np.isfinite(jeans_5)])
jeans_corr_10 = scipy.stats.spearmanr(med_spacings[np.isfinite(jeans_10)], jeans_10[np.isfinite(jeans_10)])

toomre_legend = Line2D([0], [0], marker='o', color='w', label='Toomre Length', 
                       markerfacecolor='firebrick', markersize=5, markeredgecolor='maroon')

jeans_legend = Line2D([0], [0], marker='o', color='w', label='Jeans Length', 
                       markerfacecolor='lightblue', markersize=5, markeredgecolor='teal')

#r'$\rho$ = '+str(np.around(toomre_corr[0], 2))
legend = [toomre_legend, jeans_legend]


#ax[0].legend(handles=legend)
ax[0].tick_params(direction='in', which='both')

ax[1].legend(handles=legend, loc='lower right')
ax[1].tick_params(direction='in', which='both')

ax[0].plot([0.1,10000], [0.1,10000], ls='--', color='black', lw=3)
ax[1].plot([0.1,10000], [0.1,10000], ls='--', color='black', lw=3)

ax[0].set_xlim(20,6500)
ax[1].set_xlim(20,6500)

ax[0].set_ylim(1,2500)
ax[1].set_ylim(1,4500)
ax[0].set_title(r'$\sigma$ = 5 km$s^{-1}$')
ax[1].set_title(r'$\sigma$ = 10 km$s^{-1}$')
fig.tight_layout(pad=1.0)
plt.savefig(home_dir + 'gmc_spacing/figures/toomre_jean_length.pdf', bbox_inches='tight')