###################################################################################################
###################################################################################################
### THIS SCRIPT TAKES THE CPROPS .FITS CLOUD CATALOGS AND CALCULATES THE RELEVANT SPACING STATS ###
## THIS SCRIPT RETURNS A .CSV FILE WITH CLOUD LOCATION AND PROPERTIES, AND DISTANCE TO NEIGHBORS ##
###################################################################################################
###################################################################################################

######################################## IMPORTANT NOTE ###########################################
######## This script attempts to calculate distances in cartesian space & polar space. ############
###### As of writing (6/1/22), POLAR DISTANCES ARE BROKEN. Cartesian distances are reliable. ######
###################################################################################################
###################################################################################################


#### IMPORTS ####

import numpy as np
from astropy.table import Table
import pandas as pd
import os.path
import sys
sys.path.append('/home/machado.35/projects/gmc_spacing/modules') ## SET THIS WHEREVER YOUR MODULES ARE STORED
from deprojectGalaxy import deproject
import spacing_calc
from os import listdir


### SETUP ###
sources = listdir('/home/machado.35/projects/gmc_spacing/data/cloud_cats_fits/homogen_150pc_low_noise/')
sources_full = listdir('/home/machado.35/projects/gmc_spacing/data/cloud_cats_fits/homogen_150pc_low_noise/')
for i in range(len(sources)):
    sources[i] = sources[i].split('_')[0]
tab_loc = '/home/machado.35/projects/gmc_spacing/data/phangs_sample_table_v1p6.fits'
t = Table.read(tab_loc)

## DISTANCE CALCULATIONS ##
fp = '/home/machado.35/projects/gmc_spacing/data/cloud_cats_fits/homogen_150pc_low_noise/'
for i in range(len(sources)):
    cat_fp = fp+sources_full[i]
    cat = Table.read(cat_fp)
    if len(cat) > 5:
        row = int(np.where(t['name'] == sources[i])[0])
        galRA=t[row]['orient_ra']
        galDEC=t[row]['orient_dec']
        pa=t[row]['orient_posang']
        incl=t[row]['orient_incl']
        galDist=t[row]['dist'] #Mpc
        center_coord=[galRA, galDEC]
        x = np.array(cat['XCTR_DEG'])
        y = np.array(cat['YCTR_DEG'])
        ### DEPROJECT  - necessary for all calculations ###
        deproj = deproject(center_coord=center_coord, incl=incl, pa=pa, ra=x, dec=y, return_offset=True) 

        

        ### RADIAL APERTURE COUNTING ###
        radial_density_300 = spacing_calc.rad_density(x=deproj[2], y=deproj[3], rad=300, galdist=galDist*1e6)
        radial_density_500 = spacing_calc.rad_density(x=deproj[2], y=deproj[3], rad=500, galdist=galDist*1e6)
        radial_density_700 = spacing_calc.rad_density(x=deproj[2], y=deproj[3], rad=700, galdist=galDist*1e6)


        cloud = [k for k in range(len(cat))]
        df = pd.DataFrame()
        df['cloud_num'] = cat['CLOUDNUM']
        df['RA'] = x
        df['DEC'] = y
        df['MLUM_MSUN'] = cat['MLUM_MSUN']


        ### NEAREST NEIGHBOR ###
        #ALL DISTANCE CALCULATIONS #
        dist, dr, daz, index = spacing_calc.full_dist_calc(deproj[2],deproj[3], neighbor_num=1)
        pairs = np.array(list(map(list, zip(cloud, index))))
        # add distance (cartesian)
        df['nn_dist'] = np.deg2rad(dist[pairs[:,0],pairs[:,1]]) * galDist * 1e6
        # add distances (polar)
        df['nn_dr'] = np.deg2rad(dr[pairs[:,0],pairs[:,1]]) * galDist * 1e6
        df['nn_daz'] = np.deg2rad(daz[pairs[:,0],pairs[:,1]]) * galDist * 1e6
        # add neighbor index
        df['nn_index'] = index

        ### 3RD NEAREST NEIGHBOR ###
        #ALL DISTANCE CALCULATIONS #
        dist, dr, daz, index = spacing_calc.full_dist_calc(deproj[2],deproj[3], neighbor_num=3)
        pairs = np.array(list(map(list, zip(cloud, index))))
        # add distance (cartesian)
        df['nn3_dist'] = np.deg2rad(dist[pairs[:,0],pairs[:,1]]) * galDist * 1e6
        # add distances (polar)
        df['nn3_dr'] = np.deg2rad(dr[pairs[:,0],pairs[:,1]]) * galDist * 1e6
        df['nn3_daz'] = np.deg2rad(daz[pairs[:,0],pairs[:,1]]) * galDist * 1e6
        # add neighbor index
        df['nn3_index'] = index

        ### 5TH NEAREST NEIGHBOR ###
        #ALL DISTANCE CALCULATIONS #
        dist, dr, daz, index = spacing_calc.full_dist_calc(deproj[2],deproj[3], neighbor_num=5)
        pairs = np.array(list(map(list, zip(cloud, index))))
        # add distance (cartesian)
        df['nn5_dist'] = np.deg2rad(dist[pairs[:,0],pairs[:,1]]) * galDist * 1e6
        # add distances (polar)
        df['nn5_dr'] = np.deg2rad(dr[pairs[:,0],pairs[:,1]]) * galDist * 1e6
        df['nn5_daz'] = np.deg2rad(daz[pairs[:,0],pairs[:,1]]) * galDist * 1e6
        # add neighbor index
        df['nn5_index'] = index

        ### ADD RADIAL APERTURE CLOUD COUNTS ###
        df['clouds_r300pc'] = radial_density_300
        df['clouds_r500pc'] = radial_density_500
        df['clouds_r700pc'] = radial_density_700


        out_fn = '/home/machado.35/projects/gmc_spacing/data/cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
        df.to_csv(out_fn, index=False)