import os, sys
sys.path.append('../')
sys.path.append('../aklpyutils/')

from utils_spacing import *
from key_dirs import *

### CALCULATE SPACING STATISTICS FOR SIMULATED CATALOGS ###
res = ['150']
sources = source_list
sources = ['ic1954']

#FLAT
radii = [300, 500, 700]
for k in range(len(res)):
    for i in range(len(sources)):
        image = derived_root_dir+sources[i]+'/'+sources[i]+'_12m+7m+tp_co21_'+str(res[k])+'pc_coverage2d.fits'
        cat_fp = root_dir+'cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
        sim_fp = root_dir+'sim_results/flat/'+sources[i]+'_flat_'+str(res[k])+'pc.csv'
        if os.path.isfile(image):
            pass
        else:
            image = derived_root_dir+sources[i]+'/'+sources[i]+'_7m+tp_co21_'+str(res[k])+'pc_coverage2d.fits' # ONLY FOR NGC0253
            if not os.path.isfile(image):
                image = derived_root_dir+sources[i]+'/'+sources[i]+'_7m_co21_'+str(res[k])+'pc_coverage2d.fits' # ONLY FOR NGC0300
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
        image = derived_root_dir+sources[i]+'/'+sources[i]+'_12m+7m+tp_co21_'+str(res[k])+'pc_coverage2d.fits'
        cat_fp = root_dir+'cloud_cats_csv/'+sources[i]+'_v4cprops_distcalc_46mK_150pc.csv'
        sim_fp = root_dir+'sim_results/exp/'+sources[i]+'_exp_'+str(res[k])+'pc.csv'

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
            
            
