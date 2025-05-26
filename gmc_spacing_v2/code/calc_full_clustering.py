import sys
sys.path.append('../')
sys.path.append('../utilsForTwopt/') 

from key_dirs import *
from utils_spacing import *
from astropy.utils.console import ProgressBar

corr_dir = '../corr_results/'
sim_dir = root_dir+'sim_results/'

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Define bins and other tuning parameters
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

log_xmin = 1.0
log_xbin = 0.1
log_xlo = log_xmin + np.arange(31)*log_xbin
log_xhi = log_xlo + log_xbin

bins_lo = 10.**(log_xlo)
bins_hi = 10.**(log_xhi)
bins_mid = (bins_lo+bins_hi)*0.5
bins_geo = np.sqrt(bins_lo*bins_hi)
nbins = len(bins_lo)

thresh_cloud_count = 5
cat_iter_fac = 100
exp_vs_flat_fac = 10

do_exp_vs_flat = False
do_monte = False

do_envrm_corr = False
do_select_bright = False

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Loop over sources
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
source_list=['ngc3627']
for this_source in source_list:
    
    # Define output filenames
    tabfile_real_vs_flat = corr_dir + this_source + '_realvsflat.fits'
    tabfile_real_vs_exp = corr_dir + this_source + '_realvsexp.fits'
    tabfile_exp_vs_flat = corr_dir + this_source + '_expvsflat.fits'
    tabfile_real_vs_rad = corr_dir + this_source + '_realvsrad.fits'
    tabfile_real_vs_Ico = corr_dir + this_source + '_realvsIco.fits'
    tabfile_Ico_vs_exp = corr_dir + this_source + '_Icovsexp.fits'

    # Extract source information
    this_row = deets_dict[this_source] 
    this_ra = this_row['orient_ra']
    this_dec = this_row['orient_dec']
    this_pa = this_row['orient_posang']
    this_incl = this_row['orient_incl']
    this_dist = this_row['dist']
    this_center_coord = [this_ra, this_dec]

    # Deg to physics
    deg_to_pc = this_dist*1e6*np.pi/180.
    
    # Read fiducial catalog
    this_real_cat = Table.read(fid_cat_dict[this_source].replace('_props.fits','_props_valueadded.fits'))
    if do_envrm_corr == True:
        this_real_cat = this_real_cat[this_real_cat['ENV_ARM']==1] 
    if do_select_bright == True:
        this_real_cat = this_real_cat[this_real_cat['TMAX_K']>0.5]
    if len(this_real_cat) < thresh_cloud_count:
        continue
    x_real = this_real_cat['DELTAX_DEG']*deg_to_pc
    y_real = this_real_cat['DELTAY_DEG']*deg_to_pc
    n_clouds = len(this_real_cat)

    n_sim_to_use = n_clouds*cat_iter_fac
    n_for_exp_vs_flat = n_clouds*exp_vs_flat_fac
    
    # Read the random catalogs and recast variables

    # ... exponential catalog
    exp_sim_dir = sim_dir + 'exp/'
    this_exp_cat_file = exp_sim_dir+this_source+'_'+senstvt+'_'+reslt+'_exp.fits'
    this_exp_cat = Table.read(this_exp_cat_file)

    x_exp = this_exp_cat['dx']*deg_to_pc
    y_exp = this_exp_cat['dy']*deg_to_pc

    x_exp = x_exp[0:n_sim_to_use]
    y_exp = y_exp[0:n_sim_to_use]
    
    n_exp = len(x_exp)
    
    # ... flat catalog
    flat_sim_dir = sim_dir + 'flat/'
    this_flat_cat_file = flat_sim_dir+this_source+'_'+senstvt+'_'+reslt+'_flat.fits'
    this_flat_cat = Table.read(this_flat_cat_file)

    x_flat = this_flat_cat['dx']*deg_to_pc
    y_flat = this_flat_cat['dy']*deg_to_pc

    x_flat = x_flat[0:n_sim_to_use]
    y_flat = y_flat[0:n_sim_to_use]

    n_flat = len(x_flat)

    # ... radial catalog
    rad_sim_dir = sim_dir + 'rad/'
    this_rad_cat_file = rad_sim_dir+this_source+'_'+senstvt+'_'+reslt+'_rad.fits'
    this_rad_cat = Table.read(this_rad_cat_file)

    x_rad = this_rad_cat['dx']*deg_to_pc
    y_rad = this_rad_cat['dy']*deg_to_pc

    x_rad = x_rad[0:n_sim_to_use]
    y_rad = y_rad[0:n_sim_to_use]

    n_rad = len(x_rad)
    
    # ... Ico catalog
    Ico_sim_dir = sim_dir + 'Ico/'
    this_Ico_cat_file = Ico_sim_dir+this_source+'_'+senstvt+'_'+reslt+'_Ico.fits'
    this_Ico_cat = Table.read(this_Ico_cat_file)

    x_Ico = this_Ico_cat['dx']*deg_to_pc
    y_Ico = this_Ico_cat['dy']*deg_to_pc

    x_Ico = x_Ico[0:n_sim_to_use]
    y_Ico = y_Ico[0:n_sim_to_use]

    n_Ico = len(x_Ico)

    if verbose:
        print("")
        print("... galaxy: ", this_source)
        print("... real peaks: ", n_clouds)
        print("... flat peaks: ", n_flat)
        print("... exp peaks: ", n_exp)
        print("")
    
    # Calculate fiducial correlations

    # ... real vs. flat
    print("... calculating real vs. flat")
    twopt_real_flat, twopt_flat_err = \
        calc_twopoint(
            x_real, y_real,
            x_flat, y_flat,
            bins_lo, bins_hi)
    
    # ... real vs. exponential
    print("... calculating real vs. exp")
    twopt_real_exp, twopt_exp_err = \
        calc_twopoint(
            x_real, y_real,
            x_exp, y_exp,
            bins_lo, bins_hi)
    
    # ... exponential vs. flat
    if do_exp_vs_flat:
        print("... calculating exp vs. flat")
        twopt_exp_flat = \
            calc_twopoint(
                x_exp[0:n_for_exp_vs_flat], y_exp[0:n_for_exp_vs_flat],
                x_flat, y_flat,
                bins_lo, bins_hi)
    else:
        twopt_exp_flat = twopt_real_flat*np.nan

    # ... real vs radial
    print('... calculating real vs rad')
    twopt_real_rad, twopt_rad_err = \
        calc_twopoint(
            x_real, y_real,
            x_rad, y_rad,
            bins_lo, bins_hi)

    # ... real vs Ico weighted
    print("... calculating real vs Ico")
    twopt_real_Ico, twopt_Ico_err = \
        calc_twopoint(
            x_real, y_real,
            x_Ico, y_Ico, 
            bins_lo, bins_hi)

    # Loop over realizations of random catalogs
    monte_flat_twopt = np.zeros((cat_iter_fac,nbins))*np.nan
    monte_exp_twopt = np.zeros((cat_iter_fac,nbins))*np.nan
    monte_rad_twopt = np.zeros((cat_iter_fac,nbins))*np.nan
    monte_Ico_twopt = np.zeros((cat_iter_fac,nbins))*np.nan

    print("... monte carlo calculation")
    counter = 0
    
    with ProgressBar(cat_iter_fac) as bar:
        
        for ii in range(cat_iter_fac):

            if do_monte == False:
                continue

            # Generate false data - note that the split avoid overlap
            # at low distance by construction.
            
            exp_mask = np.zeros_like(x_exp)
            exp_mask[counter:counter+n_clouds] = 1

            flat_mask = np.zeros_like(x_flat)
            flat_mask[counter:counter+n_clouds] = 1
            
            monte_exp_x_data = x_exp[exp_mask==1]
            monte_exp_y_data = y_exp[exp_mask==1]

            monte_exp_x_rand = x_exp[exp_mask==0]
            monte_exp_y_rand = y_exp[exp_mask==0]
            
            monte_flat_x_data = x_flat[flat_mask==1]
            monte_flat_y_data = y_flat[flat_mask==1]

            monte_flat_x_rand = x_flat[flat_mask==0]
            monte_flat_y_rand = y_flat[flat_mask==0]
            
            counter += n_clouds

            this_twopt_monte_flat = \
                calc_twopoint(
                    monte_flat_x_data, monte_flat_y_data,
                    monte_flat_x_rand, monte_flat_y_rand,
                    bins_lo, bins_hi, verbose=False)
        
            this_twopt_monte_exp = \
                calc_twopoint(
                    monte_exp_x_data, monte_exp_y_data,
                    monte_exp_x_rand, monte_exp_y_rand,
                    bins_lo, bins_hi, verbose=False)
            
            monte_flat_twopt[ii,:] = this_twopt_monte_flat
            monte_exp_twopt[ii,:] = this_twopt_monte_exp

            bar.update()
        
    # Produce final data table

    print("... making table and writing to disk")
    
    monte_flat_05 = np.zeros(nbins)*np.nan
    monte_flat_16 = np.zeros(nbins)*np.nan
    monte_flat_50 = np.zeros(nbins)*np.nan
    monte_flat_84 = np.zeros(nbins)*np.nan
    monte_flat_95 = np.zeros(nbins)*np.nan

    monte_exp_05 = np.zeros(nbins)*np.nan
    monte_exp_16 = np.zeros(nbins)*np.nan
    monte_exp_50 = np.zeros(nbins)*np.nan
    monte_exp_84 = np.zeros(nbins)*np.nan
    monte_exp_95 = np.zeros(nbins)*np.nan
    
    for ii in range(nbins):

        if do_monte == False:
            continue
            
        this_monte_flat_bin = monte_flat_twopt[:,ii]
        this_monte_exp_bin = monte_exp_twopt[:,ii]
        
        monte_flat_05[ii] = np.percentile(this_monte_flat_bin,5)
        monte_flat_16[ii] = np.percentile(this_monte_flat_bin,16)
        monte_flat_50[ii] = np.percentile(this_monte_flat_bin,50)
        monte_flat_84[ii] = np.percentile(this_monte_flat_bin,84)
        monte_flat_95[ii] = np.percentile(this_monte_flat_bin,95)

        monte_exp_05[ii] = np.percentile(this_monte_exp_bin,5)
        monte_exp_16[ii] = np.percentile(this_monte_exp_bin,16)
        monte_exp_50[ii] = np.percentile(this_monte_exp_bin,50)
        monte_exp_84[ii] = np.percentile(this_monte_exp_bin,84)
        monte_exp_95[ii] = np.percentile(this_monte_exp_bin,95)
              
    # Write to disk
    twopt_real_vs_flat = {
        'xlo':bins_lo,
        'xmid':bins_mid,
        'xgeo':bins_geo,
        'xhi':bins_hi,
        'twopt':twopt_real_flat,
        'twopt_err':twopt_flat_err,
        'monte05':monte_flat_05,
        'monte16':monte_flat_16,
        'monte50':monte_flat_50,
        'monte84':monte_flat_84,
        'monte95':monte_flat_95,
    }
    tab_real_vs_flat = Table(twopt_real_vs_flat)
    tab_real_vs_flat.write(tabfile_real_vs_flat,format='fits',overwrite=True)
    
    twopt_real_vs_exp = {
        'xlo':bins_lo,
        'xmid':bins_mid,
        'xgeo':bins_geo,
        'xhi':bins_hi,
        'twopt':twopt_real_exp,
        'twopt_err':twopt_exp_err,
        'monte05':monte_exp_05,
        'monte16':monte_exp_16,
        'monte50':monte_exp_50,
        'monte84':monte_exp_84,
        'monte95':monte_exp_95,
    }
    tab_real_vs_exp = Table(twopt_real_vs_exp)
    tab_real_vs_exp.write(tabfile_real_vs_exp,format='fits',overwrite=True)
    
    twopt_exp_vs_flat = {
        'xlo':bins_lo,
        'xmid':bins_mid,
        'xgeo':bins_geo,
        'xhi':bins_hi,
        'twopt':twopt_exp_flat,
    }
    tab_exp_vs_flat = Table(twopt_exp_vs_flat)
    tab_exp_vs_flat.write(tabfile_exp_vs_flat,format='fits',overwrite=True)
    
    twopt_real_vs_rad = {
        'xlo':bins_lo,
        'xmid':bins_mid,
        'xgeo':bins_geo,
        'xhi':bins_hi,
        'twopt':twopt_real_rad,
        'twopt_err':twopt_rad_err,
    }
    tab_real_vs_rad = Table(twopt_real_vs_rad)
    tab_real_vs_rad.write(tabfile_real_vs_rad,format='fits',overwrite=True)

    twopt_real_vs_Ico = {
        'xlo':bins_lo,
        'xmid':bins_mid,
        'xgeo':bins_geo,
        'xhi':bins_hi,
        'twopt':twopt_real_Ico,
        'twopt_err':twopt_Ico_err,
    }
    tab_real_vs_Ico = Table(twopt_real_vs_Ico)
    tab_real_vs_Ico.write(tabfile_real_vs_Ico,format='fits',overwrite=True)
