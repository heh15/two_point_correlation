import sys
sys.path.append('../')
sys.path.append('../utilsForTwopt/')

from scipy.interpolate import interp1d

from key_dirs import *
from utils_spacing import *

from reproject import reproject_interp

weight_map_dir = root_dir+'weight_maps_arm/'
sim_dir = root_dir+'sim_results_arm/'

dry_run = False

niters = 1e4
verbose = True

do_envrm_mask = True
do_strict_mask = False

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Run random placements
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# temporarily remove ngc0253 since it does not have environmental mask
if do_envrm_mask:
    source_list = [source for source in source_list if source not in ['ngc0253','ngc0300', 'ngc3489', 'ngc4459', 'ngc4476', 'ngc4477', 'ngc4596', 'ngc5236', 'ngc7793']] 


for this_source in source_list:

    # Extract details for this source
    this_row = deets_dict[this_source]
    this_ra = this_row['orient_ra']
    this_dec = this_row['orient_dec']
    this_pa = this_row['orient_posang']
    this_incl = this_row['orient_incl']
    this_dist = this_row['dist']
    this_lstar_deg = this_row['size_scalelength']/3600.

    # Extract the array and coverage file
    this_array = array_dict[this_source]
    this_cov_file = fid_cov_dict[this_source]
    this_fid_broad_map = broad_map_dict[this_source]
    this_fid_strict_mask = strict_mask_dict[this_source]

    # Read the radial profile
    radial_profile_table = Table.read('../radial_profile/'
            +this_source.upper()+'_base_annulus_0p5kpc.ecsv')
    CO_radial = radial_profile_table['I_CO21'] # / radial_profile_table['I_CO21'][0]
    radial_dist_pc = radial_profile_table['r_gal']*1000
    radial_profile_func = interp1d(radial_dist_pc, CO_radial, fill_value='extrapolate')
    
    # Read the kpc-scale moment 0 file
    hdulist = fits.open(this_fid_broad_map)
    kpc_mom0_map = hdulist[0].data
    hdulist.close()

    # Read the file
    hdulist = fits.open(this_cov_file)
    cov_map = hdulist[0].data

    # Make coordinates
    hdr = hdulist[0].header
    wcs = WCS(hdr, naxis=2)
    naxis = wcs._naxis

    grid  = np.indices((naxis[1],naxis[0]))
    ra_map, dec_map  = wcs.wcs_pix2world(grid[1],grid[0],0)

    center_coords = [this_ra, this_dec]

    # Deproject ra and dec to dx and dy and cylindrical coords
    radius_map, theta_map, dx_map, dy_map = deproject(
        center_coord=center_coords, incl=this_incl, pa=this_pa,
        ra=ra_map, dec=dec_map, return_offset=True)

    # import the environment mask
    if do_envrm_mask:
        hdulist = fits.open('../envrm_mask/masks_v5_simple/'+this_source.upper()+'_simple.fits')
        # reproject environment mask to the moment maps
        envrm_mask, _ = reproject_interp(hdulist, wcs, shape_out=np.shape(cov_map))
        envrm_mask = envrm_mask.astype(int)

        with fits.open(this_fid_strict_mask) as hdul:
            strict_mask = hdul[0].data
        strict_mask = np.max(strict_mask, axis=0)

    # Make the models
    mask = (cov_map > 0.25)
    if do_envrm_mask:
        mask = mask & (envrm_mask==6)
    if do_strict_mask:
        mask = mask & strict_mask 
    
    if len(np.where(mask>0)[0]) == 0:
        continue

    flat_model = mask*1.0
    flat_model = flat_model/np.nansum(flat_model)

    exp_model = mask*(np.exp(-1.0*radius_map/this_lstar_deg))
    exp_model = exp_model/np.nansum(exp_model)

    # add the model sample the radial profile
    radius_map_pc = radius_map * 3600 * 4.85 * this_dist
    rad_model = mask * radial_profile_func(radius_map_pc)
    rad_model = rad_model/np.nansum(rad_model)

    # change by Hao He (July 1st, 2024), add the model to sample the GMCs based on the kpc-scale CO emission
    Ico_model = mask*kpc_mom0_map 
    Ico_model = Ico_model/np.nansum(Ico_model)

    # Write to disk
    flat_model_file = weight_map_dir+this_source+'_flat_model.fits'
    exp_model_file = weight_map_dir+this_source+'_exp_model.fits'
    rad_model_file = weight_map_dir+this_source+'_rad_model.fits'
    Ico_model_file = weight_map_dir+this_source+'_Ico_model.fits'

    flat_hdu = fits.PrimaryHDU(data=flat_model,header=hdr)
    flat_hdu.writeto(flat_model_file, overwrite=True)

    exp_hdu = fits.PrimaryHDU(data=exp_model,header=hdr)
    exp_hdu.writeto(exp_model_file, overwrite=True)

    rad_hdu = fits.PrimaryHDU(data=rad_model,header=hdr)
    rad_hdu.writeto(rad_model_file, overwrite=True)

    Ico_hdu = fits.PrimaryHDU(data=Ico_model, header=hdr)
    Ico_hdu.writeto(Ico_model_file, overwrite=True)

    # Read the fiducial catalog
    this_cat = fid_cat_dict[this_source]

    this_nclouds = fid_ncloud_dict[this_source]
    n_clouds_to_place = this_nclouds*niters

    this_sim_dir = sim_dir+'flat/'
    flat_cat_outfile = this_sim_dir+this_source+'_'+senstvt+'_'+reslt+'_flat.fits'

    this_sim_dir = sim_dir+'exp/'
    exp_cat_outfile = this_sim_dir+this_source+'_'+senstvt+'_'+reslt+'_exp.fits'

    this_sim_dir = sim_dir+'rad/'
    rad_cat_outfile = this_sim_dir+this_source+'_'+senstvt+'_'+reslt+'_rad.fits'

    this_sim_dir = sim_dir+'Ico/'
    Ico_cat_outfile = this_sim_dir+this_source+'_'+senstvt+'_'+reslt+'_Ico.fits'

    if verbose:
        print("Placing clouds:")
        print("... galaxy:", this_source)
        print("... number of clouds:", n_clouds_to_place)
                  
    if dry_run == False:
        place_nclouds_in_map(
            exp_model_file, n_clouds_to_place, exp_cat_outfile,
            this_ra, this_dec, this_pa, this_incl)

        place_nclouds_in_map(
            flat_model_file, n_clouds_to_place, flat_cat_outfile,
            this_ra, this_dec, this_pa, this_incl)
        try:      
            place_nclouds_in_map(
                rad_model_file, n_clouds_to_place, rad_cat_outfile,
                this_ra, this_dec, this_pa, this_incl)
        except:
            print(this_source)
        place_nclouds_in_map(
            Ico_model_file, n_clouds_to_place, Ico_cat_outfile,
            this_ra, this_dec, this_pa, this_incl)
