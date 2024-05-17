import sys
sys.path.append('../')
sys.path.append('../aklpyutils/')

from key_dirs import *
from utils_spacing import *

dry_run = False

niters = 1e4
verbose = True

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Run random placements
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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

    # Make the models
    flat_model = (cov_map > 0.25)*1.0
    flat_model = flat_model/np.nansum(flat_model)

    exp_model = (cov_map > 0.25)*(np.exp(-1.0*radius_map/this_lstar_deg))
    exp_model = exp_model/np.nansum(exp_model)

    # Write to disk
    flat_model_file = mod_map_dir+this_source+'_flat_model_150pc.fits'
    exp_model_file = mod_map_dir+this_source+'_exp_model_150pc.fits'

    flat_hdu = fits.PrimaryHDU(data=flat_model,header=hdr)
    flat_hdu.writeto(flat_model_file, overwrite=True)

    exp_hdu = fits.PrimaryHDU(data=exp_model,header=hdr)
    exp_hdu.writeto(exp_model_file, overwrite=True)

    # Read the fiducial catalog
    this_cat = fid_cat_dict[this_source]

    this_nclouds = fid_ncloud_dict[this_source]
    n_clouds_to_place = this_nclouds*niters

    this_sim_dir = sim_dir+'flat/'
    flat_cat_outfile = this_sim_dir+this_source+'_flat_150pc.fits'

    this_sim_dir = sim_dir+'exp/'
    exp_cat_outfile = this_sim_dir+this_source+'_exp_150pc.fits'

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
        
