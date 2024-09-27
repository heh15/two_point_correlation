import sys
sys.path.append('../')
# sys.path.append('../aklpyutils/')
sys.path.append('../utilsForTwopt/')
from key_dirs import *
from utilsForTwopt.utils_spacing import *

dry_run = False
verbose = True

bins_lo = np.arange(80)*25.-12.5
bins_hi = bins_lo+25.
bins_mid = (bins_lo+bins_hi)*0.5
nbins = len(bins_lo)

# log_xmin = 1.0
# log_xbin = 0.1
# log_xlo = log_xmin + np.arange(31)*log_xbin
# log_xhi = log_xlo + log_xbin
# 
# bins_lo = 10.**(log_xlo)
# bins_hi = 10.**(log_xhi)
# bins_mid = (bins_lo+bins_hi)*0.5
# bins_geo = np.sqrt(bins_lo*bins_hi)
# bins_mid = bins_geo
# nbins = len(bins_lo)

mod_map_dir = '../model_maps_hi_150pc_Tmax/'
outfile_dir = '../cross_results_Tmax/'
cloud_rad_dir = '../cloud_rad_profl_Tmax/'

# assumes the astrometry of the mocks and the image are the same -
# could fix this if any issues arise

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Place gaussians at each catalog peak
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# source_list = ['ngc0253']
for this_source in source_list:

    print("")
    print("Running peak image crosses for ", this_source)
    print("")
    
    # Extract details for this source
    
    this_row = deets_dict[this_source]
    this_ra = this_row['orient_ra']
    this_dec = this_row['orient_dec']
    this_pa = this_row['orient_posang']
    this_incl = this_row['orient_incl']
    this_dist = this_row['dist']
    this_lstar_deg = this_row['size_scalelength']/3600.

    # Get the catalog
    
    this_cat_file = fid_cat_dict[this_source]
    
    # Read the catalog
    
    this_cat = Table.read(this_cat_file)
    peak_ra = np.array(this_cat['XMAX_DEG'])
    peak_dec = np.array(this_cat['YMAX_DEG'])

    # Note the moment0 file and the model file
    
    this_mom0_file = mom0_map_dict[this_source]
    model_file = mod_map_dir+this_source+'_source_model_'+reslt+'.fits'
    modelbeam_file = mod_map_dir+this_source+'_source_model_beam_'+reslt+'.fits'
    modelbakgrd_file = mod_map_dir+this_source+'_source_model_bakgrd_'+reslt+'.fits'

    # note the kpc resolution map
    this_1kpc_mom0_file = broad_map_dict[this_source]

    # Read the file 
    hdulist = fits.open(this_mom0_file)
    mom0_map = hdulist[0].data

    model_hdulist = fits.open(model_file)
    model_map = model_hdulist[0].data

    modelbeam_hdulist = fits.open(modelbeam_file)
    modelbeam_map = modelbeam_hdulist[0].data

    modelbakgrd_hdulist = fits.open(modelbakgrd_file)
    modelbakgrd_map = modelbakgrd_hdulist[0].data

    this_1kpc_map_hdulist = fits.open(this_1kpc_mom0_file)
    this_1kpc_map = this_1kpc_map_hdulist[0].data

    # Make coordinates
    
    hdr = hdulist[0].header
    wcs = WCS(hdr, naxis=2)
    naxis = wcs._naxis
    grid  = np.indices((naxis[1],naxis[0]))
    ra_map, dec_map  = wcs.wcs_pix2world(grid[1],grid[0],0)

    # Peak positions
    
    x_map, y_map  = wcs.wcs_world2pix(peak_ra,peak_dec,0)
    peak_val = mom0_map[np.round(y_map).astype(int),np.round(x_map).astype(int)]
    peak_1kpc_val = this_1kpc_map[np.round(y_map).astype(int),np.round(x_map).astype(int)]

    # Deproject the map
    
    center_coords = [this_ra, this_dec]
    map_rgal, map_theta, map_deltax, map_deltay = \
        deproject(
            center_coord=center_coords, incl=this_incl, pa=this_pa,
            ra=ra_map, dec=dec_map, return_offset=True)
    map_deltax = map_deltax/180.*np.pi*this_dist*1e6
    map_deltay = map_deltay/180.*np.pi*this_dist*1e6        

    # Deproject the catalog
    
    center_coords = [this_ra, this_dec]
    peak_rgal, peak_theta, peak_deltax, peak_deltay = \
        deproject(
            center_coord=center_coords, incl=this_incl, pa=this_pa,
            ra=peak_ra, dec=peak_dec, return_offset=True)
    peak_deltax = peak_deltax/180.*np.pi*this_dist*1e6
    peak_deltay = peak_deltay/180.*np.pi*this_dist*1e6        
    
    # Calculate the crosses for the real data
    
    all_real_crosses, all_real_counts = \
        calc_many_peak_image_crosses(
            peak_val*0+1, peak_deltax, peak_deltay,
            mom0_map, map_deltax, map_deltay,
            bins_lo, bins_hi)

    #print(all_real_crosses.shape)
    
    outfile_name = outfile_dir + this_source + '_realcross.fits'
    nan_mask = (np.isnan(all_real_crosses)) | (np.isnan(all_real_counts))
    avrg_intst_maskd = np.ma.masked_array(all_real_crosses / all_real_counts,
                                        nan_mask)
    y_norm = avrg_intst_maskd / (avrg_intst_maskd[:,0])[:,np.newaxis]
    
    real_tab = Table()
    real_tab['bins_xlo'] = bins_lo
    real_tab['bins_xmid'] = bins_mid
    real_tab['bins_xhi'] = bins_hi
    real_tab['I_mean']  = np.ma.average(avrg_intst_maskd, axis=0)
    real_tab['I_mean_weightd'] = np.ma.average(avrg_intst_maskd, 
                                            weights= peak_val, axis=0)
    real_tab['I_median']  = np.nanmedian(all_real_crosses/all_real_counts, axis=0)    
    real_tab['y_mean'] = np.nanmean(y_norm, axis=0) 
    real_tab['y_median'] = np.nanmedian(y_norm, axis=0)
    real_tab['y_mean_weightd'] = np.ma.average(y_norm, weights=peak_val, axis=0)
    real_tab.write(outfile_name, format='fits', overwrite=True)    

    # calculate the crosses for the kpc-scale model map

    all_model_crosses_1kpc, all_model_counts_1kpc = \
        calc_many_peak_image_crosses(
            peak_val*0+1, peak_deltax, peak_deltay,
            this_1kpc_map, map_deltax, map_deltay,
            bins_lo, bins_hi)

    outfile_name = outfile_dir + this_source + '_modelcross_1kpc.fits'
    nan_mask = (np.isnan(all_model_crosses_1kpc)) | (np.isnan(all_model_counts_1kpc))
    avrg_intst_maskd = np.ma.masked_array(all_model_crosses_1kpc / all_model_counts_1kpc,
                                        nan_mask)
    y_norm = avrg_intst_maskd / (avrg_intst_maskd[:,0])[:,np.newaxis]

    model_tab_1kpc  = Table()
    model_tab_1kpc['bins_xlo'] = bins_lo
    model_tab_1kpc['bins_xmid'] = bins_mid
    model_tab_1kpc['bins_xhi'] = bins_hi
    model_tab_1kpc['I_mean']  = np.ma.average(avrg_intst_maskd, axis=0)
    model_tab_1kpc['I_mean_weightd'] = np.ma.average(avrg_intst_maskd,
                                            weights= peak_val, axis=0)
    model_tab_1kpc['I_median']  = np.nanmedian(all_model_crosses_1kpc/all_model_counts_1kpc, axis=0)
    model_tab_1kpc['y_mean'] = np.nanmean(y_norm, axis=0)
    model_tab_1kpc['y_median'] = np.nanmedian(y_norm, axis=0)
    model_tab_1kpc['y_mean_weightd'] = np.ma.average(y_norm, weights=peak_val, axis=0)
    model_tab_1kpc.write(outfile_name, format='fits', overwrite=True)

    # calculate the crosses for the model image

    all_model_crosses, all_model_counts = \
        calc_many_peak_image_crosses(
            peak_val*0+1, peak_deltax, peak_deltay,
            model_map, map_deltax, map_deltay,
            bins_lo, bins_hi)
    # print(all_model_crosses)

    outfile_name = outfile_dir + this_source + '_modelcross.fits'
    nan_mask = (np.isnan(all_model_crosses)) | (np.isnan(all_model_counts))
    avrg_intst_maskd = np.ma.masked_array(all_model_crosses / all_model_counts,
                                        nan_mask)
    y_norm = avrg_intst_maskd / (avrg_intst_maskd[:,0])[:,np.newaxis]

    model_tab = Table()
    model_tab['bins_xlo'] = bins_lo
    model_tab['bins_xmid'] = bins_mid
    model_tab['bins_xhi'] = bins_hi
    model_tab['I_mean']  = np.nanmean(all_model_crosses/all_model_counts, axis=0)
    model_tab['I_mean_weightd'] = np.ma.average(avrg_intst_maskd,
                                            weights= peak_val, axis=0)
    model_tab['I_median']  = np.nanmedian(all_model_crosses/all_model_counts, axis=0)    
    model_tab['y_mean'] = np.nanmean(y_norm, axis=0)
    model_tab['y_median'] = np.nanmedian(y_norm, axis=0)
    model_tab['y_mean_weightd'] = np.ma.average(y_norm, weights=peak_val, axis=0)
    model_tab.write(outfile_name, format='fits', overwrite=True)    

    # calculate the crosses for model image convolved with beam

    all_modelbeam_crosses, all_modelbeam_counts = \
        calc_many_peak_image_crosses(
            peak_val*0+1, peak_deltax, peak_deltay,
            modelbeam_map, map_deltax, map_deltay,
            bins_lo, bins_hi)
    # print(all_model_crosses)

    outfile_name = outfile_dir + this_source + '_modelbeam.fits'
    nan_mask = (np.isnan(all_modelbeam_crosses)) | (np.isnan(all_modelbeam_counts))
    avrg_intst_maskd = np.ma.masked_array(all_modelbeam_crosses / all_modelbeam_counts, nan_mask)
    y_norm = avrg_intst_maskd / (avrg_intst_maskd[:,0])[:,np.newaxis]

    modelbeam_tab = Table()
    modelbeam_tab['bins_xlo'] = bins_lo
    modelbeam_tab['bins_xmid'] = bins_mid
    modelbeam_tab['bins_xhi'] = bins_hi
    modelbeam_tab['I_mean']  = np.nanmean(all_modelbeam_crosses/all_modelbeam_counts, axis=0)
    modelbeam_tab['I_mean_weightd'] = np.ma.average(avrg_intst_maskd,
                                            weights= peak_val, axis=0)
    modelbeam_tab['I_median']  = np.nanmedian(all_modelbeam_crosses/all_modelbeam_counts, axis=0)
    modelbeam_tab['y_mean'] = np.nanmean(y_norm, axis=0)
    modelbeam_tab['y_median'] = np.nanmedian(y_norm, axis=0)
    modelbeam_tab['y_mean_weightd'] = np.ma.average(y_norm, weights=peak_val, axis=0)
    modelbeam_tab.write(outfile_name, format='fits', overwrite=True) 

    # calculate the crosses for model image with background emission 
    all_modelbakgrd_crosses, all_modelbakgrd_counts = \
        calc_many_peak_image_crosses(
            peak_val*0+1, peak_deltax, peak_deltay,
            modelbakgrd_map, map_deltax, map_deltay,
            bins_lo, bins_hi)

    # calculate the radial profile of individual clouds without including contribution from other GMCs
    flux_val = np.array(this_cat['FLUX_KKMS_PC2'])
    this_width_sig = np.array(np.sqrt(this_cat['MOMMAJ_NODC_PC']*this_cat['MOMMIN_NODC_PC']))
    peak_values = flux_val / (2*math.pi*this_width_sig**2)
    gaussian_rad = []
    for i, peak in enumerate(peak_val):
        gaussian_rad.append(peak * np.exp(-1./2.*bins_mid**2/this_width_sig[i]**2))
    gaussian_rad = np.array(gaussian_rad)

    # write the radial profile of individual clouds into a Table
    this_cat_rad = Table()
    this_cat_rad['name'] = this_cat['GALAXY']
    this_cat_rad['arm_lab'] = this_cat['ENV_ARM']
    this_cat_rad['ring_lab'] = this_cat['ENV_RING']
    this_cat_rad['bar_lab'] = this_cat['ENV_BAR']
    this_cat_rad['peak_val'] = peak_val
    this_cat_rad['peak_1kpc_val'] = peak_1kpc_val
    this_cat_rad['radius_1sigma_pc'] = this_width_sig 
    this_cat_rad['radial_bins'] = np.tile(bins_mid, (len(this_cat),1))
    this_cat_rad['radial_real_val'] = all_real_crosses / all_real_counts
    this_cat_rad['radial_model_val'] = all_model_crosses / all_model_counts
    this_cat_rad['radial_modelbeam_val'] = all_modelbeam_crosses / all_modelbeam_counts
    this_cat_rad['radial_modelbakgrd_val'] = all_modelbakgrd_crosses / all_modelbakgrd_counts
    this_cat_rad['radial_singlgmc_val'] = gaussian_rad
    this_cat_rad['radial_val_model_1kpc'] = all_model_crosses_1kpc / all_model_counts_1kpc
    # save it to a new table
    outfile_name = cloud_rad_dir+this_source+'_'+reslt+'_'+senstvt+'_cat.fits'
    this_cat_rad.write(outfile_name, format='fits', overwrite=True)

    # Now run the 
    #model_file = mod_map_dir+this_source+'_source_model_150pc.fits'
    #model_hdu = fits.open(data=mock_map,header=mock_hdr)
    #model_hdu.writeto(model_file, overwrite=True)
