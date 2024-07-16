import sys
sys.path.append('../')
# sys.path.append('../aklpyutils/')
sys.path.append('../utilsForTwopt/')
from key_dirs import *
from utilsForTwopt.utils_spacing import *

dry_run = False
verbose = True

bins_lo = np.arange(20)*25.-12.5
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

outfile_dir = '../cross_60pc_hi_results/'

# assumes the astrometry of the mocks and the image are the same -
# could fix this if any issues arise

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Place gaussians at each catalog peak
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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
#     model_file = mod_map_dir+this_source+'_source_model_90pc.fits'    

    # Read the file
    
    hdulist = fits.open(this_mom0_file)
    mom0_map = hdulist[0].data

    model_hdulist = fits.open(model_file)
    model_map = model_hdulist[0].data

    # Make coordinates
    
    hdr = hdulist[0].header
    wcs = WCS(hdr, naxis=2)
    naxis = wcs._naxis
    grid  = np.indices((naxis[1],naxis[0]))
    ra_map, dec_map  = wcs.wcs_pix2world(grid[1],grid[0],0)

    # Peak positions
    
    x_map, y_map  = wcs.wcs_world2pix(peak_ra,peak_dec,0)
    peak_val = mom0_map[np.round(y_map).astype(int),np.round(x_map).astype(int)]

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

    real_tab = Table()
    real_tab['bins_xlo'] = bins_lo
    real_tab['bins_xmid'] = bins_mid
    real_tab['bins_xhi'] = bins_hi
    real_tab['ymean']  = np.nanmean(all_real_crosses/all_real_counts, axis=0)
    real_tab['ymedian']  = np.nanmedian(all_real_crosses/all_real_counts, axis=0)    
    
    real_tab.write(outfile_name, format='fits', overwrite=True)
    
    all_model_crosses, all_model_counts = \
        calc_many_peak_image_crosses(
            peak_val*0+1, peak_deltax, peak_deltay,
            model_map, map_deltax, map_deltay,
            bins_lo, bins_hi)
    # print(all_model_crosses)
    
    outfile_name = outfile_dir + this_source + '_modelcross.fits'

    model_tab = Table()
    model_tab['bins_xlo'] = bins_lo
    model_tab['bins_xmid'] = bins_mid
    model_tab['bins_xhi'] = bins_hi
    model_tab['ymean']  = np.nanmean(all_model_crosses/all_model_counts, axis=0)
    model_tab['ymedian']  = np.nanmedian(all_model_crosses/all_model_counts, axis=0)    

    model_tab.write(outfile_name, format='fits', overwrite=True)    
    
    # Now run the 
    #model_file = mod_map_dir+this_source+'_source_model_150pc.fits'
    #model_hdu = fits.open(data=mock_map,header=mock_hdr)
    #model_hdu.writeto(model_file, overwrite=True)
