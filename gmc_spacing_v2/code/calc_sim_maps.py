import sys
sys.path.append('../')
# sys.path.append('../aklpyutils/')
sys.path.append('../utilsForTwopt/')

import matplotlib.pyplot as plt
import pickle

from key_dirs import *
# from aklpyutils.utils_spacing import *
from utilsForTwopt.utils_spacing import *

mod_map_dir = root_dir+'model_maps_'+senstvt+'_'+reslt+'/'
mod_map_dir = root_dir+'model_maps_hi_150pc_Tmax1/'

dry_run = False
verbose = True

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Place gaussians at each catalog peak
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

flux_fracs = []
# source_list = ['ngc0253']
for this_source in source_list:

    print("")
    print("Making a mock source map for ", this_source)
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
    flux_val = np.array(this_cat['FLUX_KKMS_PC2'])
    Tmax = np.array(this_cat['TMAX_K'])
    pa = np.array(this_cat['POSANG']) * 180 / math.pi # convert to deg
    # Extract the moment0 file
    this_mom0_file = mom0_map_dict[this_source]
    # Extract the 1kpc map file
    this_1kpc_mom0_file = broad_map_dict[this_source]
    # Extract the 1kpc strict mask file
    this_1kpc_strict_mask_file = strict_1kpc_mask_dict[this_source] 

    # Read the file
    hdulist = fits.open(this_mom0_file)
    mom0_map = hdulist[0].data

    # Read the kpc file
    this_1kpc_map_hdulist = fits.open(this_1kpc_mom0_file)
    this_1kpc_map = this_1kpc_map_hdulist[0].data

    # Read the kpc strict mask
    with fits.open(this_1kpc_strict_mask_file) as hdul:
        strict_mask_1kpc = hdul[0].data
        strict_mask_1kpc = np.max(strict_mask_1kpc, axis=0)

    # Make coordinates
    hdr = hdulist[0].header
    wcs = WCS(hdr, naxis=2)
    naxis = wcs._naxis
    grid  = np.indices((naxis[1],naxis[0]))
    ra_map, dec_map  = wcs.wcs_pix2world(grid[1],grid[0],0)

    # Peak positions
    x_map, y_map  = wcs.wcs_world2pix(peak_ra,peak_dec,0)
    peak_val = mom0_map[np.round(y_map).astype(int),np.round(x_map).astype(int)]   
    # Make the model map
    beam_fwhm = hdr['BMAJ']
    peak_fwhm_pc = np.array(2.354*np.sqrt(this_cat['MOMMAJ_NODC_PC']*this_cat['MOMMIN_NODC_PC']))
    peak_fwhm = peak_fwhm_pc / 4.85 / this_dist / 3600
    mock_map = make_mock_image(
        Tmax, peak_ra, peak_dec, beam_fwhm,
        ra_map, dec_map, this_dist, peak_fwhm=peak_fwhm, mask_to_blank=(np.isfinite(mom0_map)==False))
    peak_temp = mock_map[np.round(y_map).astype(int),np.round(x_map).astype(int)]
    ratio = np.nanmedian(peak_val) / np.nanmedian(peak_temp)
    mock_map = mock_map * ratio

    # Make the model map with elliptical beam
    peak_fwhm = 2.354 * this_cat['MOMMAJ_NODC_PC'] / 4.85 / this_dist / 3600
    peak_fwhm_min = 2.354 * this_cat['MOMMIN_NODC_PC'] / 4.85 / this_dist / 3600
    mock_map_elptc = make_mock_image(
        Tmax, peak_ra, peak_dec, beam_fwhm,
        ra_map, dec_map, this_dist, peak_fwhm=peak_fwhm, mask_to_blank=(np.isfinite(mom0_map)==False), peak_fwhm_min=peak_fwhm_min, pa=(pa+90)) # conform to the function definition
    peak_temp = mock_map_elptc[np.round(y_map).astype(int),np.round(x_map).astype(int)]
    ratio = np.nanmedian(peak_val) / np.nanmedian(peak_temp)
    mock_map_elptc = mock_map_elptc * ratio

    # make the model map with fixed beam
    mock_map_beam = make_mock_image(
        Tmax, peak_ra, peak_dec, beam_fwhm,
        ra_map, dec_map, this_dist, peak_fwhm=beam_fwhm, mask_to_blank=(np.isfinite(mom0_map)==False))
    peak_temp = mock_map_beam[np.round(y_map).astype(int),np.round(x_map).astype(int)]
    ratio = np.nanmedian(peak_val) / np.nanmedian(peak_temp)
    mock_map_beam = mock_map_beam * ratio

    # Scale up the model map to conserve the flux
    flux_frac = np.nansum(mock_map_elptc) / np.nansum(mom0_map)
    print(flux_frac)
    flux_fracs.append(flux_frac)

    # Make the model map by adding the emission not captured by GMC catalogs
    ratio = np.nansum(mom0_map) * (1-flux_frac) / np.nansum(strict_mask_1kpc)
    mock_map_plus_bakgrd =  mock_map + strict_mask_1kpc * ratio
#     for i in range(2):
#         peak_temp = mock_map_plus_bakgrd[np.round(y_map).astype(int),np.round(x_map).astype(int)]
#         ratio = np.nanmedian(peak_val / peak_temp)
#         mock_map_plus_bakgrd = mock_map_plus_bakgrd * ratio
#         flux_frac = np.nansum(mock_map_plus_bakgrd) / np.nansum(mom0_map)
#         mock_map_plus_bakgrd = mock_map_plus_bakgrd + this_1kpc_map * (1-flux_frac) 

    # Write to disk
    model_file = mod_map_dir+this_source+'_source_model_'+reslt+'.fits'
    model_hdu = fits.PrimaryHDU(data=mock_map,header=hdr)
    model_hdu.writeto(model_file, overwrite=True)

    # write the model image with elliptical beam size
    model_file = mod_map_dir+this_source+'_source_model_elptc_'+reslt+'.fits'
    model_hdu = fits.PrimaryHDU(data=mock_map_elptc,header=hdr)
    model_hdu.writeto(model_file, overwrite=True)

    # write the model image with fixed beam size
    model_file = mod_map_dir+this_source+'_source_model_beam_'+reslt+'.fits'
    model_hdu = fits.PrimaryHDU(data=mock_map_beam,header=hdr)
    model_hdu.writeto(model_file, overwrite=True)

    # write the model image with backgroud emission
    model_file = mod_map_dir+this_source+'_source_model_bakgrd_'+reslt+'.fits'
    model_hdu = fits.PrimaryHDU(data=mock_map_plus_bakgrd, header=hdr)
    model_hdu.writeto(model_file, overwrite=True)

galaxy_fluxes = {}
galaxy_fluxes['Name'] = source_list
galaxy_fluxes['ratio'] = flux_fracs

with open('galaxy_fluxes.pkl', 'wb') as fp:
    pickle.dump(galaxy_fluxes, fp)
    print('dictionary saved successfully to file')
