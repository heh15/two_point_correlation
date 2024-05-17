import sys
sys.path.append('../aklpyutils/')
import glob

# import key_dirs
from key_dirs import *
from utils_spacing import *

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Define bins and other tuning parameters
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

bins = np.linspace(0,4400,44)
nbins = len(bins)

do_sims = False

# this_res = res_list[0]

thresh_cloud_count = 5

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Loop over catalogs
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Find all catalogs in the directory
cat_glob = glob.glob(fid_cat_dir+'*_value_added.fits')

for in_cat_file in cat_glob:

    # Remap to an outfile
    out_cat_file = in_cat_file.replace('_props.fits','_props_valueadded.fits')

    gal_name = in_cat_file.split(in_dir)[1].split('_')[0]
    
    print("Calculating spacing for catalogs: ")
    print("... from: ", in_cat_file)
    print("... to: ", out_cat_file)
    print("... galaxy: ", gal_name)

    sample_tab_row = int(np.where(sample_tab['name'] == gal_name)[0])

    this_row = sample_tab[sample_tab_row]

    this_ra = float(this_row['orient_ra'])
    this_dec = float(this_row['orient_dec'])
    this_pa = float(this_row['orient_posang'])
    this_incl = float(this_row['orient_incl'])
    this_dist = float(this_row['dist'])

    add_spacings_to_catalog(
        in_cat_file, out_cat_file,
        this_ra, this_dec, this_pa, this_incl, this_dist,
        use_peak=True,
        overwrite=True)
 
