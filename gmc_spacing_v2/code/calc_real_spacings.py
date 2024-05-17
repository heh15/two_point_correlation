import sys
sys.path.append('../')
sys.path.append('../aklpyutils/')
from key_dirs import *
from aklpyutils.utils_spacing import *

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Define tuning parameters
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Loop over catalogs
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

cat_dir_list = ['native_props']
cat_dir_res = {}
# cat_dir_res['native_props'] = None
 
for this_cat_dir in cat_dir_list:

    in_dir = cat_root_dir + this_cat_dir+'/'
    print(in_dir)
    
    # Find all catalogs in the directory
    cat_glob = glob.glob(in_dir+'*_props.fits')
    
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

        found_map = False
        for this_array in ['12m+7m+tp','7m+tp','7m']:
            if found_map == True:
                continue
            this_gal = gal_name
            if this_gal == 'eso097-013':
                this_gal = 'circinus'
            this_derived_dir = derived_root_dir + this_gal + '/'
            this_broad_map = this_derived_dir + this_gal + '_' + this_array + \
                '_co21_1kpc_broad_mom0.fits'
            if os.path.isfile(this_broad_map):
                found_map = True
        if found_map == False:
            print("Missing map for ", this_gal)
            raw = input("Stopped.")
        
        field_dict = {
            'ICO21_500PC': this_broad_map.replace('1kpc','500pc'),
            'EICO21_500PC': this_broad_map.replace('1kpc','500pc').replace('_mom0.fits','_emom0.fits'),
            'ICO21_1KPC': this_broad_map,
            'EICO21_1KPC': this_broad_map.replace('_mom0.fits','_emom0.fits'),
#            'ICO21_1P5KPC': this_broad_map.replace('1kpc','1p5kpc'),
#            'EICO21_1P5KPC': this_broad_map.replace('1kpc','1p5kpc').replace('_mom0.fits','_emom0.fits'),
        }
        
        add_maps_to_catalog(
            out_cat_file, out_cat_file, field_dict,
            overwrite=True)        
