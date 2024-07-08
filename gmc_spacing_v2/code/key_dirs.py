from astropy.table import Table
import glob, os
import numpy as np
from astropy.table import Table

# Directories for Machado spacing project.
root_dir = '../'
# data_dir = root_dir+'data/'
cat_root_dir =  root_dir+'cloud_cats_fits/'
fid_cat_dir =  root_dir+'cloud_cats_2024/gmc_cats_150pc_46mK/'
# fid_cat_dir =  cat_root_dir+'homogen_hi_90pc_props/'

# try to calculate the GMC auto-correlation function at 90 pc resolution (Feb. 27th, 2024)

cov_map_dir = root_dir+'data/coverage_maps/'
mod_map_dir = root_dir+'model_maps/'
# mod_map_dir = root_dir+'model_maps_90pc/'
mom0_map_dir = root_dir+'data/mom0/'
sim_dir =  root_dir+'sim_results/'
corr_dir = root_dir+'corr_results/'
# corr_dir = root_dir + 'corr_150pc_lo_results/'
derived_root_dir = '/vol/alcina/data1/hhe/two_point_correlation/delivery_v4p0/'
# Sample table
sample_tab_file = root_dir+'phangs_sample_table_v1p6.fits'
sample_tab = Table.read(sample_tab_file)

# change the name of eso097-013 to circinus
sample_tab['name'][0] = 'circinus'

# Get a list of our targets
temp_source_list = os.listdir(fid_cat_dir)
temp_source_list.sort()
source_list = []
for this_source_file in temp_source_list:
    this_source = this_source_file.split('_')[1]
    if this_source == 'ngc4945':
        continue
    if this_source == 'ngc0247':
        continue
    if this_source == 'ngc7743':
        continue
    source_list.append(this_source)

source_list = [*set(source_list)]
source_list.sort()
n_sources = len(source_list)    
print(source_list)

# Array strings
array_string_1 = '12m+7m+tp'

array_string_2 = '7m+tp'
# gals_array_string_2 = ['ngc0247','ngc0253','ngc0300','ngc4945','ngc5128','ngc7793']
gals_array_string_2 = ['ngc0247','ngc0253','ngc4945','ngc5128','ngc7793']

array_string_3 = '7m'
gals_array_string_3 = ['ngc0300']

# List of catalogs
# note lo means low sensitivity, high noise
cat_dir_list = [
    'homogen_lo_120pc_props',
    'homogen_hi_120pc_props',
    'homogen_lo_150pc_props',
    'homogen_hi_150pc_props',
    'homogen_lo_200pc_props',
    'homogen_hi_200pc_props',
    'homogen_lo_60pc_props',
    'homogen_hi_60pc_props',
    'homogen_lo_90pc_props',
    'homogen_hi_90pc_props',
    'matched_120pc_props',
    'matched_150pc_props',
#    'matched_200pc_props',
    'matched_60pc_props',
    'matched_90pc_props',
    'native_props',
]

cat_dir_res = {
    'homogen_lo_120pc_props':120.,
    'homogen_hi_120pc_props':120.,
    'homogen_lo_150pc_props':150.,
    'homogen_hi_150pc_props':150.,
#     'homogen_lo_200pc_props':200.,
#     'homogen_hi_200pc_props':200.,
    'homogen_lo_60pc_props':60.,
    'homogen_hi_60pc_props':60.,
    'homogen_lo_90pc_props':90.,
    'homogen_hi_90pc_props':90.,
    'matched_120pc_props':120.,
    'matched_150pc_props':150.,
#     'matched_200pc_props':200.,
    'matched_60pc_props':60.,
    'matched_90pc_props':90.,
    'native_props':None,
}

# Toggle output
verbose = True

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Prep items related to the fiducial catalogs
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

fid_cat_dict = {}
fid_cov_dict = {}
fid_ncloud_dict = {}
fid_ncloud_list = []

array_dict = {}
deets_dict = {}
mom0_map_dict = {}
broad_map_dict = {}

missing_cats = 0
missing_coverage = 0
    
exp_mask_ext = '_150pc_exp_mask.fits'

source_low_recfrac = []

for this_source in source_list:

    sample_tab_row = int(np.where(sample_tab['name'] == this_source)[0])
    deets_dict[this_source] = sample_tab[sample_tab_row]
                
    fidcat_glob = glob.glob(fid_cat_dir+'PHANGS_'+this_source+'*_gmccat.fits')
    if len(fidcat_glob) != 1:
        print("Unexpected number of catalogs: ", fidcat_glob)
        print("... for ", this_source)
    fid_cat_dict[this_source] = fidcat_glob[0]
    
    cloud_cat = Table.read(fidcat_glob[0], format='fits')
    nclouds = len(cloud_cat)
    fid_ncloud_list.append(nclouds)
    fid_ncloud_dict[this_source] = nclouds

    # exclude the galaxies with low GMC flux recovery fraction
    fid_flux_frac = np.nanmedian(cloud_cat['FRACF_CO21_OBJ'])
    if fid_flux_frac < 0.7:
        source_low_recfrac.append(this_source)

    this_cov_file = derived_root_dir+this_source+'/'+this_source+'_'+array_string_1+'_co21_150pc_coverage2d.fits'
#     this_cov_file = derived_root_dir+this_source+'/'+this_source+'_'+array_string_1+'_co21_90pc_coverage2d.fits'
    print(this_cov_file)
    if os.path.isfile(this_cov_file):
        array_dict[this_source] = array_string_1
        fid_cov_dict[this_source] = this_cov_file
    else:
        this_cov_file = derived_root_dir+this_source+'/'+this_source+'_'+array_string_2+'_co21_150pc_coverage2d.fits'
#         this_cov_file = derived_root_dir+this_source+'/'+this_source+'_'+array_string_2+'_co21_90pc_coverage2d.fits'
        if os.path.isfile(this_cov_file):
            array_dict[this_source] = array_string_2
            fid_cov_dict[this_source] = this_cov_file
        else:
            this_cov_file = derived_root_dir+this_source+'/'+this_source+'_'+array_string_3+'_co21_150pc_coverage2d.fits'
#             this_cov_file = derived_root_dir+this_source+'/'+this_source+'_'+array_string_3+'_co21_90pc_coverage2d.fits'
            if os.path.isfile(this_cov_file):
                array_dict[this_source] = array_string_3
                fid_cov_dict[this_source] = this_cov_file
            else:
                print("Missing coverage map for: ", this_source)
                missing_coverage += 1

    this_derived_dir = derived_root_dir + this_source + '/'
    this_fid_broad_map = this_derived_dir + this_source + '_' + array_dict[this_source] + \
        '_co21_1kpc_broad_mom0.fits'
    print(this_fid_broad_map)
    if os.path.isfile(this_fid_broad_map):
        broad_map_dict[this_source] = this_fid_broad_map
    else:
        print("Missing broad mom0 for ", this_source)

    this_fid_mom0_map = this_derived_dir + this_source + '_' + array_dict[this_source] + '_co21_150pc_broad_mom0.fits'
    # this_fid_mom0_map = this_derived_dir + this_source + '_' + array_dict[this_source] + '_co21_90pc_broad_mom0.fits'
    if os.path.isfile(this_fid_mom0_map):
        mom0_map_dict[this_source] = this_fid_mom0_map
    else:
        print("Missing fiducial mom0 for ", this_source)

source_list = [source for source in source_list if source not in source_low_recfrac]

#print(fid_ncloud_list)
ncloud_ra = np.array(fid_ncloud_list)
max_nclouds = np.nanmax(ncloud_ra)
