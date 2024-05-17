from key_dirs import *
from spacing_calcs import place_nclouds_in_map

dry_run = False

niters = 1e4
verbose = True

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Run random placements
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

this_res = res_list[0]

for this_source in source_list:
    
    for mask_type in ['flat','exp']:
            
        if mask_type == 'flat':
            this_image = mask_cov_dict[this_source]

        if mask_type == 'exp':
            this_image = mask_exp_dict[this_source]

        this_out_dir = sim_dir + mask_type + '/'
            
        this_cat = cat_dict[this_source]

        this_nclouds = ncloud_dict[this_source]
        
        this_row = deets_dict[this_source]
        this_ra = this_row['orient_ra']
        this_dec = this_row['orient_dec']
        this_pa = this_row['orient_posang']
        this_incl = this_row['orient_incl']
        this_dist = this_row['dist']
        
        this_outfile = this_out_dir+this_source+'_'+mask_type+'_'+str(this_res)+'pc.csv'

        n_clouds_to_place = this_nclouds*niters

        if verbose:
            print("Placing clouds:")
            print("... galaxy:", this_source)
            print("... mask:", this_image)
            print("... number of clouds:", n_clouds_to_place)
                  
        
        if dry_run == False:
            place_nclouds_in_map(
                this_image, n_clouds_to_place, this_outfile,
                this_ra, this_dec, this_pa, this_incl)
