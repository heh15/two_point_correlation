from astropy.table import Table
from os import listdir

# Directories for Machado spacing project.
root_dir = '../'
data_dir = root_dir+'data/'
cat_dir =  root_dir+'data/cloud_cats_csv/'

# Sample table
sample_tab_file = data_dir+'phangs_sample_table_v1p6.fits'
t = Table.read(sample_tab_file)

# Get a list of our targets
sources = listdir(cat_dir)
for ii in range(len(sources)):
    sources[ii] = sources[ii].split('_')[0]

# Note our resolutions
res = [150]
