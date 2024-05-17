import pandas as pd
import numpy as np

def clustering_uncertainty(source):
    
    # For a given galaxy, this function will retrieve the uncertainty on clustering amplitude curves#
    ### RETURNS ARRAY OF UPPER, LOWER LIMITS TO UNCERTAINTY ###
    
    fp = '/data/rubin/machado.35/phangs/products/clust_curves/' + source.lower() + '.csv'
    df = pd.read_csv(fp)
    upper = []
    lower = []
    for i in df.columns:
        upper.append(np.nanpercentile(df[i], 84))
        lower.append(np.nanpercentile(df[i], 16))
    
    return(np.array([upper,lower]))