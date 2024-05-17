import numpy as np
import astropy.units as u
import astropy.constants as c

def jeans_len(vel, rho):
    '''
    Requires inputs:
        vel: km/s
        rho: M_sun / pc^3
    '''
    G = c.G.to(u.pc**3 / u.solMass / u.s**2)
    
    lambda_j = np.sqrt((np.pi * vel.to(u.pc/u.s)**2) / (G * rho))
    
    return(lambda_j)