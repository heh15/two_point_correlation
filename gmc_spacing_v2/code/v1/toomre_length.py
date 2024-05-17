import numpy as np
import astropy.units as u
import astropy.constants as c

G = c.G.to(u.pc**3 / u.Msun / u.s**2)

def kappa2(beta, v_circ, r_gal):
    kappa = 2 * (1+beta) * (v_circ/r_gal)**2
    return(kappa)

def lambda_T(sigma_gal, kappa2):
    lT = 4*np.pi*G * sigma_gal / kappa2
    return(lT)