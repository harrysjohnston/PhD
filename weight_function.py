from __future__ import print_function,division
import numpy as np
from astropy.cosmology import Planck13
import scipy.integrate as scint
import operator

def fivept_stencil(func,x,h):
    # returns f'(x), via 5pt stencil, for grid-spacing h
    return (-func(x+2*h)+8*func(x+h)-8*func(x-h)+func(x-2*h))/(12*h)

def compute_Wz(z,nofz):
    # Wz = [p^2 / X^2*X'] / int[p^2 / X^2*X' dz]

    # compute p(z) = unconditional pdf
    pz = nofz/sum(nofz)
    assert pz.shape==z.shape, "p(z) vs. z mismatch"

    # compute X(z) = comoving coordiante
    P13comov = lambda x: Planck13.comoving_distance(x)
    Xz = P13comov(z)
    Xz2 = Xz**2

    # compute X'(z) = first deriv.
    h = z[1]-z[0]
    Xprime = fivept_stencil(P13comov,z,h)

    # combine & integrate (Riemann sum) over z
    Wz_nom = (pz**2)/(Xz2*Xprime)
    Wz_dom = sum(Wz_nom)*h
    Wz = Wz_nom/Wz_dom
    return Wz # WILL NEED TO EDIT FOR HIGH- & LOW-Z 

def compute_wgp(Wz,wgp_rz,nbin,dz):
    # Riemann sum over W(z) for wgp(r,z) -> wgp(r)
    for i in range(nbin):
        wgp_rz[i] *= Wz[i]
    wgp_r = np.sum(wgp_rz,axis=0)*dz
    return wgp_r



