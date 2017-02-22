from __future__ import print_function,division
import numpy as np
from astropy.cosmology import Planck13
import scipy.integrate as scint
import operator

def fivept_stencil(func,x,h):
    # returns f'(x), via 5pt stencil, for grid-spacing h
    return (-func(x+2*h)+8*func(x+h)-8*func(x-h)+func(x-2*h))/(12*h)

def compute_Wz(z,nofz_s,nofz_d,eta,beta,Rmag):
    # Wz = [p^2 / X^2*X'] / int[p^2 / X^2*X' dz]

    # compute p(z) = unconditional pdf
    pz_s = nofz_s/sum(nofz_s)
    pz_d = nofz_d/sum(nofz_d)
    assert pz_s.shape==z.shape, "p(z) vs. z mismatch"
    assert pz_d.shape==z.shape, "p(z) vs. z mismatch"

    # compute X(z) = comoving coordiante
    P13comov = lambda x: Planck13.comoving_distance(x)
    Xz = P13comov(z)
    Xz2 = Xz**2

    # compute X'(z) = first deriv.
    h = z[1]-z[0]
    Xprime = fivept_stencil(P13comov,z,h)

    # combine & integrate (Riemann sum) over z
    Wz_nom = (pz_s*pz_d)/(Xz2*Xprime)
    Wz_dom = sum(Wz_nom)*h
    Wz = Wz_nom/Wz_dom
    # redshift & luminosity dep. -> eta,beta free params when fitting!!
    # defaults; eta=-0.27, beta=1.13 (Joa+'11, all-sample-fit)
    zfactor = ((1+z)/(1+0.3))**eta # z0=0.3
    Lfactor = 10**(0.4*(-22-Rmag)*beta) # R0=-22
    Wz *= zfactor*Lfactor
    return Wz

def compute_wgp(Wz,wgp_rz,nbin,dz):
    # Riemann sum over W(z) for wgp(r,z) -> wgp(r)
    for i in range(nbin):
        wgp_rz[i] *= Wz[i]
    wgp_r = np.sum(wgp_rz,axis=0)*dz
    return wgp_r



