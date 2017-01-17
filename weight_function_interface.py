from cosmosis.datablock import names, option_section
import numpy as np
import weight_function
from weight_function import compute_Wz
from weight_function import compute_wgp

# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.
cosmo = names.cosmological_parameters

def setup(options):
    # this function is called ONCE per processor per chain
    # use to load fixed parameters, perform one-off calculations (i.e. not changing
    # with sampling) etc.

    nz_section = options.get_string(option_section,'nz_section',default='wl_number_density')
    nbin = options[option_section,'nbin']
    wgp_section = options[option_section,'wgp_section']
    # INCLUDE A ZCUT ARG FOR HIGH- & LOW-Z WEIGHT-FNS

    # return the config for execute fn
    return nz_section,nbin,wgp_section

def execute(block, config):
    # this function is called every time you have a new sample of cosmological and other parameters

    nz_section,nbin,wgp_section = config

    # load from block
    z = block[nz_section,'z']
    dz = z[1]-z[0]
    wgp_rz = np.array([block[wgp_section,'bin_%s_%s'%(i+1,i+1)] for i in range(nbin)])
    nofz = block[nz_section,'nofz',nofz]
    # bin_popns = [block[nz_section,"bin_%s"%(i+1)] for i in range(nbin)]

    # compute W(z), wgp(r)
    Wz = compute_Wz(z,nofz)
    wgp_r = compute_wgp(Wz, wgp_rz)
    block.put_double_array_1d(wgp_section,'wgp(r_p)',wgp_r)
    block.put_double_array_1d(wgp_section,'W(z)',Wz)

    # return zero == all done, no probs
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass