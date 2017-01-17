from cosmosis.datablock import names, option_section
import numpy as np
from PkHankel import read_z
from PkHankel import create_nz
from PkHankel import cut_krange

# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.
cosmo = names.cosmological_parameters

def setup(options):
    # this function is called ONCE per processor per chain
    # use to load fixed parameters, perform one-off calculations (i.e. not changing
    # with sampling) etc.

    power_section = options[option_section,'power_section'] # options[] must be in .ini
    nbin = options[option_section,'nbin']
    make_nz = options.get_bool(option_section,'make_nz',default=False)
    nz_section = options.get_string(option_section,'nz_section',default=None)
    cat_path = options.get_string(option_section,'catalog_path',default=None)
    nz = options.get_int(option_section,'nz',default=None)
    outfile = options.get_string(option_section,'outfile',default=None)

    # optionally, generate ascii & save info on n(z) into datablock
    if make_nz:
        gal_z = read_z(cat_path)
        z_mids,nofz = create_nz(gal_z,nz,nbin,outfile)
    else:
        z_mids,nofz,nz = (None,None,None)

    # return the config for execute fn
    return (power_section,nbin,make_nz,z_mids,nofz,nz_section,nz)

def execute(block, config):
    # this function is called every time you have a new sample of cosmological and other parameters
    # collect config variables
    power_section,nbin,make_nz,z_mids,nofz,nz_section,nz = config

    # load from datablock
    k_h = block[power_section,'k_h']
    p_k = block[power_section,'p_k']

    # execute main function
    # aim for simplicity
    # here, cutting k,P(k) & saving back to db in format for Hankel transfm
    k_h,p_k = cut_krange(k_h,p_k, kmin=10**-2.2, kmax=10**1.2) # k-limits in h/Mpc
    block.put_double_array_1d(power_section,'ell',k_h)
    block.put_int(power_section,'nbin',nbin)
    for i in range(len(p_k)):
        block.put_double_array_1d(power_section,'bin_%s_%s'%(i+1,i+1),p_k[i])

    if make_nz:
        block.put_double_array_1d(nz_section,'z',z_mids)
        block.put_int(nz_section,'nz',nz)
        block.put_int(nz_section,'nbin',nbin)
        for i in range(len(nofz)):
            bin_nz = np.zeros_like(nofz)
            bin_nz[i] = nofz[i]
            block.put_double_array_1d(nz_section,'bin_%s'%(i+1),bin_nz)

    # return zero == all done, no probs
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass