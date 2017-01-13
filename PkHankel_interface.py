from cosmosis.datablock import names, option_section
from PkHankel import read_z
from PkHankel import create_nz

# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.
cosmo = names.cosmological_parameters

def setup(options):
    # This function is called ONCE per processor per chain
    # use to load fixed parameters, perform one-off calculations (i.e. not changing
    # with sampling) etc.

    power_section = options[option_section,'power_section'] # options[] must be in .ini
    nbin = options[option_section,'nbin']
    make_nz = options.get_bool(option_section,'make_nz',default=False)
    cat_path = options.get_string(option_section,'catalog_path',default=None)
    nz = options.get_int(option_section,'nz',default=None)
    outfile = options.get_string(option_section,'outfile',default=None)

    # optionally, generate n(z) ascii file for later use
    if make_nz:
        gal_z = read_z(cat_path)
        create_nz(gal_z,nz,nbin,outfile)

    # return the config for execute fn
    return (power_section,nbin)

def execute(block, config):
    # This function is called every time you have a new sample of cosmological and other parameters
    
    # redefine config variables
    power_section,nbin = config

    # load from datablock
    k_h = block[power_section,'k_h']
    p_k = block[power_section,'p_k']

    # execute main function
    # aim for simplicity
    # here, saving k,P(k) back to db in format for Hankel transfm
    block.put_double_array_1d(power_section,'ell',k_h)
    block.put_int(power_section,'nbin',nbin)
    for i in xrange(p_k):
        block.put_double_array_1d(power_section,'bin_%s_%s'%(i,i),p_k[i])

    # return zero == all done, no probs
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass