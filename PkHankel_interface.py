from cosmosis.datablock import names, option_section
import numpy as np
from PkHankel import read_z
from PkHankel import create_nz
from PkHankel import cut_krange
from PkHankel import zero_pad

# We have a collection of commonly used pre-defined block section names.
# If none of the names here is relevant for your calculation you can use any
# string you want instead.
cosmo = names.cosmological_parameters

def setup(options):
    # this function is called ONCE per processor per chain
    # use to load fixed parameters, perform one-off calculations (i.e. not changing
    # with sampling) etc.

    # prep P(k) for hankel transformation
    coerce_pk = options.get_bool(option_section,'coerce_pk',default=False)
    power_section = options.get_string(option_section,'power_section',default='dummy')
    nbin = options.get_int(option_section,'nbin',default=50)

    # create n(z) files for weight functions
    make_nz = options.get_bool(option_section,'make_nz',default=False)
    nz_section = options.get_string(option_section,'nz_section',default='dummy')
    shapz_path = options.get_string(option_section,'path1',default='dummy')
    densz_path = options.get_string(option_section,'path2',default='dummy')
    nz = options.get_int(option_section,'nz',default=50)
    outfile = options.get_string(option_section,'outfile',default='dummy')

    if make_nz:
        shapes_z,density_z = read_z(shapz_path,densz_path)
        z_mids,nofz_shap = create_nz(shapes_z,nz,nbin,outfile+'shapes.asc')
        z_mids,nofz_dens = create_nz(density_z,nz,nbin,outfile+'density.asc')
    else:
        z_mids,nofz,nz,nofz_shap,nofz_dens = (None,None,None,None,None)

    if not coerce_pk | make_nz:
        print('PKHANKEL IS DOING NOTHING !!')
        print('PkHankel can (1) coerce P(k) into hankel-transformable (cl_to_xi) format in datablock')
        print('and/or (2) create 2x n(z) 1d-arrays (2 samples) for computing weight fns')


    # return the config for execute fn
    return (power_section,nbin,make_nz,coerce_pk,z_mids,nofz_shap,nofz_dens,nz_section,nz)

def execute(block, config):
    # this function is called every time you have a new sample of cosmological and other parameters
    # collect config variables
    power_section,nbin,make_nz,coerce_pk,z_mids,nofz_shap,nofz_dens,nz_section,nz = config

    if coerce_pk:
        print('PkHankel: readying P(k) for hankel transform..')
        # load from datablock
        k_h = block[power_section,'k_h']
        p_k = block[power_section,'p_k']

        # take measures against possible ringing - cut k-range or zero-pad

        k_h,p_k = cut_krange(k_h, p_k, kmin=10**-3, kmax=10**3) 
        # start zero-padding at effklims, and stop at zeroklims
        k_h,p_k = zero_pad(k_h, p_k, effkmin=10**-3, effkmax=10**2, zerokmin=1e-4, zerokmax=10**5.7842)

        block[power_section,'ell'] = k_h
        block[power_section,'nbin'] = nbin
        for i in range(len(p_k)):
            block[power_section,'bin_%s_%s'%(i+1,i+1)] = p_k[i]

    if make_nz:
        print('PkHankel: making n(z) arrays for weight functions..')
        block.put_double_array_1d(nz_section,'z',z_mids)
        block.put_int(nz_section,'nz',nz)
        block.put_int(nz_section,'nbin',nbin)
        block.put_double_array_1d(nz_section,'nofz_shapes',nofz_shap)
        block.put_double_array_1d(nz_section,'nofz_density',nofz_dens)
        # for i in range(len(nofz)):
        #     bin_nz = np.zeros_like(nofz)
        #     bin_nz[i] = nofz[i]
        #     block.put_double_array_1d(nz_section,'bin_%s'%(i+1),bin_nz)

    # return zero == all done, no probs
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass
