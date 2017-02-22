from cosmosis.datablock import names, option_section
import numpy as np
import weight_function
from weight_function import compute_Wz
from weight_function import compute_wgp

cosmo = names.cosmological_parameters

def setup(options):
    # this function is called ONCE per processor per chain
    # use to load fixed parameters, perform one-off calculations (i.e. not changing
    # with sampling) etc.

    nz_section = options.get_string(option_section,'nz_section',default='wl_number_density')
    IA_section = options.get_string(option_section,'IA_section',default='intrinsic_alignment_parameters')
    wgp_section = options[option_section,'wgp_section']
    nbin = options[option_section,'nbin']
    Rmag = options.get_double(option_section,'Rmag')

    # return the config for execute fn
    return nz_section,IA_section,wgp_section,nbin,Rmag

def execute(block, config):
    # this function is called every time you have a new sample of cosmological and other parameters

    nz_section,IA_section,wgp_section,nbin,Rmag = config

    # load from block
    z = block[nz_section,'z']
    dz = z[1]-z[0]
    wgp_rz = np.array([block[wgp_section,'bin_%s_%s'%(i+1,i+1)] for i in range(nbin)])
    nofz_shap = block[nz_section,'nofz_shapes']
    nofz_dens = block[nz_section,'nofz_density']
    eta = block.get_double(IA_section,'eta',default=0)
    beta = block.get_double(IA_section,'beta',default=1.13)
    # bin_popns = [block[nz_section,"bin_%s"%(i+1)] for i in range(nbin)]

    # compute W(z), wgp(r)
    Wz = compute_Wz(z,nofz_shap,nofz_dens,eta,beta,Rmag)
    wgp_r = compute_wgp(Wz,wgp_rz,nbin,dz)
    block.put_double_array_1d(wgp_section,'wgp_r',wgp_r)
    block.put_double_array_1d(wgp_section,'wgp_r_minus',-wgp_r)
    block.put_double_array_1d(wgp_section,'W_z',Wz)

    # return zero == all done, no probs
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass