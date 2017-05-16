from cosmosis.datablock import names, option_section
import numpy as np
import weight_function
from weight_function import compute_Wz
from weight_function import compute_wgp

cosmo = names.cosmological_parameters
bias_section = names.bias_field

def setup(options):
    # this function is called ONCE per processor per chain
    # use to load fixed parameters, perform one-off calculations (i.e. not changing
    # with sampling) etc.

    # sections from which to read n(z), proj-fns, NLA power-laws(opt)
    nz_section = options.get_string(option_section,'nz_section',default='wl_number_density')
    IA_section = options.get_string(option_section,'IA_section',default='intrinsic_alignment_parameters')
    hkl_section = options.get_string(option_section,'projn_section',default='hankel_out')
    nbin = options.get_int(option_section,'nbin',default=50)

    nla = options.get_bool(option_section,'NLA',default=False) # only need bias factors if using NLA
    bias = options.get_string(option_section,'sample_bias',default='dummy') # None/dummy == using HOD
    Rmag = options.get_double(option_section,'Rmag',default=0.) # if fitting luminosity scaling
    wgg = options.get_bool(option_section,'do_wgg',default=False)

    wg_section = options[option_section,'wg_section'] # for weight/corrn fn outputs

    # return the config for execute fn
    return nz_section,IA_section,hkl_section,wg_section,bias,nbin,Rmag,nla,wgg

def execute(block, config):
    # this function is called every time you have a new sample of cosmological and other parameters

    nz_section,IA_section,hkl_section,wg_section,bias,nbin,Rmag,nla,wgg = config

    # load from block
    z = block[nz_section,'z']
    dz = z[1]-z[0]
    wg_rz = np.array([block[hkl_section,'bin_%s_%s'%(i+1,i+1)] for i in range(nbin)])
    nofz_shap = block[nz_section,'nofz_shapes']
    nofz_dens = block[nz_section,'nofz_density']
    if wgg:
        nofz_shap = nofz_dens
    eta = block.get_double(IA_section,'eta',default=0)
    beta = block.get_double(IA_section,'beta',default=0)

    # compute W(z), wg+(r)/wgg(r)
    Wz,Wz_scaled = compute_Wz(z,nofz_shap,nofz_dens,eta,beta,Rmag,wgg) # wgg switch prevents application of power-law scalings to gg
    wg_r = compute_wgp(Wz_scaled,wg_rz,nbin,dz)

    tags = ['wgp','wgg']
    if (bias!=None)&(nla):
        bg = block[bias_section,'b_g_%s'%bias]
        wg_r *= bg
        if wgg:
            print('COMPUTING w_gg; bias factor squared')
            wg_r *= bg
    block.put_double_array_1d(wg_section,'%s_r'%tags[int(wgg)],wg_r)
    block.put_double_array_1d(wg_section,'%s_r_minus'%tags[int(wgg)],-wg_r)
    block.put_double_array_1d(wg_section,'W_z',Wz)

    # return zero == all done, no probs
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass
