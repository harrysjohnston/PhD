from cosmosis.datablock import names, option_section
import sys
import gc
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

    # compute W(z) and skip computation of wp(rp)
	compute_wp = options.get_bool(option_section, 'compute_wp', default=True)
	wz_section = options.get_string(option_section, 'wz_section', default='')
	if not compute_wp: assert wz_section != '', "specify W(z) out-section if skipping wp computation; wz_section=______"

    # sections from which to read n(z), proj-fns, NLA power-laws(opt)
	nz_section = options.get_string(option_section,'nz_section',default='wl_number_density')
	IA_section = options.get_string(option_section,'IA_section',default='intrinsic_alignment_parameters')
	hkl_section = options.get_string(option_section,'projn_section',default='hankel_out')
	nbin = options.get_int(option_section,'nbin',default=50)

	nla = options.get_bool(option_section,'NLA',default=False) # only need bias factors if using NLA
	skey = options.get_string(option_section,'shapes_key',default='dummy')
	dkey = options.get_string(option_section,'dens_key',default=skey)
	HOD_bias_loc = options.get_string(option_section,'HOD_bias_loc',default='dummy/dummy') # hod_section/name of HOD generated bias(z) for manual appln to matter-int spectra
	HOD_bias_loc = HOD_bias_loc.split('/')
	if nla & (HOD_bias_loc[0]!='dummy'):
		print('HOD/HM & NLA conflict !!')
		print('check NLA flags and HOD_bias_loc args')
		sys.exit()

	wgg = options.get_bool(option_section,'do_wgg',default=False)
	wg_section = options[option_section,'wg_section'] # for weight/corrn fn outputs

	# return the config for execute fn
	config = {
		'nz_section':nz_section,
		'IA_section':IA_section,
		'hkl_section':hkl_section,
		'wg_section':wg_section,
		'skey':skey,
		'dkey':dkey,
		'nbin':nbin,
		'nla':nla,
		'wgg':wgg,
		'HOD_bias_loc':HOD_bias_loc,
		'compute_wp':compute_wp,
		'wz_section':wz_section
		}

	return config

def execute(block, config):
    # this function is called every time you have a new sample of cosmological and other parameters

	for key in config.keys():
		exec '%s = config["%s"]'%(key, key)

    # load from block
	z = block[nz_section,'z']
	dz = z[1]-z[0]
	nofz_shap = block[nz_section,'nofz_shapes']
	nofz_dens = block[nz_section,'nofz_density']
	if wgg:
		nofz_shap = nofz_dens

	eta = block.get_double(IA_section,'eta',default=0)

	# run beta_red/blue
	if skey.endswith('R'):
		beta = block[IA_section, 'beta_R']
	elif skey.endswith('B'):
		beta = block[IA_section, 'beta_B']
	else:
		print('skey must be colour-specific!!')
		sys.exit()

	Rmag = block.get_double(IA_section, 'Rmag_%s'%skey, default=0.)

	# compute W(z)
	Wz,Wz_scaled = compute_Wz(z,nofz_shap,nofz_dens,eta,beta,Rmag,wgg)
	# ^wgg switch prevents application of power-law scalings to gg

	if not compute_wp:
		if nla:
			block[wz_section, 'w_z'] = Wz_scaled
		else:
			block[wz_section, 'w_z'] = Wz
		block[wz_section, 'z'] = z

	else:
		# compute wp(rp)
		wg_rz = np.array([block[hkl_section,'bin_%s_%s'%(i+1,i+1)] for i in range(nbin)])
		if HOD_bias_loc[0] != 'dummy':
			b_z = block[HOD_bias_loc[0], HOD_bias_loc[1]]
			assert len(wg_rz)==len(b_z), 'wgp(r, z) vs. HOD bias(z) length mismatch'
			wg_rz = (wg_rz.T * b_z).T
			if wgg:
				wg_rz = (wg_rz.T * b_z).T
		wg_r = compute_wgp(Wz_scaled,wg_rz,nbin,dz)

		tags = ['wgp','wgg']
		if nla:
			if skey.endswith('R'):
				A_i =  block[IA_section,'A_R']
			elif skey.endswith('B'):
				A_i =  block[IA_section,'A_B']
			else:
				print('skey must be colour-specific!!')
				sys.exit()

			if not wgg:
				wg_r *= A_i

			bg = block[bias_section,'b_g_%s'%dkey]
			wg_r *= bg

			if wgg:
				C = block[bias_section,'ic_%s'%dkey]
				wg_r *= bg
				wg_r += C

		elif wgg:
				C = block[bias_section,'ic_%s'%dkey]
				wg_r += C
		block.put_double_array_1d(wg_section,'%s_r'%tags[int(wgg)],wg_r)
		block.put_double_array_1d(wg_section,'%s_r_minus'%tags[int(wgg)],-wg_r)
		block.put_double_array_1d(wg_section,'W_z',Wz)

	# return zero == all done, no probs
	return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass




