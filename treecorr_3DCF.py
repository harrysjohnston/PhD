from __future__ import division
import numpy as np
from scipy.integrate import simps
#from progress.bar import ChargingBar as Bar
from tqdm import tqdm
import treecorr
import pickle
midpoints = lambda x: (x[1:] + x[:-1]) / 2.

def compute_w(dataf, randf, config, estimator='PW1', compensated=1, nbins_rpar=30, random_oversampling=10., verbosity=1, largePi=0, **kwargs):
	"""
	dataf: paths to galaxy samples
	randf: paths to random points corresponding to galaxy samples
	config: path to config file, or dict specifying file types; column names/numbers, etc. a la TreeCorr configuration
	estimator: 'PW1', 'PW2', 'AS' or 'wgg' to specify correlation & estimator
	nbins_rpar: number of line-of-sight bins for 3D correlation function -- specify limits in config arg
	"""

	assert estimator in ['PW1', 'PW2', 'AS', 'wgg'], "for IA: estimator must be 'PW1/2' (pair_weighted 1=RDs, 2=RRs norm) or 'AS' (average shear), for clustering: 'wgg'"
	assert hasattr(dataf, '__iter__'), "dataf must be list/tuple of 2x paths; density, shapes for IA, or density1, density2 for clustering (these can be the same!)"
	assert hasattr(randf, '__iter__'), "randf must be list/tuple of 2x paths for randoms corresponding to each of dataf (these can be the same!)"
	nbins_rpar = int(nbins_rpar)
	random_oversampling = float(random_oversampling)
	verbosity = int(verbosity)
	largePi = int(largePi)

	if type(config) == str:
		config = treecorr.read_config(config)
	if config['file_type'] == 'ASCII':
		config['ra_col'] =  int(config['ra_col'])
		config['dec_col'] = int(config['dec_col'])
		config['r_col'] =   int(config['r_col'])
		config['g1_col'] =  int(config['g1_col'])
		config['g2_col'] =  int(config['g2_col'])
	config['verbose'] = verbosity

	if not largePi:
		Pi = np.linspace(config['min_rpar'], config['max_rpar'], nbins_rpar + 1)
	else:
		dPi = (config['max_rpar'] - config['min_rpar']) / nbins_rpar
		Pi = np.arange(config['min_rpar']*1.5, config['max_rpar']*1.5 + dPi, step=dPi, dtype=float)

	if estimator in ['PW1', 'PW2', 'AS']:
		corr = 'ng'
		gt_3D = np.zeros([len(Pi)-1, config['nbins']])
		gx_3D = np.zeros([len(Pi)-1, config['nbins']])
		varg_3D = np.zeros([len(Pi)-1, config['nbins']])
	elif estimator == 'wgg':
		corr = 'nn'
		wgg_3D = np.zeros([len(Pi)-1, config['nbins']])
	else:
		raise ValueError, "unsupported estimator choice"

	config_r = config.copy()
	config_r['flip_g1'] = False
	config_r['flip_g2'] = False
	data1 = treecorr.Catalog(dataf[0], config) # 1 = density/lenses
	data2 = treecorr.Catalog(dataf[1], config) # 2 = shapes
	rand1 = treecorr.Catalog(randf[0], config_r, is_rand=1)
	rand2 = treecorr.Catalog(randf[1], config_r, is_rand=1)
	f1 = data1.ntot * random_oversampling / float(rand1.ntot)
	f2 = data2.ntot * random_oversampling / float(rand2.ntot)
	rand1.w = np.array(np.random.rand(rand1.ntot) < f1, dtype=float)
	rand2.w = np.array(np.random.rand(rand2.ntot) < f2, dtype=float)
	varg = treecorr.calculateVarG(data2)

	for p in tqdm(range(len(Pi)-1), ascii=True, desc='Correlating'):

		if largePi & any(abs(Pi[p:p+2]) < config['max_rpar']):
			continue

		conf_pi = config.copy()
		conf_pi['min_rpar'] = Pi[p]
		conf_pi['max_rpar'] = Pi[p+1]

		if corr == 'ng':
			ng = treecorr.NGCorrelation(conf_pi)
			rg = treecorr.NGCorrelation(conf_pi)
			ng.process_cross(data1, data2)
			rg.process_cross(rand1, data2)

			if estimator == 'PW1': # RDs norm
				f = data1.ntot / rand1.w.sum()
				norm1 = rg.weight * f
				norm2 = rg.weight
			if estimator == 'PW2': # RRs norm
				RRs = get_RRs(rand1, rand2, conf_pi, **kwargs)
				f1 = data1.ntot / rand1.w.sum()
				f2 = data2.ntot / rand2.w.sum()
				norm1 = RRs * f1 * f2
				norm2 = RRs * f2
			elif estimator == 'AS': # DDs norm
				norm1 = ng.weight
				norm2 = rg.weight

			if int(compensated):
				gt_3D[p] += (ng.xi / norm1) - (rg.xi / norm2)
				gx_3D[p] += (ng.xi_im / norm1) - (rg.xi_im / norm2)
				varg_3D[p] += (varg / norm1) + (varg / norm2)
			else:
				gt_3D[p] += ng.xi / norm1
				gx_3D[p] += ng.xi_im / norm1
				varg_3D[p] += varg / norm1

		elif corr == 'nn':
			nn = treecorr.NNCorrelation(conf_pi)
			rr = treecorr.NNCorrelation(conf_pi)
			nr = treecorr.NNCorrelation(conf_pi)
			rn = treecorr.NNCorrelation(conf_pi)

			if dataf[0] == dataf[1]:
				nn.process(data1)
				rr.process(rand1)
				nr.process(data1, rand1)
				xi, varxi = nn.calculateXi(rr, nr)
			else:
				nn.process(data1, data2)
				rr.process(rand1, rand2)
				nr.process(data1, rand2)
				rn.process(rand1, data2)
				xi, varxi = nn.calculateXi(rr, nr, rn)

			wgg_3D[p] += xi

	if corr == 'ng':
		#gt = np.trapz(gt_3D, x=midpoints(Pi), axis=0)
		#gx = np.trapz(gx_3D, x=midpoints(Pi), axis=0)
		#gt = simps(gt_3D, x=midpoints(Pi), axis=0)
		gt = np.sum(gt_3D * (Pi[1] - Pi[0]), axis=0)
		gx = np.sum(gx_3D * (Pi[1] - Pi[0]), axis=0)
		varg = np.sum(varg_3D, axis=0)
		r = ng.rnom
		return r, gt, gx, varg**0.5
	elif corr == 'nn':
		wgg = np.sum(wgg_3D * (Pi[1] - Pi[0]), axis=0)
		r = nn.rnom
		return r, wgg

def get_RRs(R_cat, Rs_cat, config, load_RRs=None, save_RRs=None, **kwargs):
	minmax = (config['min_rpar'], config['max_rpar'])

	if load_RRs is not None:
		# read paircounts from pickle file
		RRs_dict = pickle.load(open(load_RRs, 'r'))
		rrs_weight = RRs_dict[minmax]
		assert len(rrs_weight) == config['nbins'], "RRs pairs from file do not match specified transverse binning"
		return rrs_weight
	else:
		# compute paircounts and create pickle file
		rrs = treecorr.NNCorrelation(config)
		rrs.process_cross(R_cat, Rs_cat)
		if save_RRs is not None:
			try:
				RRs_dict = pickle.load(open(save_RRs, 'r'))
			except IOError:
				RRs_dict = {}
			RRs_dict[minmax] = rrs.weight
			pickle.dump(RRs_dict, open(save_RRs, 'w'))
		return rrs.weight




