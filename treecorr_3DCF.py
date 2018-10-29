from functions import *
import treecorr

def treecorr_3DCF(dataf, randf, config, estimator, nbins_rpar=30):
	"""
	dataf: paths to galaxy samples
	randf: paths to random points corresponding to galaxy samples
	config: dict specifying file types, column names/numbers, etc. a la TreeCorr configuration
	estimator: 'PW', 'AS' or 'wgg' to specify correlation & estimator
	nbins_rpar: number of line-of-sight bins for 3D correlation function -- specify limits in config arg
	"""
	assert estimator in ['PW', 'AS', 'wgg'], "for IA: estimator must be 'PW' (pair_weighted) or 'AS' (average shear), for clustering: 'wgg'"
	assert hasattr(dataf, '__iter__'), "dataf must be list/tuple of 2x paths; density, shapes for IA, or density1, density2 for clustering (these can be the same!)"
	assert hasattr(randf, '__iter__'), "randf must be list/tuple of 2x paths for randoms corresponding to each of dataf (these can be the same!)"

	if estimator in ['PW', 'AS']:
		corr = 'ng'
		XXCorrelation = treecorr.NGCorrelation
		gt_3D = np.zeros([len(Pi)-1, config['nbins']])
		gx_3D = np.zeros([len(Pi)-1, config['nbins']])
	elif estimator == 'wgg':
		corr = 'nn'
		XXCorrelation = treecorr.NNCorrelation
		wgg_3D = np.zeros([len(Pi)-1, config['nbins']])
	else:
		raise ValueError, "something went wrong with treecorr_3D!"

	Pi = np.linspace(config['min_rpar'], config['max_rpar'], nbins_rpar + 1)
	config_r = config.copy()
	config_r['flip_g1'] = config_r['flip_g2'] = False
	data1 = treecorr.Catalog(dataf[0], config)
	data2 = treecorr.Catalog(dataf[1], config)
	rand1 = treecorr.Catalog(randf[0], config_r, is_rand=1)
	rand2 = treecorr.Catalog(randf[1], config_r, is_rand=1)

	for p in range(len(Pi)-1):
		conf_pi = config.copy()
		conf_pi['min_rpar'] = Pi[p]
		conf_pi['max_rpar'] = Pi[p+1]

		# HOW to handle shot-noise errors?
		if corr == 'ng':
			ng = XXCorrelation(conf_pi)
			rg = XXCorrelation(conf_pi)
			ng.process_cross(data1, data2)
			rg.process_cross(rand1, data2)

			if estimator == 'PW':
				# these are RD_s pairs, not R_sR
				# may need to get R_sR from clustering...
				norm = rg.weight * float(data2.ntot) / float(rand1.ntot)
			elif estimator == 'AS':
				# norm by D_sD pairs
				norm = ng.weight

			gt_3D += -ng.xi / norm
			gx_3D += -ng.xi_im / norm

		elif corr == 'gg':
			nn = XXCorrelation(conf_pi)
			rr = XXCorrelation(conf_pi)
			nr = XXCorrelation(conf_pi)
			rn = XXCorrelation(conf_pi)

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
		gt = np.trapz(gt_3D, x=midpoints(Pi), axis=0)
		gx = np.trapz(gx_3D, x=midpoints(Pi), axis=0)
		r = ng.rnom
		return r, gt, gx
	elif corr == 'nn':
		wgg = np.trapz(wgg_3D, x=midpoints(Pi), axis=0)
		r = nn.rnom
		return r, wgg






