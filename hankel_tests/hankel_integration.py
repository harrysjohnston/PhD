from hankel_functions import *
from cosmosis.datablock import names, option_section
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import matplotlib.colors as colors
from argparse import Namespace
from time import time
import os

def setup(options):
# FOR MAIN ANALYSIS, NEED TO INCLUDE APPLICATION OF COLOUR_SPEC A_IA
	in_section = options.get_string(option_section, 'in_section')
	out_section = options.get_string(option_section, 'out_section')
	bessel_nu = options.get_int(option_section, 'bessel_nu', default=2)
	#mode = options.get_string(option_section, 'mode', default='3D')
	stepsize_h = options.get_double(option_section, 'stepsize_h', default=0.)
	r_min, r_max = [float(ir) for ir in options.get_string(option_section, 'r_lims', default='0.07, 80.').strip(' ').split(',')]
	num_r = options.get_int(option_section, 'num_r', default=20)
	abs_Pi_max = options.get_double(option_section, 'abs_Pi_max', default=60.)
	num_Pi = options.get_int(option_section, 'num_Pi', default=30)
	Pisymm = options.get_bool(option_section, 'Pisymm', default=True)
	norm = options.get_string(option_section, 'norm', default='-1./(2*np.pi**2)')
	nzeros = options.get_int(option_section, 'nzeros', default=10)
	kpar_max = options.get_double(option_section, 'kpar_max', default=0.)
	f_kpar = options.get_double(option_section, 'f_kpar', default=1.)
	f_kperp = options.get_double(option_section, 'f_kperp', default=1.)
	intz = options.get_bool(option_section, 'intz', default=False)
	nz_sec = options.get_string(option_section, 'nz_section')
	if intz:
		assert nz_sec is not None, ("MUST provide n(z) section name, containing "
							"n(z) per sample, if integrating over redshift (intz=T)")
	div_sec = options.get_string(option_section, 'div_section', default='none')

	config = {
		'in_sec':in_section,
		'out_sec':out_section,
		'nu':bessel_nu,
		'h':stepsize_h,
		'r_min':r_min,
		'r_max':r_max,
		'num_r':num_r,
		'abs_Pi_max':abs_Pi_max,
		'num_Pi':num_Pi,
		'Pisymm':Pisymm,
		'norm':norm,
		'nzeros':nzeros,
		'kpmax':kpar_max,
		'fkpar':f_kpar,
		'fkperp':f_kperp,
		'intz':intz,
		'nz_sec':nz_sec,
		'div_sec':div_sec
			}
	return config

def execute(block, config):
	# Namespace object for options: cf.<arg>
	cf = Namespace(**config)
	# take input power spectrum
	k = block[cf.in_sec, 'k_h']
	pk = block[cf.in_sec, 'p_k']
	z = block[cf.in_sec, 'z']
	if cf.nzeros != 0:
		k, pk = zero_pad(k, pk, nzeros=cf.nzeros)

	# prep integrator
	if cf.h == 0:
		if block.has_value('hankel', 'h'):
			cf.h = block['hankel', 'h']
		else:
			# try to optimise step-size h
			# this method is NOT always reliable, or even good...
			fpk = interp1d(k, pk[0], fill_value=0., bounds_error=0)
			cf.h = choose_h(fpk, cf.nu, cf.r_min, cf.r_max)
			block.put('hankel', 'h', cf.h)
			print('optimal hankel stepsize h determined = %f'%cf.h)
	HT = hankel.HankelTransform(nu=cf.nu, h=cf.h)#, N=10000)

	if cf.intz:
		# compute redshift weight-fn, skip W(z) = 0 loops
		nz1 = block[cf.nz_sec, 'nofz_shapes']
		nz2 = block[cf.nz_sec, 'nofz_density']
		nz_z = block[cf.nz_sec, 'z']
		W = compute_Wz(nz_z, nz1, nz2)
		#W = np.interp(z, nz_z, W)

	# define real- & fourier-space boundaries
	rp = np.logspace(log10(cf.r_min), log10(cf.r_max), cf.num_r)
	Pi = np.linspace(-cf.abs_Pi_max, cf.abs_Pi_max, cf.num_Pi+1)
	Pi = (Pi[1:] + Pi[:-1]) / 2.
	if cf.Pisymm:
		Pi = Pi[len(Pi)/2:]

	kpar, kperp = k.copy(), k.copy()
	#kpar = upsample_k(kpar, 1, 'lin', n=4e4)
	kperp = k_scale_cuts(kperp, 1e-3, 1e5)
	kperp = upsample_k(kperp, 1, n=500)
	kpar = upsample_k(kpar, 1, n=500)
	#kpar = upsample_k(kpar, 1, 'lin', n=2**16 + 1)
	#kgrid, cosine = prep_kpar_integrand(kpar, kperp, Pi)

	if cf.kpmax == 0:
		# set kpar_max = inf
		cf.kpmax = choose_kpar_max(Pi)
	else:
		kpar = k_scale_cuts(kpar, 0, cf.kpmax)

	if cf.fkpar != 1:
		kpar = upsample_k(kpar, cf.fkpar, 'lin')
	if cf.fkperp != 1:
		kperp = upsample_k(kperp, cf.fkperp)

	xi = np.zeros([len(pk), len(Pi), len(rp)])
	# looping over redshifts
	for iz, ipk in enumerate(pk):
		# if W(z) == 0; skip
		if cf.intz:
			if W[iz] == 0:
				xi[iz] = np.zeros_like(xi[0])
				continue
			else:
				print('z=%.2f'%z[iz])

		ti = time()
		#os.system('rm /share/data1/hj/millennium/hankel_tests/order%s_Pi*_curve.txt'%cf.nu)
			#tk1 = time()

		for ip, pi in enumerate(Pi):
			#fname = '/share/data1/hj/millennium/hankel_tests/order%s_Pi%i_curve.txt'%(cf.nu,pi)
			Ptilde_kperp = []
			for ikpp, kpp in enumerate(kperp):
				#tp1 = time()
				kpar_integrand = make_kpar_integrand_1d(k, ipk, kpp)#, pi)
				Ptilde, err = quad(kpar_integrand, 0., np.inf, weight='cos', wvar=pi)#, points=np.arange(100000)[1:]/(pi*2.))#, full_output=1)
				#f = lambda kpp_: interp1d(k, ipk, fill_value=0., bounds_error=0)(kpp_) * Ptilde
				#with open(fname, 'a') as scr:
				#	scr.write('%f\t%f\n'%(kpp, Ptilde))
				Ptilde_kperp.append(Ptilde)

			Ptilde_kperp = np.array(Ptilde_kperp)
			f = interp1d(kperp, Ptilde_kperp, fill_value='extrapolate', bounds_error=0)
			xi[iz, ip] += HT.transform(f, rp, ret_err=False)
				#tp2 = time() - tp1
				#print '{kperp}Pi-loop time: ', tp2
			#tk2 = time() - tk1
			#print 'kperp-loop time: %.3f'%tk2, '-> ~%.1f for %s kperp values'%(tk2*len(kperp), len(kperp))
		#ttot = time() - ti
		#print 'total integration time: %.3f'%ttot
				

#		ti = time() ; ttot = 0.
#		kpar_integrand = make_kpar_integrand_2d(k, ipk, *kpar_integrand_args)
#		tf = time() ; ttot += tf - ti
#		print('time to compute kpar_integrand: %.3f'%(tf-ti))
#		kperp_integrand = np.trapz(kpar_integrand, x=kpar, axis=1)
#		ti = tf ; tf = time() ; ttot += tf - ti
#		print('time to compute (trapz) kperp_integrand: %.3f'%(tf-ti))

		# build interpolator
		#KPIinterp = interp2d(kperp, Pi, kperp_integrand.T, fill_value=0., bounds_error=0.)

		# integrate over k_perpendicular
#		def LookUp(kperp, pi):
#			if kperp.ndim == 2:
#				# hankel module may want to pass
#				# 2darray with shape = (rp-out, kperp-in)
#				return np.array([KPIinterp(kpi, pi) for kpi in kperp]).squeeze()
#			else:
#				return KPIinterp(kperp, pi)
#
#		#for ip in range(len(Pi)):
#			f = lambda kperp: LookUp(kperp, pi)
#			xi[iz, ip] += HT.transform(f, rp, ret_err=False)
#
#		ti = tf ; tf = time() ; ttot += tf - ti
#		print('time to hankel transform xi-per-Pi: %.3f'%(tf-ti))
#		print('total per z-loop: %.3f'%ttot)

	# save 3D correlation function to block
	norm = eval(cf.norm)
	xi *= norm# / rpow
	if cf.Pisymm:
		xi_negPi = xi.copy()
		xi = np.concatenate((xi[:, ::-1], xi_negPi), axis=1)
		Pi = np.concatenate((-Pi[::-1], Pi))
	block.put(cf.out_sec, 'xi', xi)

	# divide CF by 1+xi_gg clustering CF, if appl.
	if cf.div_sec != 'none':
		assert block.has_value(cf.div_sec, 'xi'), "MUST have (sec='div_sec', val='xi') in the block for division by 1 + xi -- did you run this module already?"
		xi_gg = block[cf.div_sec, 'xi']
		assert xi_gg.shape == xi.shape, "(z, Pi, rp) dimensions mismatched between spectra!"
		xi = xi / (1. + xi_gg)
		block.put(cf.out_sec, 'xi_div', xi)

	# sum over Pi & save -- any divisions (above) are present henceforth
	xi_proj = np.trapz(xi, x=Pi, axis=1)
	block.put(cf.out_sec, 'rp', rp)
	block.put(cf.out_sec, 'pi', Pi)
	block.put(cf.out_sec, 'w_rpz', xi_proj)

	if cf.intz:
		# integrate over redshift & save
		xi_proj = np.trapz(xi_proj.T * W, x=z, axis=1)
		block.put(cf.out_sec, 'w_rp', xi_proj)
		block.put(cf.out_sec, 'w_z', W)
		block.put(cf.out_sec, 'z', z)

	return 0

def cleanup(config):
	pass

#MODE = ['IA', 'CLUS'] [0]
#if MODE=='CLUS':
#	rid = 'ells'
#	section = 'galaxy_power'
#	NU = 0
#	AR = 1.
#	prefactor = 1. / (2.*pi**2.) * AR
#if MODE=='IA':
#	rid = 'v17_ells'
#	section = 'galaxy_intrinsic_power'
#	NU = 2
#	AR = 7.
#	prefactor = -1. / (2.*pi**2.) * AR
#
#h = 0.01
#ht = hankel.HankelTransform(nu=NU, h=h)
#k = np.loadtxt('./%s_%s/k_h.txt'%(section, rid))
#pkz = np.loadtxt('./%s_%s/p_k.txt'%(section, rid))
#zbase = np.loadtxt('./%s_%s/z.txt'%(section, rid))
#
#Wz = np.loadtxt('./wgp_z2_r/w_z.txt')[:len(pkz)]
#true_wgp = np.loadtxt('./wgp_z2_r/wgp_r_minus.txt')
#m06_wgp = np.loadtxt('./wgp_z2_b/wgp_r_minus.txt')
#th = np.loadtxt('./theta.txt')
#
#rp = np.logspace(-1.5, 2, 30)
#Pi = np.linspace(-60., 60., 31)
#Pi = (Pi[1:] + Pi[:-1]) / 2.
#
#kcut = (-np.inf, np.inf)
#
#ZERO_PAD = 1
#
#if ZERO_PAD:
#	nzeros = 10
#	kl = 10**(np.log10(k[0]) - np.arange(1, nzeros+1) * np.log10(k[1]/k[0]))[::-1]
#	ku = 10**(np.log10(k[-1]) + np.arange(1, nzeros+1) * np.log10(k[-1]/k[-2]))
#	p0 = np.zeros_like(ku)
#	k = np.concatenate((kl, k, ku))
#	pkz_ = pkz.copy()
#	pkz = []
#	for i in range(len(pkz_)):
#		pkz.append(np.concatenate((p0, pkz_[i], p0)))
#	pkz = np.array(pkz)
#
#HTEST = 0
#KPTEST = 0
#
#def k_scale_cuts(kx, low, high):
#	return kx[(kx >= low) & (kx <= high) & (kx >= k.min()) & (kx <= k.max())]
#def upsample_k(kx, factor):
#	return np.logspace(np.log10(kx.min()), np.log10(kx.max()), int(factor * len(kx)))
#
#def run(kpar_max=0.031, plots=(1, 1)):
#	kperp = k.copy()
#	kpar = kperp.copy()
#	kperp, kpar = map(lambda x: k_scale_cuts(x, kcut[0], kcut[1]), [kperp, kpar])
#	# cut large-ks
#	kpar = k_scale_cuts(kpar, 0., kpar_max)
#	# upsample kpar for approx integrators
#	#kpar = upsample_k(kpar, 0.5)
#	#kperp = upsample_k(kperp, 0.2)
#	kpar_m, Pi_m = np.meshgrid(kpar, Pi)
#	xpar = kpar_m * Pi_m
#	k_grid = np.sqrt(np.add.outer(kperp**2, (xpar/Pi_m)**2))
#
#	pkzz = pkz.copy()
#	xi = np.zeros([len(pkz), len(Pi), len(rp)])
#
#	for z, pk in enumerate(pkzz):
#
#		if Wz[z] == 0:
#			xi[z] = np.zeros_like(xi[0])
#			continue
#
#		PKinterp = interp1d(k, pk, bounds_error=0, fill_value=0.)
#		pk_grid = PKinterp(k_grid)
#		#import pdb ; pdb.set_trace()
#		#cosine_zeros = np.abs(pi / (2. * Pi))
#		#kpar_inds = [np.argwhere(kpar > cz)[0][0] for cz in cosine_zeros]
#		#kperp_term_grid = np.array([simps(kpar_term_grid[:, p, kpi:], x=xpar[p][kpi:], axis=1) for p, kpi in zip(range(len(Pi)), kpar_inds)])
#
#		kpar_term_grid = np.cos(xpar) * pk_grid
#		kperp_term_grid = np.array([simps(kpar_term_grid[:, p, :], x=xpar[p], axis=1) for p in range(len(Pi))])
#		kperp_term_grid = (kperp_term_grid.T / Pi).T
#		interp_kperp_ints = interp2d(kperp, Pi, kperp_term_grid)
#
#		def LookUp(kperp, pi):
#			if kperp.ndim == 2:
#				return np.array([interp_kperp_ints(kpi, pi) for kpi in kperp]).squeeze() # this will loop over the output dimension = rp
#			else:
#				return interp_kperp_ints(kperp, pi)
#
#		for j, p in enumerate(Pi):
#			PerpInt = lambda kpp: LookUp(kpp, p)
#			#import pdb ; pdb.set_trace()
#			xi[z, j] += prefactor * ht.transform(PerpInt, rp, ret_err=False)
#			#print('%i / %i' % ((len(Pi)*z + j + 1), (len(pkzz) * len(Pi))))
#
#	#xi_tot = (xi.T * Wz).T.sum(axis=0) * np.diff(zbase)[0]
#	#wgp = xi_tot.sum(axis=0) * np.diff(Pi)[0]
#	xi_tot = np.trapz((xi.T * Wz), x=zbase, axis=2).T
#	wgp = np.trapz(xi_tot, x=Pi, axis=0)
###################################################
#### DIVISION BY 1+xi_gg BEFORE SUMMING OVER PI ###
###################################################
#
#	if plots[0]:
#		f, ax = plt.subplots(2, sharex=True)
#		plt.sca(ax[0])
#		l, = plt.loglog(rp, wgp, label='my wg+ out')
#		plt.loglog(rp, abs(wgp), c=l.get_color(), ls='--')
#		plt.plot(th, true_wgp, ls=':', label='cosmosis v17-fit out')
#		plt.plot(th, m06_wgp, ls=':', label='cosmosis m06-fit out')
#		for a in ax:
#			a.axvline(0.1,c='k',lw=0.7)
#			a.axvline(60,c='k',lw=0.7)
#		plt.ylabel('$w_{g+}$')
#		#plt.title('h=%f'%h, fontsize=35)
#		plt.xlim(1e-4,1e5)
#		plt.ylim(1e-4,1e1)
#		plt.legend(fontsize=14)
#		fx = interp1d(th, true_wgp)
#		plt.sca(ax[1])
#		plt.plot(rp, wgp/fx(rp) - 1., 'C0-')
#		plt.axhline(0, c='C1', ls=':')
#		plt.axhspan(-0.1, 0.1, color='grey', alpha=0.1)
#		plt.axhspan(-0.05, 0.05, color='grey', alpha=0.1)
#		plt.axhspan(-0.01, 0.01, color='grey', alpha=0.1)
#		plt.ylabel('mine/cosmosis - 1.')
#		plt.xlabel('$r_{p}$')
#		plt.tight_layout()
#		plt.show()
#
#	if plots[1]:
#		plt.figure()
#		pcm1 = plt.imshow(np.log10(xi_tot), cmap='Greens', extent=[np.log10(rp.min()),np.log10(rp.max()),Pi.min(),Pi.max()], aspect='auto')#, vmin=-6, interpolation='gaussian')
#		plt.colorbar(pcm1, extend='both')
#		if any(xi_tot.flatten() < 0):
#			pcm2 = plt.imshow(np.log10(-xi_tot), cmap='Greys', extent=[np.log10(rp.min()),np.log10(rp.max()),Pi.min(),Pi.max()], aspect='auto')#, vmin=-6, interpolation='gaussian')
#			plt.colorbar(pcm2, extend='both')
#		plt.contour(np.log10(np.abs(xi_tot)), cmap='Purples_r', origin='lower', extent=[np.log10(rp.min()),np.log10(rp.max()),Pi.min(),Pi.max()], aspect='auto')
#		plt.ylabel('$\Pi$', rotation='horizontal') ; plt.xlabel(r'${\rm{log}}(r_{p})$') ; plt.xlim(-1, 2) ; plt.tight_layout()
#		plt.show()
#	print("stepsize h = %f"%h)
#	print("max k_para = %f"%kpar_max)
#	print("%% diff cosmosis vs. mine:")
#	f = interp1d(th, true_wgp)
#	print([float("%.1f"%i) for i in 100*(f(rp) / wgp - 1.)])
#
#	return xi_tot
#
#if HTEST:
#	#hvals = np.linspace(0.07, 1, 5)
#	hvals = np.logspace(-3, -1, 5)
#	for h in hvals:
#		ht = hankel.HankelTransform(nu=2, h=h)
#		run()
#		plt.title('h=%f'%h, fontsize=35)
#		plt.tight_layout()
#elif KPTEST:
#	xi_tots = []
#	kpm_vals = np.logspace(-2, 0, 7)
#	for kpm in kpm_vals:
#		xt = run(kpm)
#		plt.title('max $k_{\parallel}$=%f'%kpm, fontsize=35)
#		plt.tight_layout()
#		xi_tots.append(xt)
#else:
#	xi_tot = run()
#
#
#
#
#
