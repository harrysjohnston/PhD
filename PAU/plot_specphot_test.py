# coding: utf-8
from functions import *
fdir = sys.argv[1]
flist = glob('%s/wg[gp]*dat'%fdir)
flist.sort()

wgp_ref = ascii.read('OUTPUTS_specphot_test_unifPiBins/wgp_DspecRspec.dat')
wgg_ref = ascii.read('OUTPUTS_specphot_test_unifPiBins/wgg_DspecRspec.dat')

f, ax = plt.subplots(2, sharex=True, figsize=(7,7))
ax[0].errorbar(wgp_ref['rnom'], wgp_ref['rnom']**0.8 * wgp_ref['wgplus'], wgp_ref['rnom']**0.8 * wgp_ref['wgplus_jackknife_err'], fmt='ro:', label='GAMA spec-z (truth)')
ax[1].errorbar(wgg_ref['rnom'], wgg_ref['wgg'], wgg_ref['wgg_jackknife_err'], fmt='ro:', label='GAMA spec-z (truth)')

s1 = iter(splitN(len(flist)/2))
s2 = iter(splitN(len(flist)/2))
for fl in flist:
	if ('bj' in fl or
		'pair' in fl or 'ic.dat' in fl or
		'Dzphot_1Rspec' in fl or
		'Rspec' in fl or
		1 != 1): continue
	if 'wgp' in fl:
		s = next(s1)
		a = ax[0]
		dat = ascii.read(fl)
		r, w = dat['rnom'], dat['wgplus']
		try:
			e = dat['wgplus_jackknife_err']
		except: e = np.zeros_like(w)
		fac = r**0.8
	elif 'wgg' in fl:
		s = next(s2)
		a = ax[1]
		if 'bj' not in fl:
			dat = ascii.read(fl)
			r, w = dat['rnom'], dat['wgg']
			try:
				e = dat['wgg_jackknife_err1']
			except: e = np.zeros_like(w)
		else:
			r, w, e = np.loadtxt(fl).T
		fac = 1.
	print 'not plotting erorrs'
	e *= 0
	eb = a.errorbar(r, fac*w, fac*e, fmt='.-', capsize=1, elinewidth=1.3, alpha=0.7,
				label=r'$%s$'%fl.replace('OUTPUTS_specphot_test', '').replace('wgg_','').replace('wgp_','').replace('.dat','').replace('_','').replace('zphot1','_{zphot1}').replace('zphot2','_{zphot2}').replace('spec','_{spec}'))
	if 'Rzphot_2' in fl:
		eb[0].set_c('C2')
for a in ax:
	a.axhline(0, ls=':', c='k')
	a.set_xscale('log')
	a.legend(loc='best')
ax[1].set_yscale('log')
ax[1].set_ylim(1, 1e3)
ax[0].set_ylabel(r'$r_{p}^{0.8}w_{\rm{g+}}\,[h^{-1}\rm{Mpc}]^{1.8}$')
ax[1].set_ylabel(r'$r_{p}w_{\rm{gg}}\,[h^{-1}\rm{Mpc}]^{2}$')
ax[1].set_xlabel(r'$r_{p}\,[h^{-1}\rm{Mpc}]$')
plt.tight_layout()
plt.show()


