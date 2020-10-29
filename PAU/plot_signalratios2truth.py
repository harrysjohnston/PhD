# coding: utf-8
from functions import *

if len(sys.argv) < 2:
	print '\
Must give:\n\
1: mode [old, red, blue, total]\n\
2+: directories containing outputs\n'
	sys.exit()

plot_errors = 0
mode = sys.argv[1]
fdirs = sys.argv[2:]
fdirs.sort()

kstr = '_Dzphot_2Rzphot_2.dat'
if mode == 'total' or mode == 'old':
	wgp_dat = ascii.read('OUTPUTS_specphot_test_unifPiBins/wgp_DspecRspec.dat')
	wgg_dat = ascii.read('OUTPUTS_specphot_test_unifPiBins/wgg_DspecRspec.dat')
elif mode == 'red':
	wgp_dat = ascii.read('OUTPUTS_specphot_test_unifPiBins/wgp_DspecRspec_red.dat')
	wgg_dat = ascii.read('OUTPUTS_specphot_test_unifPiBins/wgg_DspecRspec_red.dat')
elif mode == 'blue':
	wgp_dat = ascii.read('OUTPUTS_specphot_test_unifPiBins/wgp_DspecRspec_blue.dat')
	wgg_dat = ascii.read('OUTPUTS_specphot_test_unifPiBins/wgg_DspecRspec_blue.dat')
elif mode == 'mixed':
	wgp_dat = ascii.read('OUTPUTS_specphot_test_unifPiBins/wgp_DspecRspec_red.dat')
	wgg_dat = ascii.read('OUTPUTS_specphot_test_unifPiBins/wgg_DspecRspec.dat')

if mode == 'blue' or mode == 'mixed':
	f, ax = plt.subplots(2, sharex=True, figsize=(7.8,9))
else:
	f, ax = plt.subplots(2, sharex=True, figsize=(7.3,9))
r_ref = wgp_dat['rnom']
fac = r_ref ** 0.8
wgp_ref = wgp_dat['wgplus']
wgp_err_ref = wgp_dat['wgplus_jackknife_err']
wgg_ref = wgg_dat['wgg']
wgg_err_ref = wgg_dat['wgg_jackknife_err']

spec_col = {'old':'k','total':'k','mixed':'k','red':'r','blue':'b'}
if mode == 'mixed':
	ax[0].fill_between(r_ref, (1+wgp_err_ref/wgp_ref), (1-wgp_err_ref/wgp_ref), alpha=0.2, lw=0, color='r', label='spectroscopic-$z$\njackknife $1\sigma$', zorder=0)
	ax[0].axhline(1, c='r', ls='--')
else:
	ax[0].fill_between(r_ref, fac*(wgp_ref+wgp_err_ref), fac*(wgp_ref-wgp_err_ref), alpha=0.2, lw=0, color=spec_col[mode], zorder=0)
	ax[0].plot(r_ref, fac*wgp_ref, ls='--', c=spec_col[mode], zorder=0)
ax[1].fill_between(r_ref, (1+wgg_err_ref/wgg_ref), (1-wgg_err_ref/wgg_ref), alpha=0.2, lw=0, color=spec_col[mode], label='spectroscopic-$z$\njackknife $1\sigma$', zorder=0)
ax[0].axhline(0, c=spec_col[mode], ls=':')
ax[1].axhline(1, c=spec_col[mode], ls='--')

for fdir in fdirs:
	label = {
		'OUTPUTS_specphot_test_unifPiBins':'uniform $\Pi$-bins\nunwindowed randoms',
		'OUTPUTS_specphot_test_unifPiBins_windowed':'uniform $\Pi$-bins\nwindowed randoms',
		'OUTPUTS_specphot_test_unifPiBins_zph':'uniform $\Pi$-bins\nunwindowed zph-randoms',
		'OUTPUTS_specphot_test_unifPiBins_windowed_zph':'uniform $\Pi$-bins\nwindowed zph-randoms',
		'OUTPUTS_specphot_test_dynPiBins':'dynamic $\Pi$-bins\nunwindowed randoms',
		'OUTPUTS_specphot_test_dynPiBins_windowed':'dynamic $\Pi$-bins\nwindowed randoms',
		'OUTPUTS_specphot_test_dynPiBins_zph':'dynamic $\Pi$-bins\nunwindowed zph-randoms',
		'OUTPUTS_specphot_test_dynPiBins_windowed_zph':'dynamic $\Pi$-bins\nwindowed zph-randoms',
		'OUTPUTS_GAMA_zphot_2_uniform_windowed':'uniform $\Pi$-bins\nwindowed randoms',
		'OUTPUTS_GAMA_zphot_2_fibonacci_windowed':'Fibonacci $\Pi$-bins\nwindowed randoms',
		'OUTPUTS_GAMA_zphot_2_dynamic_windowed':'dynamic $\Pi$-bins\nwindowed randoms',
		'OUTPUTS_GAMA_zphot_2_uniform_unwindowed':'uniform $\Pi$-bins\nunwindowed randoms',
		'OUTPUTS_GAMA_zphot_2_fibonacci_unwindowed':'Fibonacci $\Pi$-bins\nunwindowed randoms',
		'OUTPUTS_GAMA_zphot_2_dynamic_unwindowed':'dynamic $\Pi$-bins\nunwindowed randoms',
		'OUTPUTS_GAMA_zphot_2_uniform_zph-windowed':'uniform $\Pi$-bins\nzph-windowed randoms',
		'OUTPUTS_GAMA_zphot_2_fibonacci_zph-windowed':'Fibonacci $\Pi$-bins\nzph-windowed randoms',
		'OUTPUTS_GAMA_zphot_2_dynamic_zph-windowed':'dynamic $\Pi$-bins\nzph-windowed randoms',
		'OUTPUTS_GAMA_zphot_2_uniform_zph-unwindowed':'uniform $\Pi$-bins\nzph-unwindowed randoms',
		'OUTPUTS_GAMA_zphot_2_fibonacci_zph-unwindowed':'Fibonacci $\Pi$-bins\nzph-unwindowed randoms',
		'OUTPUTS_GAMA_zphot_2_dynamic_zph-unwindowed':'dynamic $\Pi$-bins\nzph-unwindowed randoms',
		'OUTPUTS_GAMA_zphot_2_fibonacci_gama-windowed':'Fibonacci $\Pi$-bins\ngama-windowed randoms',
		'OUTPUTS_GAMA_Z_TONRY_uniform_gama-windowed':'spec-z uniform $\Pi$-bins\ngama-windowed randoms',
	}[fdir.strip('/')]
	#if mode == 'mixed':
	label = label.replace('Fibonacci', 'dynamic')
	if 'uniform' in label:
		fmt = 'v'
	if 'dynamic' in label:
		fmt = 'o'
	if 'Fibonacci' in label:
		fmt = 'd'
	if 'unw' in label:
		mfc = 'w'
	else:
		mfc = None
	if 'zph' in label:
		fmt += '-.'
	else:
		fmt += '-'
	if mode == 'old' or mode =='total':
		wgp_dati = ascii.read(join(fdir, 'wgp_total.dat'))
		wgg_dati = ascii.read(join(fdir, 'wgg_total.dat'))
	elif mode == 'red':
		wgp_dati = ascii.read(join(fdir, 'wgp_red.dat'))
		wgg_dati = ascii.read(join(fdir, 'wgg_red.dat'))
	elif mode == 'blue':
		wgp_dati = ascii.read(join(fdir, 'wgp_blue.dat'))
		wgg_dati = ascii.read(join(fdir, 'wgg_blue.dat'))
	elif mode == 'mixed':
		wgp_dati = ascii.read(join(fdir, 'wgp_red.dat'))
		wgg_dati = ascii.read(join(fdir, 'wgg_total.dat'))
	r = wgp_dati['rnom']
	fac = r ** 0.8
	wgp = wgp_dati['wgplus']
	wgg = wgg_dati['wgg']
	if 'wgplus_jackknife_err' in wgp_dati.keys() and plot_errors:
		err = wgp_dati['wgplus_jackknife_err']
		if mode != 'mixed':
			l1 = ax[0].errorbar(r, fac*wgp, fac*err, fmt=fmt, label=label, mfc=mfc, capsize=1.5)[0]
		else:
			l1 = ax[0].errorbar(r, wgp/wgp_ref, err/wgp_ref, fmt=fmt, label=label, mfc=mfc, capsize=1.5)[0]
	else:
		if mode != 'mixed':
			l1, = ax[0].plot(r, fac*wgp, fmt, label=label, mfc=mfc)
		else:
			l1, = ax[0].plot(r, wgp/wgp_ref, fmt, label=label, mfc=mfc)
	if 'wgg_jackknife_err' in wgg_dati.keys() and plot_errors:
		err = wgg_dati['wgg_jackknife_err']
		l2 = ax[1].errorbar(r, wgg/wgg_ref, err/wgg_ref, fmt=fmt, label=label, c=l1.get_c(), mfc=mfc, capsize=1.5)[0]
	else:
		l2, = ax[1].plot(r, wgg/wgg_ref, fmt, label=label, c=l1.get_c(), mfc=mfc)

for a in ax:
	a.set_xscale('log')
if mode == 'red' or mode == 'mixed':
	if mode == 'red': fs = 12
	else: fs = 10
	ax[1].legend(loc='best',ncol=2,fontsize=fs,frameon=0)
if mode == 'mixed':
	ax[0].set_ylim(None, 2.1)
else:
	ax[0].set_ylim(None, 1.5)
ax[1].set_ylim(0.3, 2)
if mode == 'blue':
	ax[1].set_ylabel(r'$w_{\rm{gg}} \, / \,  w_{\rm{gg,spec.}}$')
if mode == 'blue':
	ax[0].set_ylabel(r'$r_{p}^{0.8}w_{\rm{g+}}\,[h^{-1}\rm{Mpc}]^{1.8}$')
if mode == 'mixed':
	ax[1].set_ylabel(r'$w^{\rm{tot.}}_{\rm{gg}} \, / \,  w^{\rm{tot.}}_{\rm{gg,spec.}}$')
	ax[0].set_ylabel(r'$w^{\rm{red}}_{\rm{g+}} \, / \,  w^{\rm{red}}_{\rm{g+,spec.}}$')
ax[1].set_xlabel(r'$r_{p}\,[h^{-1}\rm{Mpc}]$')
plt.tight_layout()
plt.show()

plt.savefig('signalratios_gama_zphot_2_%s.png'%mode, bbox_inches='tight')
plt.savefig('signalratios_gama_zphot_2_%s.pdf'%mode, bbox_inches='tight')


