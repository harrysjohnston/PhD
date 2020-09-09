# coding: utf-8
from functions import *
flist = glob('../OUTPUTS_GAMA*/wgg_*dat') + glob('OUTPUTS_PAUS_KSB*/wgg_*dat')
#flist = glob('OUTPUTS_PAUS_KSB*/wgg_*dat')
flist = np.sort(flist)
#cut = np.array(['../OUTPUTS_GAMA_f0/' in i for i in flist]) | np.array(['qzselec' in j for j in flist]) | np.array(['pz' in k for k in flist])
#flist = flist[~cut]

rpp = 0

f, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(13,5))
for i, datf in enumerate(flist):
	#if 'PAUS' in datf: continue
	dat = ascii.read(datf)
	try:
		r, w, e = dat['rnom'].data, dat['wgg'].data, dat['wgg_jackknife_err'].data
	except KeyError:
		r, w, e = dat['rnom'].data, dat['wgg'].data, np.zeros_like(dat['wgg'].data)
		print 'no err for ', datf
	if 'red' in datf: c = 'r'
	elif 'blue' in datf: c = 'b'
	else: c = 'c'
	#if 'pz' in datf: m = 's'
	if 'PAU' in datf: m = 's'
	elif 'GAMA' in datf: m = '.'
	else: m = '_-.'
	s = splitN(len(flist), 0.01)[i]
	#if 'pz' in datf: a = [ax[1]]
	if 'PAU' in datf: a = [ax[1]]
	elif 'GAMA' in datf: a = [ax[0]]
	else: a = ax
	for ai in a:
		ai.errorbar(r*s, r**rpp*w, r**rpp*e, fmt='%s%s'%(c, m), lw=0.7, ms=3,
				label=datf.replace('OUTPUTS_','').replace('_unwindowed.dat','').replace('_KSBtotal_fg1fg2','').replace('_f0','').replace('_clust','0.01'))
for a in ax:
	a.grid(which='both', lw=0.3, alpha=0.3, c='grey')
	a.set_xscale('log')
	if rpp == 0:
		a.set_yscale('log')
	a.legend(ncol=2)
	a.set_xlabel(r'$r_{p}\,[h^{-1}\rm{Mpc}]$')
ax[0].set_ylabel(r'$w_{\rm{gg}}(r_{p})\,[h^{-1}\rm{Mpc}]$')
#plt.axhline(0, c='k', lw=0.7)
plt.tight_layout()
plt.show()



