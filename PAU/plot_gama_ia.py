# coding: utf-8
from functions import *
xx = 'pz01'
gdat = read_files('highZ_vs', dire='../Wcorr_%sfullGAMA_rb_c0.66/to_plot/'%xx)
gdat['total'] = read_files('highZ_vs', dire='../Wcorr_%sfullGAMA_allz/to_plot/'%xx).values()[0]
for gdir in glob('../OUTPUTS_GAMA_f*PW1'):
	if 'pz' in gdir: continue
	dd = read_files('wgp*dat', dire=gdir, asc=1)
	f, (a1,a2) = plt.subplots(2, sharex=True, sharey=True, figsize=(8,8))
	a1.set_title(gdir)
	a1.set_xscale('log')
	a2.set_xscale('log')
	a1.axhline(0, c='k', lw=0.7)
	a2.axhline(0, c='k', lw=0.7)
	for k in np.sort(gdat.keys()):
		if 'Red' in k: c = 'r'
		elif 'Blue' in k: c = 'b'
		else: c = 'goldenrod'
		r, wp, wx, wpe, wxe = gdat[k].T[[0,1,4,3,6]]
		#wpe = wxe = np.zeros_like(r)
		l1 = a1.errorbar(r, r**0.8*wp, r**0.8*wpe, fmt='o:', c=c, ms=8, mfc='none', alpha=0.6, label=basename(k).replace('.dat',''))
		a2.errorbar(r, r**0.8*wx, r**0.8*wxe, fmt='o:', c=l1[0].get_c(), ms=8, mfc='none', alpha=0.6)
	for k in np.sort(dd.keys()):
		try:
			r, wp, wx, wpe, wxe = dd[k]['rnom'], dd[k]['wgplus'], dd[k]['wgcross'], dd[k]['wgplus_jackknife_err'], dd[k]['wgcross_jackknife_err']
		except KeyError:
			r, wp, wx, wpe, wxe = dd[k]['rnom'], dd[k]['wgplus'], dd[k]['wgcross'], np.zeros_like(dd[k]['rnom']), np.zeros_like(dd[k]['rnom'])
		if 'red' in k: c = 'r'
		elif 'blue' in k: c = 'b'
		else: c = 'goldenrod'
		l1 = a1.errorbar(r, r**0.8*wp, r**0.8*wpe, fmt='+-', c=c, label=basename(k).replace('_unwindowed.dat',''))
		a2.errorbar(r, r**0.8*wx, r**0.8*wxe, fmt='x-', c=l1[0].get_c())
	a1.legend(ncol=2, loc='best')
	plt.tight_layout()
	#break
    
plt.show()
