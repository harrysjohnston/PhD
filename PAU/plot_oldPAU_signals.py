# coding: utf-8
from functions import *
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10

ks = ['pau', 'odds0.4_pau', 'red_pau', 'red_odds0.4_pau', 'blue_pau', 'blue_odds0.4_pau']
cols = dict(zip(ks, ['k','k','r','r','b','b']))
fmt = dict(zip(ks, ['-', '*--']*3))
lw = dict(zip(ks, [1.3, 1]*3))
aid = dict(zip(ks, [0, 0, 1, 1, 2, 2]))
labs = dict(zip(ks, ['All PAU', 'PAU odds>perc50', 'Red PAU', 'Red PAU odds>perc50', 'Blue PAU', 'Blue PAU odds>perc50']))
split = dict(zip(ks, [0.94, 1.06]*3))
kw1 = {'capsize':2, 'elinewidth':1.6, 'capthick':1.6}

err_ratio = 0

dd = read_files('dat', asc=1)
f, ax = plt.subplots(2, 3, sharex=True, sharey='row', figsize=(14,6))
for k in np.sort(dd.keys()):
	key = k.replace('wgg_','').replace('wgp_','').replace('.dat','')
	if 'wgg' in k:
		a = ax[0, :]
		rpp = 0
		wp = dd[k]['wgg']
	if 'wgp' in k:
		a = ax[1, :]
		rpp = 0.8
		wp = dd[k]['wgplus']
		wx = dd[k]['wgcross']
	try:
		e = dd[k]['jackknife_err']
	except KeyError:
		print "no jackknife error for %s"%k
		e = dd[k]['noise']
	r = dd[k]['rnom']
	nans = ~np.isnan(wp) & ~np.isnan(e)
	r, wp, e = map(lambda x: x[nans], [r, wp, e])

#	if 'odds' in k:
#		jkdat = read_files(k.replace('.dat','.jk'), d=0)
#		for jk in jkdat:
#			jkr, jkwp = jk.T[0], jk.T[3]
#			a.plot(jkr, jkr**rpp*jkwp, '%s--'%kw['c'], lw=0.1)

	a = a[aid[key]]
	s = split[key]
	lab = labs[key]
	if not err_ratio:
		kw = kw1.copy()
		kw['fmt'] = fmt[key]
	else:
		kw = {}
		if len(fmt[key])>1:
			kw['ls'] = fmt[key][1:]
		else:
			kw['ls'] = fmt[key]

	kw['c'] = cols[key]
	kw['lw'] = lw[key]

	if not err_ratio:
		err = a.errorbar(r*s, r**rpp*wp, r**rpp*e, alpha=0.6, label=lab, **kw)
		if 'odds' in k: err[-1][0].set_linestyle('--')
		if 'wgp' in k:
			kw['c'] = 'gray'
			nans &= ~np.isnan(wx)
			wx = wx[nans]
			errx = a.errorbar(r*s*0.85, r**rpp*wx, r**rpp*e, alpha=0.35, zorder=0, **kw)
			if 'odds' in k: errx[-1][0].set_linestyle('--')
	else:
		a.plot(r*s, e / dd[k]['noise'][nans], label=lab, **kw)

for i, a in enumerate(ax.flatten()):
	a.set_xscale('log')
	#a.set_yscale('log')
	if i > 2:
		a.axhline(0, ls=':', c='k', zorder=10)
	a.set_xlim(0.05,10)
	#a.set_ylim(0.01,None)
ax[0,0].set_ylabel(r'$w_{\rm{gg}}$')
ax[1,0].set_ylabel(r'$r_{p}^{0.8}w_{\rm{g+}}$')
for i in range(3):
	#ax[0,i].set_ylim(1,None)
	ax[0,i].set_yscale('log')
	ax[0,i].legend(loc='best')
	ax[1,i].set_xlabel(r'$r_{p} \, [h^{-1}\rm{Mpc}]$')
plt.show()
plt.tight_layout()



