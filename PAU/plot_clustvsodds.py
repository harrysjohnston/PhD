# coding: utf-8
from functions import *

#dirs = ['OUTPUTS_wgg_vs_qz_3d', 'OUTPUTS_wgg_qzweights_3d_uncomp']
#labs = ['wgg_qzltperc', 'wgg_wqz_qzltperc']
dirs = ['OUTPUTS_wgg_qzweights', 'OUTPUTS_wgg_qzweights_3d']
labs = ['wgg_wqz_qzltperc']*2
wdir = os.getcwd()
cmaps = [plt.cm.viridis_r, plt.cm.spring]

f, ax = plt.subplots(2, figsize=(8,10))
for j in range(len(dirs)):
	os.chdir(dirs[j])
	lab = labs[j]
	dat = read_files('.dat', asc=1)
	dkeys = list(np.sort(dat.keys()))
	try:
		dkeys.insert(len(dkeys), dkeys.pop(dkeys.index('%s100.dat'%lab)))
	except:
		del dat[lab+'0.dat'], dkeys[dkeys.index(lab+'0.dat')]
	cmap = cmaps[j](np.linspace(0.1, 0.9, len(dkeys)))

	for i, k in enumerate(dkeys):
		r = dat[k]['rnom']
		w = dat[k]['wgg']
		e = dat[k]['noise']
		try:
			e = dat[k]['jackknife_err']
			#jm = dat[k]['jackknife_mean']
			#plt.plot(r*s, jm, '*', c='gray', alpha=0.4, label=None)
		except:
			print k, '-- no jackknife'
		s = splitN(len(dkeys), gap=0.02)[i]
		plt.sca(ax[0])
		plt.errorbar(r*s, w, e, fmt='o-', c=cmap[i], alpha=0.7)
		plt.sca(ax[1])
		plt.errorbar(r*s, w/e, fmt='^-', c=cmap[i], alpha=0.7)
	os.chdir(wdir)
ax[0].set_title('Projected clustering, weighted by 1 / Qz')
h = [new_handle(marker='o', ls='-', c=cmap[i], label='Qz < p%s'%((i+1)*10)) for i in range(len(dkeys))]
ax[0].legend(handles=h, loc='best', fontsize=11, ncol=2)
ax[1].set_xlabel(r'$r_{p}\,[h^{-1}\rm{Mpc}]$')
ax[0].set_ylabel(r'$w_{\rm{gg}}\,[h^{-1}\rm{Mpc}]$')
ax[1].set_ylabel('S/N')
ax[0].set_ylim(0.5, 6e3)
ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[1].set_xscale('log')
plt.tight_layout()
plt.show()
