# coding: utf-8
from functions import *

f, ax = plt.subplots(2, 2, sharex=True, figsize=(14,9))
ax = ax.flatten()
ls = iter(['-','--',':','-.'])
spl = iter(splitN(len(sys.argv[1:]), 0.05))
for gdir in sys.argv[1:]:
	dd = read_files('.dat', dire=gdir, asc=1)
	lsi = next(ls)
	spli = next(spl)
	lab = gdir.replace('OUTPUTS_GAMA_','').replace('Z_TONRY','spec-z')

	axd = {'wgp_red.dat':ax[0], 'wgp_blue.dat':ax[1], 'wgg_red.dat':ax[2], 'wgg_blue.dat':ax[3]}
	rpd = {'wgp_red.dat':0.8, 'wgp_blue.dat':0.8, 'wgg_red.dat':0.6, 'wgg_blue.dat':0.4}
	cold = {'wgp_red.dat':'r', 'wgp_blue.dat':'b', 'wgg_red.dat':'r', 'wgg_blue.dat':'b'}
	for k in dd.keys():
		if k not in axd.keys():
			continue
		if 'wgp' in k:
			r, w = dd[k]['rnom'], dd[k]['wgplus']
			try: e = dd[k]['wgplus_jackknife_err']
			except: e = np.zeros_like(r)
		if 'wgg' in k:
			r, w = dd[k]['rnom'], dd[k]['wgg']
			try: e = dd[k]['wgg_jackknife_err']
			except: e = np.zeros_like(r)
		a = axd[k]
		rpp = rpd[k]
		col = cold[k]
		eb = a.errorbar(r*spli, r**rpp*w, r**rpp*e, fmt='%s%s'%(col,lsi), label=lab)
		eb[-1][0].set_linestyle(lsi)
		
for a in ax:
	a.axhline(0, lw=0.7, c='k')
ax[0].set_ylim(-0.5, 5)
ax[1].set_ylim(-1, 2)
#ax[2].set_ylim(-50, 600)
#ax[3].set_ylim(-50, 400)

ax[0].legend(loc='upper center', fontsize=10, frameon=0)
ax[0].set_xscale('log')
ax[0].set_ylabel('wg+,red')
ax[1].set_ylabel('wg+,blue')
ax[2].set_ylabel('wgg,red')
ax[3].set_ylabel('wgg,blue')
ax[3].set_xlabel('rp')
plt.tight_layout()
plt.show()

