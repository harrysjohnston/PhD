# coding: utf-8
import matplotlib
matplotlib.use('Agg')
from functions import *

cpath = sys.argv[1]
zmax_path = sys.argv[2]
zsp_col = sys.argv[3]
zph_col = sys.argv[4]
idcol = sys.argv[5]
Ndraws = int(sys.argv[6])
zlimit = float(sys.argv[7])

cat = fopen(cpath)
cat = cat[np.argsort(cat[idcol])]
zmaxtable_h5 = h5py.File(zmax_path, 'r')

zmaxtable = zmaxtable_h5['zmax'][:]
zspectable = zmaxtable_h5['zspec'][:]
zmax_id = zmaxtable_h5['ID'][:]
true_zmax = zmaxtable_h5['true_zmax'][:][:,1]
if zmaxtable.shape[0] < zmaxtable.shape[1]:
	zmaxtable = zmaxtable.T
if zspectable.shape[0] < zspectable.shape[1]:
	zspectable = zspectable.T

cat = cat[np.isin(cat[idcol], zmax_id)]
have_zs = np.isfinite(cat[zsp_col]) & (cat[zsp_col] > 0)

#indices = np.random.choice(max(zmaxtable.shape), replace=False, size=Ndraws)
indices = range(len(cat))[::1000][:Ndraws]
if not all(have_zs[indices]):
	indices = np.argwhere(have_zs)[::40][:Ndraws].squeeze()

f, ax = plt.subplots(int(Ndraws**0.5),int(Ndraws**0.5),sharex=True,sharey=False,figsize=(16,9))
#f.text(0.5, 0.01, zmax_path, va='center', ha='center')
ai = iter(ax.flatten())

for x in indices:
    zs, zm, zmax = zspectable[x], zmaxtable[x], true_zmax[x]
    z, zp = cat[zsp_col][x], cat[zph_col][x]
    plt.sca(next(ai))
    plt.axvline(z, c='k', label=r'$z_{\rm{spec.}}$')
    plt.axvline(zp, c='r', label=r'$z_{\rm{phot.}}$')
    plt.axvline(zmax, c='b', ls=':', label=r'true $z_{\rm{max}}$')
    plt.hist(zs, color='k', ls='--', histtype='step', density=1, bins='auto', label=r'$z_{\rm{spec.}}$ draw')
    plt.hist(zm, color='b', ls=':', histtype='step', density=1, bins='auto', label=r'$z_{\rm{max}}$ distribution')

for a in ax.flatten():
	a.tick_params(labelleft=0)#, labelbottom=0)
	if a in ax[:,0]:
		a.set_ylabel('PDF')
	if a in ax[-1,:]:
		a.tick_params(labelbottom=1)
		a.set_xlabel('$z$')

ax.flatten()[0].legend(fontsize=14, loc='best', frameon=0)
plt.xlim(-0.02, zlimit)
plt.tight_layout()
plt.subplots_adjust(wspace=0.03,hspace=0.05)
plt.show()
plt.savefig(zmax_path.replace('.zmaxtable','.pdf'), bbox_inches='tight')
plt.savefig(zmax_path.replace('.zmaxtable','.png'), bbox_inches='tight')


