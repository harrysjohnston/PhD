# coding: utf-8
import matplotlib
matplotlib.use('Agg')
from functions import *
catpath = 'K1000_SML.fits'
catpath = '../SMLambdarApMatchedPhotom.fits'
cat = fopen(catpath)
iband = cat['obsmag_i']
sigma = 0.02 * (1 + cat['Z_TONRY']) * iband / np.percentile(iband[(iband>14)&(iband<20.5)], 50.)
sigma_ = 0.02 * 1.3 * iband / np.percentile(iband[(iband>14)&(iband<20.5)], 50.)
outl_ = 0.15 * np.ones_like(iband) * iband / np.percentile(iband[(iband>14)&(iband<20.5)], 50.)
outl = outl_ + stats.norm.rvs(loc=0, scale=0.003, size=len(iband))
cut = (iband > 14) & (sigma > 0) & (outl > 0) & (iband < 20.5)

f, ax = plt.subplots(2, sharex=True, figsize=(7,9))
iband1, sigma1, sigma_1, outl1, outl_1 = iband[cut], sigma[cut], sigma_[cut], outl[cut], outl_[cut]
ax[0].plot(np.sort(iband1), sigma_1[np.argsort(iband1)], 'k--', lw=2, label='$z=0.3$')
#ax[1].plot(np.sort(iband1), outl_1[np.argsort(iband1)], 'k--', lw=2)
h1 = ax[0].hexbin(iband1, sigma1,
	#bins='log',
	C=cat['Z_TONRY'][cut],
	cmap='Spectral')
h2 = ax[1].hexbin(iband1, outl1,
	#bins='log',
	C=cat['Z_TONRY'][cut],
	cmap='Spectral')
cb1 = plt.colorbar(h1, ax=ax[0])
cb2 = plt.colorbar(h2, ax=ax[1])
#cb1.ax.set_ylabel('log(count)')
#cb2.ax.set_ylabel('log(count)')
cb1.ax.set_ylabel('spectroscopic redshift', fontsize=14)
cb2.ax.set_ylabel('spectroscopic redshift', fontsize=14)
ax[0].legend(loc='upper left', fontsize=14, frameon=0)
ax[0].set_ylabel(r'$\sigma_{z_{\rm{phot.}}}$')
ax[1].set_ylabel(r'$P_{\rm{kick}}$')
ax[1].set_xlabel('apparent $i$ magnitude')
plt.tight_layout()
plt.savefig('gama_outlscatt_test.png', bbox_inches='tight')

mask = cut
zphot_1 = np.ones_like(iband) * -99.
zphot_2 = np.ones_like(iband) * -99.
zphot_1[mask] = cat['Z_TONRY'][mask] + stats.norm.rvs(loc=0, scale=sigma[mask])
#zphot_2[mask] = cat['Z_TONRY'][mask] + stats.norm.rvs(loc=0, scale=sigma[mask])
zphot_2 = zphot_1.copy()
kick = 0.07 + np.random.rand(len(cat[mask]))/10.
kick *= np.sign(np.random.rand(len(cat[mask]))-0.6)
kick[np.random.rand(len(kick)) > outl[mask]] = 0.
zphot_2[mask] += kick
cat_outl1 = np.random.rand(len(zphot_2)) < 0.03
cat_outl2 = stats.norm.pdf(cat['Z_TONRY'], loc=0.35, scale=0.025)
cat_outl2 = 0.4 * cat_outl2 / cat_outl2.max()
cat_outl2 = np.random.rand(len(zphot_2)) < cat_outl2
zphot_2[cat_outl1 | cat_outl2] = 0.27 + stats.norm.rvs(loc=0, scale=0.025, size=(cat_outl1|cat_outl2).sum())
while any(zphot_1 < 0) or any(zphot_2 < 0):
    zphot_1[zphot_1 < 0] += 3. * (0. - zphot_1[zphot_1 < 0])
    zphot_2[zphot_2 < 0] += 3. * (0. - zphot_2[zphot_2 < 0])
zphot_1[zphot_1==99] = -99.
zphot_2[zphot_2==99] = -99.

plt.figure()
plt.hist(cat['Z_TONRY'][mask], bins=40, alpha=0.5, density=1, label='GAMA spec-$z$')
#myhist(zphot_1[mask], bins=40, lw=2, label='zphot_1: scattered')
myhist(zphot_2[mask], bins=40, lw=2, color='r', label='PAUS-like photo-$z$')
plt.xlabel('$z$')
plt.xlim(-0.02, 0.75)
plt.ylabel('PDF')
plt.legend(loc='upper right', fontsize=14, frameon=0)
plt.tight_layout()
plt.show()
plt.savefig('gama_specphot_nz_test.pdf', bbox_inches='tight')
plt.savefig('gama_specphot_nz_test.png', bbox_inches='tight')

t = Table(cat)
t['zphot_1'] = np.where(t['zphot_1'] > 0.65, -99, t['zphot_1'])
t['zphot_2'] = np.where(t['zphot_2'] > 0.65, -99, t['zphot_2'])
t['comoving_mpc_zphot_1'] = MICEcosmo.comoving_distance(t['zphot_1']) * MICEcosmo.h
t['comoving_mpc_zphot_2'] = MICEcosmo.comoving_distance(t['zphot_2']) * MICEcosmo.h
#t.write(catpath, overwrite=1)






