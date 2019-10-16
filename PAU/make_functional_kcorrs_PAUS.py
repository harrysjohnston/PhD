# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii

#cat = fits.open('PAU_BCNZ_W3KSB.fits')[1].data
#gmr = cat['lp_mg'] - cat['lp_mr']
#idcol = 'numeric_id'
#outname = 'PAUS_functional.kcorrs'
cat = fits.open('MICE_forPAUSrandomstesting_dec20.fits')[1].data
gmr = cat['gr_gal']
idcol = 'unique_gal_id'
outname = 'MICE_functional.kcorrs'

gmr_bins = np.percentile(gmr, [0, 25, 50, 75, 100])
#gmr_bins[0] = gmr.min()
#gmr_bins[-1] = gmr.max()

plt.scatter(cat['sdss_i_abs_mag'], gmr, s=4,
            c=cat['z_cgal'], alpha=0.4)
plt.xlim(-26, -13)
plt.colorbar()
for b in gmr_bins:
    plt.axhline(b, c='r')
plt.show()

zgrid = np.linspace(0., 1.2, 121)
k1 = 1.2*zgrid**2.
k2 = 0.9*zgrid**3.
k3 = 0.5*zgrid**3. - 0.1
k4 = -0.25*np.ones_like(zgrid)
k1[zgrid<0.4] = zgrid[zgrid<0.4] * (k1[zgrid>0.4][0]/0.4)
k2[zgrid<0.4] = zgrid[zgrid<0.4] * (k2[zgrid>0.4][0]/0.4)
k3[zgrid<0.4] = zgrid[zgrid<0.4] * (k3[zgrid>0.4][0]/0.4)
k4[zgrid<0.4] = zgrid[zgrid<0.4] * (k4[zgrid>0.4][0]/0.4)
plt.figure()
for k in [k1,k2,k3,k4]:
    plt.plot(zgrid, k)
plt.show()

print 'making table..'
kc = None
for zg in zgrid:
	kcorrs = np.zeros(len(cat))
	for b in range(len(gmr_bins)-1):
		cut = (gmr > gmr_bins[b]) & (gmr <= gmr_bins[b+1])
		kofz = [k1,k2,k3,k4][b]
		kcorrs[cut] = kofz[zgrid==zg]
	if kc is None:
		kc = np.column_stack((cat[idcol], [zg]*len(cat), 10**(-kcorrs/2.5)))
	else:
		kc = np.concatenate((kc, np.column_stack((cat[idcol], [zg]*len(cat), 10**(-kcorrs/2.5)))))
print 'saving..'
ascii.write(kc, outname, names=['ID','z','maggy_ratio'], delimiter=' ', overwrite=1)
print 'saving k=0 table..'
kc1 = kc.copy()
kc1[:, 2] = 1.
ascii.write(kc1, 'MICE_zero.kcorrs', names=['ID','z','maggy_ratio'], delimiter=' ', overwrite=1)




