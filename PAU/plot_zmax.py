# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
import sys
from astropy.io import fits

if len(sys.argv) < 2:
	print '\nGive catalogue path (column name hard-coded)\n'
	sys.exit(1)

cat = fits.open(sys.argv[1])[1].data

plt.figure()

zmax = cat['PETROMAG_R_fl19.8_zmax']
zmax_gama = cat['zmax_19p8']
ccoln = 'logage'
ccol = cat[ccoln]
p10, p90 = np.percentile(ccol, [10., 90.])
c = plt.scatter(zmax_gama, zmax, vmin=p10, vmax=p90,
            c=ccol, s=1, alpha=0.3)
cbar = plt.colorbar(c, extend='both')
cbar.ax.set_ylabel(ccoln)
plt.xlim(-0.05, 0.7)
plt.ylim(-0.05, 0.7)
plt.plot([0,0.5],[0,0.5],'r-')
plt.xlabel('GAMA zmax')
plt.ylabel('my zmax')
plt.show()
plt.tight_layout()

