# coding: utf-8
import matplotlib
matplotlib.use('Agg')
from functions import *
cat = fopen('PAUS_KSB.fits')
f, ax = plt.subplots(2, sharex=True, figsize=(6,9))
plt.sca(ax[0])
hb = plt.hexbin(cat['ZBEST'], log10(cat['qz']),
    C=cat['red_sequence_Cigale_2cluster'],
    extent=[0,1.8,-1.6,2.3],
    gridsize=200,
    cmap='coolwarm')
lab = 'PAUS/DEEP2 W3 matched galaxies'
plt.title(lab)
#cb = plt.colorbar(hb)
plt.ylabel('log(Qz)')
#plt.xlabel(r'$z_{\rm{spec}}$')
plt.tight_layout()

plt.sca(ax[1])
hb = plt.scatter(cat['ZBEST'], cat['bcnz_zb'],
    c=cat['red_sequence_Cigale_2cluster'],
    s=0.2, alpha=1,
    cmap='coolwarm')
plt.plot([0,1.6],[0,1.6],'y--',lw=1.2,alpha=0.7)
plt.ylim(-0.02,1.2)
plt.xlim(-0.02,1.8)
plt.ylabel(r'$z_{\rm{phot.}}$')
plt.xlabel(r'$z_{\rm{spec.}}$')

plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.show()

plt.savefig('photoz_performance.png', bbox_inches='tight')

