# coding: utf-8
import matplotlib
matplotlib.use('Agg')
from functions import *
gcat = fopen('../SMLambdarApMatchedPhotom.fits')
gcat = gcat[gcat['PETROMAG_R']<=19.8]
pcat = fopen('PAUS_KSB.fits')
gz = gcat['Z_TONRY']
pz = pcat['bcnz_zb']
gr = gcat['red_sequence']
ccol = 'red_sequence_LePhare'
pr = pcat[ccol]
kw = dict(normed=0,alpha=0.8,histtype='step')
plt.figure()
plt.hist(gz[gr], range=(0,0.7), bins=40, color='r', label='GAMA red', **kw)
plt.hist(gz[~gr], range=(0,0.7), bins=40, color='b', label='GAMA blue', **kw)
plt.hist(pz[pr], range=(0,1.2), bins=60, color='r', ls='--', label='PAUS W3 red (%s)'%ccol.replace('red_sequence_',''), **kw)
plt.hist(pz[~pr], range=(0,1.2), bins=60, color='b', ls='--', label='PAUS W3 blue (%s)'%ccol.replace('red_sequence_',''), **kw)
plt.xlabel(r'$z_{\rm{spec./phot.}}$')
plt.ylabel(r'$n(z_{\rm{spec./phot.}})$')
plt.legend(loc='best', fontsize=13, frameon=0)
plt.tight_layout()
plt.savefig('GAMA_PAU_nz.png', bbox_inches='tight')
plt.show()

