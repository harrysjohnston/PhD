# coding: utf-8
import matplotlib
matplotlib.use('Agg')
from functions import *
cat = fopen('PAUS_KSB.fits')
zp = cat['bcnz_zb']
zs = cat['ZBEST']
reds = cat['red_sequence_Cigale_2cluster']
c = (zs > 0) & (zp > 0)
col = cat['u_abs'] - cat['i_abs']
c1 = c & (cat['qz'] < np.percentile(cat['qz'][c],50))
c2 = c & (cat['qz'] >= np.percentile(cat['qz'][c],50))
delta_norm = (zp-zs)/(1+zs)

f, ax = plt.subplots()
hb = plt.hexbin(col[c1], delta_norm[c1],
    alpha=0.6,
	extent=[0,4.5,-0.05,0.05],
	gridsize=100,
	linewidths=0,
    #C=reds[c1],
	bins='log',
	cmap='BuPu', zorder=1)

delta_red = delta_norm[c1 & reds]
delta_blue = delta_norm[c1 & ~reds]
red_68ci = np.percentile(delta_red, [15.9, 84.1])
blue_68ci = np.percentile(delta_blue, [15.9, 84.1])
for pi in red_68ci:
	plt.axhline(pi, c='C3', ls='-', lw=2, alpha=0.8)
for pi in blue_68ci:
	plt.axhline(pi, c='C0', ls='-', lw=2, alpha=0.8)

plt.annotate(r'$\sigma_{68}\sim%.4f$'%(np.diff(red_68ci)[0]/2.), xy=(0.75,0.8), xycoords='axes fraction', color='C3', fontsize=12)
plt.annotate(r'$\sigma_{68}\sim%.4f$'%(np.diff(blue_68ci)[0]/2.), xy=(0.05,0.8), xycoords='axes fraction', color='C0', fontsize=12)

plt.ylim(-0.05, 0.05)
plt.axhline(0, c='k', ls='--', lw=1.3, zorder=10)
plt.xlabel(r'$u-i$')
plt.ylabel(r'$(z_{\rm{phot.}}-z_{\rm{spec.}}) / (1 + z_{\rm{spec.}})$', fontsize=15)
plt.tight_layout()
plt.show()

plt.savefig('pau_photoz_vs_colour.pdf', bbox_inches='tight')
plt.savefig('pau_photoz_vs_colour.png', bbox_inches='tight')



