# coding: utf-8
import matplotlib
matplotlib.use('Agg')
from functions import *
gcat = fopen('../SMLambdarApMatchedPhotom.fits')
gcat = gcat[gcat['PETROMAG_R']<=19.8]
pcat = fopen('PAUS_KSB.fits')
pcat = pcat[(pcat['bcnz_zb'] > 0.1) & (pcat['bcnz_zb'] < 0.8)]
gz = gcat['Z_TONRY']
pz = pcat['bcnz_zb']
gr = gcat['red_sequence']
pr = pcat['red_sequence_LePhare']
pb = pcat['blue_cloud_Cigale_2cluster']
p_absi = pcat['lp_mi']
g_absi = gcat['absmag_i']
p_appi = pcat['mag_i']
p_appr = pcat['mag_r']
#g_appi = np.where((g_absi>0)|(g_absi<-26),np.nan,g_absi) + MICEcosmo.distmod(gz).value
g_appi = gcat['PETROMAG_I']
selec1 = p_appr < 19.8
selec2 = pcat['Qz'] < p50(pcat['Qz'])

kw = dict(density=1,alpha=0.8,histtype='step',range=(14.5,23), bins=50, lw=1.5)
kwf = kw.copy()
kw1 = kw.copy()
kwf['histtype'] = 'stepfilled'
kw1f = kwf.copy()
kw1['range'] = kw1f['range'] = (-25, -16)
kwf['alpha'] = kw1f['alpha'] = 0.2
f, ax = plt.subplots(2, figsize=(7,9))
ax[0].hist(g_appi[gr], color='r', label='GAMA red', **kwf)
ax[0].hist(g_appi[~gr], color='b', label='GAMA blue', **kwf)
ax[0].hist(p_appi[pr], color='r', ls='-', label='PAUS W3 red', **kw)
ax[0].hist(p_appi[pb], color='b', ls='-', label='PAUS W3 blue', **kw)
ax[0].hist(p_appi[pr&selec2], color='r', ls='--', label=r'W3 red, Qz$_{50}$', **kw)
ax[0].hist(p_appi[pb&selec2], color='b', ls='--', label=r'W3 blue, Qz$_{50}$', **kw)
ax[0].hist(p_appi[pr&selec1], color='r', ls=':', label='W3 red, $r<19.8$', **kw)
ax[0].hist(p_appi[pb&selec1], color='b', ls=':', label='W3 blue, $r<19.8$', **kw)
ax[0].set_xlabel(r'apparent $i$ magnitude')
ax[0].set_ylabel(r'PDF')

ax[1].hist(g_absi[gr], color='r', **kw1f)
ax[1].hist(g_absi[~gr], color='b', **kw1f)
ax[1].hist(p_absi[pr], color='r', ls='-', **kw1)
ax[1].hist(p_absi[pb], color='b', ls='-', **kw1)
ax[1].hist(p_absi[pr&selec2], color='r', ls='--', **kw1)
ax[1].hist(p_absi[pb&selec2], color='b', ls='--', **kw1)
ax[1].hist(p_absi[pr&selec1], color='r', ls=':', **kw1)
ax[1].hist(p_absi[pb&selec1], color='b', ls=':', **kw1)
ax[1].set_xlabel(r'absolute $i$ magnitude')
ax[1].set_ylabel(r'PDF')

ax[0].legend(loc='upper left', fontsize=14, frameon=0)
plt.tight_layout()
plt.savefig('GAMA_PAU_samplecomp.png',bbox_inches='tight')
plt.savefig('GAMA_PAU_samplecomp.pdf',bbox_inches='tight')
plt.show()





