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
pr = pcat['red_sequence_Cigale_2cluster']
p_absi = pcat['lp_mi']
g_absi = gcat['absmag_i']
p_appi = pcat['mag_i']
p_appr = pcat['mag_r']
#g_appi = np.where((g_absi>0)|(g_absi<-26),np.nan,g_absi) + MICEcosmo.distmod(gz).value
g_appi = gcat['PETROMAG_I']
rlim = p_appr < 19.8

kw = dict(normed=1,alpha=0.8,histtype='step',range=(14.5,23), bins=50, lw=1.5)
kwf = kw.copy()
kw1 = kw.copy()
kwf['histtype'] = 'stepfilled'
kw1f = kwf.copy()
kw1['range'] = kw1f['range'] = (-25, -18)
kwf['alpha'] = kw1f['alpha'] = 0.2
f, ax = plt.subplots(2, figsize=(7,9))
ax[0].hist(g_appi[gr], color='r', label='GAMA red', **kwf)
ax[0].hist(g_appi[~gr], color='b', label='GAMA blue', **kwf)
ax[0].hist(p_appi[pr], color='r', ls='-', label='PAUS W3 red', **kw)
ax[0].hist(p_appi[~pr], color='b', ls='-', label='PAUS W3 blue', **kw)
ax[0].hist(p_appi[pr&rlim], color='r', ls='--', label='W3 red, $r<19.8$', **kw)
ax[0].hist(p_appi[~pr&rlim], color='b', ls='--', label='W3 blue, $r<19.8$', **kw)
ax[0].set_xlabel(r'apparent $i$ magnitude')
ax[0].set_ylabel(r'PDF')

ax[1].hist(g_absi[gr], color='r', **kw1f)
ax[1].hist(g_absi[~gr], color='b', **kw1f)
ax[1].hist(p_absi[pr], color='r', ls='-', **kw1)
ax[1].hist(p_absi[~pr], color='b', ls='-', **kw1)
ax[1].hist(p_absi[pr&rlim], color='r', ls='--', **kw1)
ax[1].hist(p_absi[~pr&rlim], color='b', ls='--', **kw1)
ax[1].set_xlabel(r'absolute $i$ magnitude')
ax[1].set_ylabel(r'PDF')

ax[0].legend(loc='upper left', fontsize=12)
plt.tight_layout()
plt.savefig('GAMA_PAU_samplecomp.png',bbox_inches='tight')
plt.show()





