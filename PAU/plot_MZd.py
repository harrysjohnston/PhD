# coding: utf-8
from functions import *
from clone_randoms import *
gama = fopen(PHD+'DEI_SMLv20_allbands_match.fits')
cat = fopen('PAU_BCNZ_W3KSB.fits')
f, ax = plt.subplots(2, sharex=True, sharey=True, figsize=(8,8))
for a in ax:
    plt.sca(a)
    plt.ylim(-12.5, -25)
    plt.xlim(0, 4000)
sargs = {'s':0.1,'alpha':1,'facecolors':'none','marker':'o'}
plt.sca(ax[0])
plt.scatter(dmax(19.8, gama['absmag_r'], ceiling=1900), gama['absmag_r'], s=0.01, zorder=10, c='k')
mr = gama['absmag_r']
chi = cosmo.comoving_distance(gama['Z_TONRY'])
g_scatt = plt.scatter(chi, mr, zorder=9, c=gama['fitphot_r_1'], vmin=15, vmax=20, **sargs)
g_cbar = plt.colorbar(g_scatt, extend='both')
plt.sca(ax[1])
mask = cat['lp_mi'] > -30
zcol = 'z_b'
ceil = (3770, 4300)[zcol == 'z_b']
plt.scatter(dmax(22.5, cat['lp_mi'], ceiling=ceil)[mask], cat['lp_mi'][mask], s=0.1, zorder=10, c='k')
pau_chi = cosmo.comoving_distance(cat[zcol])
p_scatt = plt.scatter(pau_chi[mask], cat['lp_mi'][mask], zorder=9, c=cat['mag_i'][mask], vmin=18, vmax=22.7, **sargs)
p_cbar = plt.colorbar(p_scatt, extend='both')
h = [new_handle(ls='-',lw=3,c='k'),
     new_handle(ls='',marker='o',markersize=3,c='C0')]
largs = {'fontsize':14, 'loc':'lower right'}
ax[0].legend(h[::-1], ['GAMA', r'$d_{\rm{max}}$ for app. r < 19.8'], **largs)
ax[1].legend(h[::-1], ['PAUS', r'$d_{\rm{max}}$ for app. i < 22.5'], **largs)
ax[0].set_ylabel('absolute R-mag')
ax[1].set_ylabel('absolute I-mag (LePhare)')
ax[1].set_xlabel('comoving distance [Mpc]')
g_cbar.ax.tick_params(labelsize=14)
p_cbar.ax.tick_params(labelsize=14)
g_cbar.ax.set_ylabel('app. R-mag')
p_cbar.ax.set_ylabel('app. I-mag')
plt.show()
plt.tight_layout()

