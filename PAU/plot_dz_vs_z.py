# coding: utf-8
from functions import *
from matplotlib.colors import LogNorm
from astropy.cosmology import z_at_value
cmd = lambda x: cosmo.comoving_distance(x).value * cosmo.h

# redshifts and comoving distances
zgrid = np.linspace(0., 10., 20000)
dgrid = cmd(zgrid)
get_z = lambda d: interp1d(dgrid, zgrid)(d)
get_d = lambda z: interp1d(zgrid, dgrid)(z)
zmax = 1.2
zrange = np.linspace(0., zmax, 100)
drange = get_d(zrange)

# areas and volumes
ster = lambda x: 4.*np.pi * x / 41253.
fsky = lambda Omega: Omega / (4.*np.pi)
def invvol(O, r1, vol):
	# return the 2sigma upper/lower distance limits
	# of a fixed 1sigma volume, for a survey area O[ster],
	# centered at comoving distance r1
	r2 = ((3. * vol / O) + r1**3.)**(1./3.)
	width = r2 - r1
	new_r1 = r1 - width/2.
	new_r2 = ((3. * vol / O) + new_r1**3.)**(1./3.) # 1s limits
	#new_width = new_r2 - new_r1
	#new_r1 = new_r1 - new_width/2.
	#new_r2 = new_r2 + new_width/2. # 2s limits
	return [new_r1, new_r2]

sqdeg = np.array([180])
Omegas = ster(sqdeg)

Nr = 50
volumes = np.logspace(6, 7.5, Nr)

widths = np.array([[[invvol(O, d, v*2)
					for d in drange]
					for v in volumes]
					for O in Omegas]).squeeze()

# convert window width (depth along LoS) into a dz(z)
dmat = widths.copy()
dmat[dmat < 0] = 0.
dmat[dmat > dgrid.max()] = None
zg = get_z(dmat)
dzmat = np.diff(zg, axis=np.where(np.array(zg.shape) == 2)[0][0]).squeeze()
dzmat[dzmat > zmax] = None
dzmat[dzmat < 0] = None

f, ax = plt.subplots()
c = plt.pcolormesh(zrange, log10(volumes), dzmat, norm=LogNorm())
cbar = plt.colorbar()#, norm=LogNorm())
cbar.ax.minorticks_off()
C = plt.contour(zrange, log10(volumes), dzmat, colors='w',
				levels=np.percentile(dzmat.flatten(), np.arange(10)*10))
plt.clabel(C)
cbar.ax.set_ylabel('$\delta{z}$ ($\pm2\sigma$ window)')
plt.ylabel('$(1\sigma)$ log$_{10}$(Volume [Mpc/h]$^{3}$)')
plt.xlabel('$z$')
plt.axhline(log10(3.5e6), c='r', ls=':', label='GAMA')
plt.axvline(0.4, c='k', lw=1)
plt.axvline(0.2, c='k', lw=1)
plt.axvline(0.15, c='k', lw=1)
plt.legend(loc='lower right')
plt.show()
plt.tight_layout()





