# coding: utf-8
from functions import *
from clone_randoms import get_k_z
from numpy import log10
zg = np.linspace(0., 1., 200)
dg = cosmo.comoving_distance(zg)
dm_interp = interp1d(zg, dg, bounds_error=0, fill_value='extrapolate')

kc_gama = fopen('../GAMA_Kcorrs/kcorr_auto_z00v05.fits')
kc_me = pd.read_csv('../SMLambdarApMatchedPhotom.kcorrs', delimiter=' ')
uID = np.unique(kc_me['ID'])
kz, zg, m = get_k_z('../SMLambdarApMatchedPhotom.kcorrs')
#kz += 0.1
#zg = np.unique(kc_me['z'])
#t_me = np.array(kc_me['maggy_ratio']).reshape(len(uID), len(zg))
#idsort = np.argsort(kc_me['ID'][:len(uID)])

#t_me = t_me[np.isin(uID, kc_gama['CATAID'])]
#t_me = t_me[idsort]
t_ga = kc_gama[np.isin(kc_gama['CATAID'], uID)]
t_ga = t_ga[np.argsort(t_ga['CATAID'])]
assert all(uID == t_ga['CATAID']), "ID mismatch"

z = t_ga['Z_TONRY']
k_ga = t_ga['KCORR_R']
k_me = np.array([np.interp(z[i], zg, kz[i]) for i in range(len(kz))])
plt.figure()
plt.scatter(k_me, k_ga,
			s=1, alpha=0.1)
plt.plot([-1,2],[-1,2],'r-')
plt.axvline(0., c='k', ls='--')
plt.axhline(0., c='k', ls='--')
plt.xlabel('my k')
plt.ylabel('GAMA k')
plt.tight_layout()
plt.show()

#for redshift in [0.0, 0.1, 0.2]:
#	s = str(redshift).replace('.','')
#	z_gama = kc_gama['Z_TONRY']
#	zerocorr = -2.5 * log10(kc['maggy_ratio'][kc['z']==redshift])
#	zerocorr = np.array(zerocorr)
#	ids = kc['ID'][kc['z']==redshift]
#	sort = np.array(np.argsort(ids))
#	gamacorr = kc_gama['KCORR_R']
#	mask = np.isin(kc_gama['CATAID'], ids)
#	gsort = np.argsort(kc_gama['CATAID'][mask])
#	plt.figure()
#	dm = (5. * log10(dm_interp(kc_gama['Z_TONRY'])[mask][gsort]) - 5.)
#	dm -= 5. * log10(dm_interp(redshift)) - 5.
#	x = dm + (gamacorr[mask][gsort] - zerocorr[sort])
#	#c = plt.scatter(zerocorr[sort], gamacorr[mask][gsort],
#	c = plt.scatter(dm, x,
#					c=z_gama[mask][gsort], s=1, alpha=0.1)
#	cbar = plt.colorbar(c)
#	#plt.plot([-2,2],[-2,2],'r-')
#	#plt.xlim(-1,2.3)
#	#plt.ylim(-1,2.3)
#	#plt.xlabel('my kcorr')
#	#plt.ylabel('GAMA kcorr')
#	plt.xlabel('distmod')
#	plt.ylabel('m_obs - m_z')
#	plt.title('GAMA_KCORR z=%s'%redshift)
#	plt.tight_layout()

