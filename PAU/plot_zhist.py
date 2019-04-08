# coding: utf-8
from functions import *
if len(sys.argv) < 6:
	print '\n1: catalogue\
		   \n2: randoms\
		   \n3: redshift column\
		   \n4: magnitude column\
		   \n5: ID column\n'
	sys.exit()

cat = fopen(sys.argv[1])
rand = fopen(sys.argv[2])
zcol = sys.argv[3]
magcol = sys.argv[4]
idcol = sys.argv[5]

hargs = {'alpha':0.2, 'normed':1}#, 'bins':'auto'}
cmask = (cat[magcol] > -30) & (cat[zcol] > 0.02)
rmask = rand[magcol+'_cloneZ'] > 0.02

sels = [(-19 < cat[magcol]) & (cat[magcol] < -18.9) & cmask,
		(-22 < cat[magcol]) & (cat[magcol] < -21.9) & cmask,
		(cat[magcol] > -17) & cmask]
labs = ['$-19 < M < -18.9$',
		'$-21 < M < -20.9$',
		'$-17 < M$']
#labs = ['M bin %s'%(i+1) for i in range(len(sels))]
#sels = [cmask]
#labs = ['all']
if cat[zcol].max() > 1.2:
	zmax = 1.25
else:
	zmax = cat[zcol].max()
nbin = np.linspace(0, zmax, 100)
f, ax = plt.subplots(figsize=(10,8))
plt.hist(cat[zcol][cmask], bins=nbin, **hargs)
plt.hist(rand[magcol+'_cloneZ'][rmask], bins=nbin, **hargs)

for i in range(len(sels)):
	scat = cat[sels[i]]
	idcut = np.isin(rand[magcol+'_cloneID'], scat[idcol])
	rcat = rand[idcut & rmask]

	lab = '%s: %s'%(i+1, labs[i])
	rlab = 'randoms %s'%(i+1)
	h = myhist(scat[zcol], bins=nbin,
				ls=':', lw=1.6, label=lab)
	myhist(rcat[magcol+'_cloneZ'], bins=nbin,
			lw=0.8, label=rlab, color=h[-1][0].get_edgecolor())

plt.xlim(-0.03, zmax)
plt.legend(fontsize=14, frameon=0)
plt.xlabel('$z$')
plt.ylabel('$P(z)$')
plt.show()
plt.tight_layout()


