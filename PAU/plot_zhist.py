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

#bg = np.loadtxt('badgalaxies.txt')

hargs = {'alpha':0.2, 'normed':1}#, 'bins':'auto'}
cmask = ((cat['absmag_r'] > -26) &
		(cat['absmag_r'] < -16)	&
		(cat[zcol] > 0.02) &
		(cat[magcol+'_fl19.8_zmax'] > 0.))
		#~np.isin(cat[idcol], bg.T[0]))
rmask = rand[magcol+'_cloneZ'] > 0.02


cutcol = 'logmstar'
percs = np.percentile(cat[cutcol], [10., 30., 60., 80.])
sels = [(cat[cutcol] < percs[0]) & cmask,
		(cat[cutcol] > percs[1]) & (cat[cutcol] < percs[2]) & cmask,
		(cat[cutcol] > percs[3]) & cmask]
labs = ['%s < %.1f' % (cutcol, percs[0]),
		'%.1f < %s < %.1f' % (percs[1], cutcol, percs[2]),
		'%s > %.1f' % (cutcol, percs[3])]
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
plt.hist(rand[magcol+'_cloneZ'][rmask], bins='auto', **hargs)

for i in range(len(sels)):
	scat = cat[sels[i]]
	idcut = np.isin(rand[magcol+'_cloneID'], scat[idcol])
	rcat = rand[idcut & rmask]
	if len(rcat)>1e5:
		rcat = down(rcat, 1e5/len(rcat))

	lab = '%s: %s'%(i+1, labs[i])
	rlab = 'randoms %s'%(i+1)
	h = myhist(scat[zcol], bins=nbin,
				ls=':', lw=1.6, label=lab)
	myhist(rcat[magcol+'_cloneZ'], bins='auto',
			lw=0.8, label=rlab, color=h[-1][0].get_edgecolor())

plt.xlim(-0.03, zmax)
plt.legend(fontsize=14, frameon=0)
plt.xlabel('$z$')
plt.ylabel('$P(z)$')
plt.show()
plt.tight_layout()


