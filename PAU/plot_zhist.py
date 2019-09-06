# coding: utf-8
from functions import *
run_startup()
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
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

#hargs = {'alpha':0.2, 'normed':1}#, 'bins':'auto'}
cmask = ((cat['absmag_r'] > -26) &
		(cat['absmag_r'] < -17)	&
		(cat[zcol] > 0.0) &
		(cat[magcol+'_fl19.8_zmax'] > 0.))
		#~np.isin(cat[idcol], bg.T[0]))
rmask = rand[magcol+'_cloneZ'] > 0.0


cutcol = 'logmstar'
Nbin = 3
percs = np.percentile(cat[cutcol], np.linspace(2., 98., Nbin+1))
sels = [(cat[cutcol] > percs[i]) & (cat[cutcol] < percs[i+1]) & cmask for i in range(Nbin)]
labs = ['%.1f < %s < %.1f' % (percs[i], cutcol, percs[i+1]) for i in range(Nbin)]
if cat[zcol].max() > 1.2:
	zmax = 1.25
else:
	zmax = cat[zcol].max()
nbin = 50
area = 180.
f, ax = plt.subplots(figsize=(8,6))
barplot = lambda x, **args:   plt.bar(midpoints(x[1])[:-1], np.diff(np.cumsum(1.*x[0]/x[0].sum())), width=np.diff(x[1])[0], **args)
stepplot = lambda x, **args: plt.plot(midpoints(x[1])[:-1], np.diff(np.cumsum(1.*x[0]/x[0].sum())), **args)

allcat_hist =  np.histogram(cat[zcol][cmask], bins=nbin*2)
allrand_hist = np.histogram(rand[magcol+'_cloneZ'][rmask], bins=nbin*2)
barplot(allcat_hist , alpha=0.2)
barplot(allrand_hist, alpha=0.2)
#stepplot(allcat_hist , c='grey', ls='-', alpha=0.4)
#stepplot(allrand_hist, c='grey', ls='--', lw=1.8, alpha=0.4)

for i in range(len(sels)):
	scat = cat[sels[i]]
	idcut = np.isin(rand[magcol+'_cloneID'], scat[idcol])
	rcat = rand[idcut & rmask]
	if len(rcat)>1e5:
		rcat = down(rcat, 1e5/len(rcat))

	lab = '%s: %s'%(i+1, labs[i])
	rlab = 'randoms %s'%(i+1)
	h = np.histogram(scat[zcol], bins=nbin)
	l = stepplot(h, ls='-',  label=lab)
	rh = np.histogram(rcat[magcol+'_cloneZ'], bins=nbin)
	stepplot(rh, ls='--', lw=1.8, label=rlab, color=l[0].get_color())

plt.xlim(-0.03, zmax)
plt.legend(fontsize=14, frameon=0)
plt.xlabel('$z$')
plt.ylabel('$P(z)$')
plt.show()
plt.tight_layout()


