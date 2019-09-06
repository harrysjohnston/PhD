# coding: utf-8
from functions import *
#run_startup()
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
if len(rand) > 1e7:
	rand = down(rand, 1e7/len(rand))
zcol = sys.argv[3]
magcol = sys.argv[4]
idcol = sys.argv[5]

#bg = np.loadtxt('badgalaxies.txt')

#hargs = {'alpha':0.2, 'normed':1}#, 'bins':'auto'}
cmask = ((cat['absmag_r'] > -26) &
		(cat['absmag_r'] < -17)	&
		(cat['PETROMAG_R'] < 19.8)	&
		(cat[zcol] > 0.02) &
		(cat[magcol+'_fl19.8_zmax'] > 0.))
		#~np.isin(cat[idcol], bg.T[0]))
rmask = rand[magcol+'_cloneZ'] > 0.0

zbin = np.linspace(0, 0.5, 6)
fmt = ['r-','g--','b-','m--','k-']
f, ax = plt.subplots()
for i in range(len(zbin)-1):
	scat = cat[cmask & (cat[zcol] > zbin[i]) & (cat[zcol] <= zbin[i+1])]
	sran = rand[np.isin(rand[magcol+'_cloneID'], scat[idcol])]
	h, b = np.histogram(sran[magcol+'_cloneZ'], bins=100)
	h = np.asarray(h, dtype=np.float32)
	h /= (h.sum() * np.diff(b)[0])
	plt.plot(midpoints(b), h, fmt[i], label='$%.1f<z<%.1f$'%(zbin[i], zbin[i+1]))

plt.xlim(0., 0.5)
plt.ylim(0., None)
plt.legend(fontsize=14, frameon=0)
plt.xlabel('$z$')
plt.ylabel('$P(z)$')
plt.show()
plt.tight_layout()


