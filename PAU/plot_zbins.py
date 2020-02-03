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
cmask = np.ones(len(cat), dtype=bool)

zcol = sys.argv[2]
magcol = sys.argv[3]
idcol = sys.argv[4]
randlist = sys.argv[5:]

zbin = np.linspace(0, cat[zcol].max(), 6)
if 'SMLambdarApMatchedPhotom.fits' in sys.argv[1]:
	zbin = np.linspace(0, 0.5, 6)
	tag = 'GAMA'
else:
	tag = 'PAUS'

f, ax = plt.subplots(len(randlist), sharex=True, figsize=(6,len(randlist)*4.5))
axit = iter(ax)

fmt = ['r-','g--','b-','m--','k-']
for ranf in randlist:
	rand = fopen(ranf)
	if len(rand) > 1e7:
		rand = down(rand, 1e7/len(rand))

	rmask = (rand[magcol+'_cloneZ'] > 0.0) & (np.isin(rand[magcol+'_cloneID'], cat[cmask][idcol]))
	axi = next(axit)

	for i in range(len(zbin)-1):
		scat = cat[cmask & (cat[zcol] > zbin[i]) & (cat[zcol] <= zbin[i+1])]
		randcut = np.isin(rand[magcol+'_cloneID'], scat[idcol])
		print frac(randcut)
		sran = rand[randcut]

		h, b = np.histogram(sran[magcol+'_cloneZ'], range=(0, cat[zcol].max()), bins='auto')
		h = np.asarray(h, dtype=np.float32)
		h /= (h.sum() * np.diff(b)[0])
		axi.plot(midpoints(b), h, fmt[i], label=r'$%.1f<z_{\rm{gal}}<%.1f$'%(zbin[i], zbin[i+1]))

	axi.set_ylim(0., None)
	lab = ranf.replace(sys.argv[1].replace('.fits',''), '').replace('_CloneZIDRandoms_','').replace('.fits','')
	if 'Unwindowed' in lab:
		lab = 'unwindowed randoms'
	else:
		lab = 'windowed randoms'
	axi.set_title(lab, loc='right')

plt.xlim(0., cat[zcol].max())
ax[0].legend(fontsize=14, frameon=0)
ax[1].set_xlabel(r'$z_{\rm{rand}}$')
for a in ax:
	a.set_ylabel('PDF')
plt.show()
plt.tight_layout()

plt.savefig('%s_zbin_farrowplot.pdf'%tag, bbox_inches='tight')


