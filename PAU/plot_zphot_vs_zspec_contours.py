# coding: utf-8
import matplotlib
matplotlib.use('Agg')
from functions import *
from glob import glob
#from chainconsumer import ChainConsumer

if len(sys.argv) < 6:
	print '\ngive:\n\
1: catalog path\n\
2: z_spec column\n\
3: z_phot column\n\
4: detection band column\n\
5: ID column\n\
6+: randoms files to plot\n'
	sys.exit()

catpath = sys.argv[1]
if 'PAU' in catpath:
	name = 'PAUS'
	ll = 0.02
else:
	name = 'GAMA'
	ll = 0.
outfile = name+'_zszp_contours'
zspcol = sys.argv[2]
zphcol = sys.argv[3]
magcol = sys.argv[4]
idcol = sys.argv[5]
ranlist = np.sort(sys.argv[6:])

cat = fopen(catpath)
cat = cat[np.argsort(cat[idcol])]
cat_zsp = cat[zspcol]
cat_zph = cat[zphcol]
cat_phsp = np.column_stack((cat_zsp, cat_zph))
cat_id = cat[idcol]

remove_patterns = [
	'PAUS_KSB_CloneZIDRandoms_',
	'_uniqK_Q-0.2',
	'_.pergaldz0.04',
	'K1000_SML_CloneZIDRandoms_',
	'_gama_McNaughtK_Q-0.7',
	'_pergaldz0.03',
	'_zph', '_', '.fits','5e6'
]

ran_phsp_list = []
ran_names = []
print 'readying..'
for rfi, rf in enumerate(ranlist):
	rf_ = rf
	for patt in remove_patterns:
		rf_ = rf_.replace(patt, '')
	rf_ = rf_.lower()
	if rf_ != 'unwindowed':
		rf_ += 'windowed randoms'
	else:
		rf_ += ' randoms'
	ran_names.append(rf_)

	ran = fopen(rf)
	if len(ran) == 0:
		print '%s is empty?'%rf
		continue
	if len(ran) > 1e7:
		ran = down(ran, 1e7/len(ran))

	ran = ran[np.argsort(ran['%s_cloneID'%magcol])]
	ran_id = ran['%s_cloneID'%magcol]
	ran_ids_bc = np.bincount(ran_id)
	ids = range(len(ran_ids_bc))
	ids = np.array(ids)[ran_ids_bc!=0]
	ran_ids_bc = ran_ids_bc[ran_ids_bc!=0]
	cat_zsp_i = cat_zsp[np.isin(cat_id, ids)]
	ran_ids_bc = ran_ids_bc[np.isin(ids, cat_id)]
	ran_zsp = np.repeat(cat_zsp_i, ran_ids_bc)
	ran_zph = ran['%s_cloneZ'%magcol][np.isin(ran_id, cat_id)]

	ran_phsp_list.append(np.column_stack((ran_zsp, ran_zph)))

print 'plotting..'
#c = ChainConsumer().add_chain(cat_phsp, parameters=['zspec','zphot'], name=name)
#for i, ran_phsp in enumerate(ran_phsp_list):
#	c.add_chain(ran_phsp, name=ranlist[i].replace('_',''))
#fig = c.plotter.plot()

#cm = plt.cm.nipy_spectral(np.linspace(0,1,len(ranlist)))
f, ax = plt.subplots(figsize=(7,6))
cat_h, bin_x, bin_y = np.histogram2d(cat_phsp[:,0], cat_phsp[:,1], bins=30, range=[[ll, nn(cat_zsp).max()],[ll, np.percentile(nn(cat_zph), 99.)]])
cat_h = 1.*cat_h/cat_h.sum()
cat_h = cat_h.T

fracs = (np.arange(13)/12.)[1:-1][::2]
def get_levels(Z, fracs):
	fZ = Z.flatten()
	fZ.sort()
	fZ = fZ[::-1]
	cumsum = np.cumsum(fZ)
	# interpolate cumulative probability back to grid-value at confidence fraction
	interp = scint.PchipInterpolator(cumsum, fZ)
	levels = [float(interp(frac)) for frac in fracs]
	levels.sort()
	return levels

c = ax.contour(midpoints(bin_x), midpoints(bin_y), cat_h,
	colors='k', label='galaxies', levels=get_levels(cat_h, fracs),
	zorder=10, alpha=1, lw=0.8)
h = [c.collections[0]]
l = ['galaxies']
for i in range(len(ranlist)):
	ran_phsp = ran_phsp_list[i]
	ran_h, bin_x, bin_y = np.histogram2d(ran_phsp[:,0], ran_phsp[:,1], bins=(bin_x,bin_y))
	ran_h = 1.*ran_h/ran_h.sum()
	ran_h = ran_h.T
	c = ax.contour(midpoints(bin_x), midpoints(bin_y), ran_h,
		colors='C%s'%(i), label=ranlist[i].replace('_',''),
		levels=get_levels(ran_h, fracs), alpha=0.8, linestyles='--', lw=1.2)
	h.append(c.collections[0])
	l.append(ran_names[i])
plt.xlabel(r'$z_{\rm{spec.}}$')
plt.ylabel(r'$z_{\rm{phot./rand.}}$')
zmax = dict(PAUS=1.2, GAMA=0.55)[name]
plt.xlim(-0.03, zmax)
plt.ylim(-0.03, zmax)
plt.plot([0,nn(cat_zsp).max()], [0,nn(cat_zsp).max()], 'r-', lw=0.6)
ax.legend(h, l, loc='upper left', frameon=0, fontsize=14)

plt.show()
plt.tight_layout()
#plt.subplots_adjust(top=0.8)
plt.savefig(outfile+'.png', bbox_inches='tight')
plt.savefig(outfile+'.pdf', bbox_inches='tight')



