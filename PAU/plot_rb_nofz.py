# coding: utf-8
import matplotlib
matplotlib.use('Agg')
from functions import *
from glob import glob
if len(sys.argv) < 6:
	print '\ngive:\n\
1: catalog path\n\
2: redshift column\n\
3: detection band column\n\
4: ID column\n\
5: red sequence (bool) column\n\
6: randoms paths\n'
	sys.exit()

catpath = sys.argv[1]
zcol = sys.argv[2]
magcol = sys.argv[3]
idcol = sys.argv[4]

cat = fopen(catpath)
is_red = cat[sys.argv[5]]

ranlist = sys.argv[6:]

#ranlist = glob('%s_Clone*'%catpath.replace('.fits',''))
#for rf in ranlist:
#	if ('quad' in rf
#		#or '6e6' in rf
#		#or '1e6' in rf
#		#or 'Unw' not in rf
#		#or '8bin' not in rf
#		#or '4e6' not in rf
#		) :
#		ranlist = np.delete(ranlist, list(ranlist).index(rf))
#		continue
#	if 'Q' in rf:
#		Q = float(rf.split('Q')[1].replace('.fits', '').replace('_cz',''))
#		if 'Lambdar' in catpath and Q != -0.7:
#			#pass
#			ranlist = np.delete(ranlist, list(ranlist).index(rf))
#		elif 'PAUS' in catpath and Q != -0.2:
#			#pass
#			ranlist = np.delete(ranlist, list(ranlist).index(rf))

if 'PAU' in catpath:
	catlab = 'PAUS W3'
	zlab = r'$z_{\rm{phot.}}$'
elif 'SML' in catpath or 'GAMA' in catpath:
	catlab = 'GAMA'
	zlab = r'$z_{\rm{spec.}}$'

f, ax = plt.subplots()
plt.hist(cat[zcol][is_red], color='r', alpha=0.3, bins=50, density=1, label=catlab+' red')
plt.hist(cat[zcol][~is_red], color='b', alpha=0.3, bins=50, density=1, label=catlab+' blue')

for rf in ranlist:
	ran = fopen(rf)
	r_clones = np.isin(ran['%s_cloneID'%magcol], cat[idcol][is_red])
	b_clones = np.isin(ran['%s_cloneID'%magcol], cat[idcol][~is_red])

	for col in ['%s_cloneZ'%magcol]:

		if 'Unwindowed' in rf:
			lst, lw = '--', 1.6
			lab = 'unwindowed'
		else:
			lst, lw = '-', 1.4
			lab = 'windowed'
		#lst, lw = '-', 1.4
		#lab = rf.replace(catpath.replace('.fits',''),'').replace('_CloneZIDRandoms_','').replace('.fits','')

		hr, br = np.histogram(ran[col][r_clones], bins=100, density=1, range=(0.02,ran[col][r_clones].max()))
		lr, = plt.plot(midpoints(br), hr, c='r', lw=lw, ls=lst, label=lab+' red',
						#label='%s--red'%rf.replace('%s_CloneZIDRandoms_'%catpath.replace('.fits',''),'').replace('.fits','')
					)

		hb, bb = np.histogram(ran[col][b_clones], bins=100, density=1, range=(0.02,ran[col][b_clones].max()))
		lb, = plt.plot(midpoints(bb), hb, c='b', lw=lw, ls=lst, label=lab+' blue',
						#label='%s--blue'%rf.replace('%s_CloneZIDRandoms_'%catpath.replace('.fits',''),'').replace('.fits','')
					)

h,l = ax.get_legend_handles_labels()
sort = list(np.argsort(l))
h = list(np.array(h)[sort])
l = list(np.array(l)[sort])
plt.legend(h, l, loc='best', ncol=1, fontsize=11, frameon=0)
plt.xlim(-0.03, None)
plt.xlabel(zlab)
plt.ylabel('PDF')
plt.show()
plt.tight_layout()
#plt.subplots_adjust(top=0.8)
catlab1 = catlab.split(' ')[0]
plt.savefig('%s_rb_nofz_figure.pdf'%catlab1, bbox_inches='tight')
plt.savefig('%s_rb_nofz_figure.png'%catlab1, bbox_inches='tight')



