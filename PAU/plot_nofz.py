# coding: utf-8
from functions import *
from glob import glob
if len(sys.argv) < 5:
	print '\ngive:\n\
1: catalog path\n\
2: redshift column\n\
3: detection band column\n\
4: ID column\n'
	sys.exit()

catpath = sys.argv[1]
zcol = sys.argv[2]
magcol = sys.argv[3]
idcol = sys.argv[4]

cat = fopen(catpath)
subcat = cat[cat[zcol] > 0.1]
f, ax = plt.subplots()
plt.hist(cat[zcol], alpha=0.4, bins=60, normed=1)
ranlist = np.sort(glob('%s_Clone*'%catpath.replace('.fits','')))

for rf in ranlist:
	if ('quad' in rf
		or '6e6' in rf
		or '1e6' in rf
		#or 'Unw' not in rf
		#or '8bin' not in rf
		#or '4e6' not in rf
		) :
		ranlist = np.delete(ranlist, list(ranlist).index(rf))
		continue
	if 'Q' in rf:
		Q = float(rf.split('Q')[1].replace('.fits', '').replace('_cz',''))
		if Q != -0.2:
			#pass
			ranlist = np.delete(ranlist, list(ranlist).index(rf))
cm = plt.cm.coolwarm(np.linspace(0,1,len(ranlist)))

for rfi, rf in enumerate(ranlist):
	ls = '-'
	if 'Q' in rf:
		Q = float(rf.split('Q')[1].replace('.fits', '').replace('_cz',''))
		if Q != 0.0:
			pass
			#ranlist = np.delete(ranlist, list(ranlist).index(rf))
		if Q < 0:
			ls = '--'

	ran = fopen(rf)
	if len(ran) == 0:
		print '%s is empty?'%rf
		continue
	c = cm[rfi]
	c = None
	#highz_clones = np.isin(ran['%s_cloneID'%magcol], subcat[idcol])
	for col in ['%s_cloneZ'%magcol]:
		try:
			h, b = np.histogram(ran[col], bins=100, normed=1, range=(0,ran[col].max()))
			l, = plt.plot(midpoints(b), h, label=rf.replace('%s_CloneZIDRandoms_'%catpath.replace('.fits',''),'').replace('.fits',''), lw=1.6, c=c, ls=ls)
			print rf
		except:
			pass
h,l = ax.get_legend_handles_labels()
if len(ranlist) > 6:
	ncol = len(ranlist)//6
else:
	ncol = 1
plt.figlegend(h, l, loc='upper center', ncol=ncol)
#plt.xlim(-0.03, ran[col].max())
plt.xlabel('z')
plt.show()
plt.tight_layout()
plt.subplots_adjust(top=0.8)

