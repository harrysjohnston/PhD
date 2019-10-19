# coding: utf-8
from functions import *
from glob import glob
if len(sys.argv) < 4:
	print '\ngive:\n\
1: catalog path\n\
2: redshift column\n\
3: detection band column\n'
	sys.exit()

catpath = sys.argv[1]
zcol = sys.argv[2]
magcol = sys.argv[3]

cat = fopen(catpath)
f, ax = plt.subplots()
plt.hist(cat[zcol], alpha=0.4, bins='auto', normed=1)
for rf in glob('%s_Clone*'%catpath.replace('.fits','')):
    if 'quad' in rf: continue
    ran = fopen(rf)
    for col in ['%s_cloneZ'%magcol]:
        try:
            h, b = np.histogram(ran[col], bins=100, normed=1, range=(0,1.2))
            plt.plot(midpoints(b), h, label=rf.replace('%s_CloneZIDRandoms_'%catpath.replace('.fits',''),''), lw=1.6)
        except:
            pass
h,l = ax.get_legend_handles_labels()
plt.figlegend(h, l, loc='upper center', ncol=2)
plt.xlim(-0.03, 1.2)
plt.xlabel('z')
plt.show()
plt.tight_layout()
plt.subplots_adjust(top=0.8)

