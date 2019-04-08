# coding: utf-8
from functions import *

if len(sys.argv[1:]) < 4:
	print "\n1: catalogue\
		   \n2: magnitude col\
		   \n3: magnitude lim\
		   \n4: redshift col"
	sys.exit()

cat = fopen(sys.argv[1])
mcol = sys.argv[2]
mlim = sys.argv[3]
zcol = cat[sys.argv[4]]

d = cosmo.comoving_distance(zcol)
dm = cat['%s_fl%s_dmax'%(mcol, mlim)]
zm = cat['%s_fl%s_zmax'%(mcol, mlim)]
mask = (dm > 0) & (zm > 0)
f, ax = plt.subplots()
axx = ax.twinx()
#axy = ax.twiny()
plt.sca(ax)
dmh = myhist(dm[mask], normed=0, color='C1', label='comoving distance')
#plt.sca(axy)
zmh = myhist(cosmo.comoving_distance(zm[mask]), normed=0, label='redshift')
ax.set_xlabel('$d_{max}$')
#axy.set_xlabel('$z_{max}$')
#h = [a.get_legend_handles_labels()[0][0] for a in [ax, axy]]
#plt.legend(handles=h, loc='best')
plt.tight_layout()
myhist(d, histtype='bar', alpha=0.5, normed=0)
plt.sca(axx)
dl = np.sort(dm[(dm < (d.max().value)*3.)])
plt.plot(dl, (dl/dl.max())**2., 'r:')
plt.xlim(None, d.max().value*3.)


