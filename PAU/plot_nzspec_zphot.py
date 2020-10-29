# coding: utf-8
import matplotlib
matplotlib.use('Agg')
from functions import *
if len(sys.argv) < 6:
	print '\nGive\n\
1: catalog path \n\
2: photo-z column name \n\
3: spec-z column name \n\
4: randoms path \n\
5: ID column name \n\
6: magnitude column name \n'
	sys.exit()

cat = fopen(sys.argv[1])
zphcol = sys.argv[2]
zspcol = sys.argv[3]
ran = fopen(sys.argv[4])
idcol = sys.argv[5]
magcol = sys.argv[6]

if 'PAU' in sys.argv[1]:
	bins = np.linspace(0.1, 1.2, 7)
	name = 'PAUS' + sys.argv[4].replace(sys.argv[1].replace('.fits',''), '').replace('CloneZIDRandoms_','').replace('.fits','')
	zmax = 1.4
else:
	bins = np.linspace(0, 0.5, 6)
	name = 'GAMA' + sys.argv[4].replace(sys.argv[1].replace('.fits',''), '').replace('CloneZIDRandoms_','').replace('.fits','')
	zmax = 0.6

f, ax = plt.subplots(len(bins)-1, sharex=True, figsize=(12,9.5))
plt.xlim(None, zmax)
#f, ax = plt.subplots(figsize=(12,9))
#ax = [ax]*(len(bins)-1)

for i in range(len(bins)-1):
	ax[i].axvline(bins[i], c='k', lw=0.5)
	ax[i].axvline(bins[i+1], c='k', lw=0.5)

if len(bins) < 11:
	cm = ['C%s'%i for i in range(10)]
else:
	cm = plt.cm.viridis(np.linspace(0.1,0.9,len(bins)-1))
for i in range(len(bins)-1):
	plt.sca(ax[i])
	cut = (cat[zphcol] > bins[i]) & (cat[zphcol] <= bins[i+1]) & (cat[zspcol] >= 0)
	subcat = cat[cut]
	subran = down(ran[np.isin(ran[magcol+'_cloneID'], subcat[idcol])], 0.1)
	h1, b1 = np.histogram(nn(subcat[zspcol]), bins='auto', density=1)
	h2, b2 = np.histogram(nn(subran[magcol+'_cloneZ']), bins='auto', density=1)
	l1, = plt.plot(midpoints(b1), h1, lw=0.8, c=cm[i])
	l2, = plt.plot(midpoints(b2), h2, c=l1.get_c(), ls='--')

labels = [r'$z_{\rm{spec.}}$', r'$z_{\rm{rand.}}$']
handles = [new_handle(ls='-', c='C0', lw=1.2),
		   new_handle(ls='--', c='C0', lw=1.4)]
ax[0].legend(handles, labels, loc='upper center', fontsize=20, frameon=0)
#ax[0].set_title(sys.argv[4])
ax[-1].set_xlabel(r'$z_{\rm{spec./rand.}}$', fontsize=20)
f.text(0.02, 0.5, r'$n(z_{\rm{spec.}} | z_{\rm{phot.}})$',
	va='center', ha='center', rotation='vertical', fontsize=20)
plt.tight_layout(pad=3)
plt.subplots_adjust(hspace=0)
plt.show()

plt.savefig(name+'_nzszp_marginals.png', bbox_inches='tight')
plt.savefig(name+'_nzszp_marginals.pdf', bbox_inches='tight')



