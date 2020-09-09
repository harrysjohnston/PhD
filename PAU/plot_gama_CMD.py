# coding: utf-8
from functions import *
cat = fopen('../SMLambdarApMatchedPhotom.fits')
cat = cat[cat['PETROMAG_R']<=19.8]
g = cat['absmag_g']
r = cat['absmag_r']
u = cat['absmag_u']
i = cat['absmag_i']
t = Table(cat)
y = np.ones_like(g)*0.66
t['red_sequence'] = g-r > y
tr = t['red_sequence']
plt.figure()
#plt.plot(np.sort(r), y[np.argsort(r)], 'C1-', lw=0.7)
#plt.scatter(i, u-i, s=0.01, alpha=0.1)
plt.scatter(i[tr], (u-i)[tr], c='r', s=0.01, alpha=0.1)
plt.scatter(i[~tr], (u-i)[~tr], c='b', s=0.01, alpha=0.1)
plt.xlabel('$i$')
plt.ylabel('$u-i$')
plt.xlim(-27, -12)
plt.ylim(-0.1, 4.1)
h = [
    new_handle(marker='.',c='r',ls='',label='GAMA red'),
    new_handle(marker='.',c='b',ls='',label='GAMA blue')]
    #new_handle(c='C1',ls='-',label='$g-r=0.66$')]
l = ['GAMA red', 'GAMA blue']#, '$g-r=0.66$']
plt.legend(h,l,fontsize=14)
plt.tight_layout()
plt.show()
plt.savefig('GAMA_CMD.png', bbox_inches='tight')

