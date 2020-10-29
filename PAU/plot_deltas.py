# coding: utf-8
from functions import *
from clone_randoms import *
if len(sys.argv) >= 2:
	deltas = read_files(sys.argv[1])
else:
	deltas = read_files('Delta', idnot='Unw')

plt.figure(figsize=(12,5))
cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(deltas.keys())))
for i, k in enumerate(np.sort(deltas.keys())):
    #if i < 1: continue
    d, D = deltas[k].T
    z = get_z(d)
    c, ls = cmap[i], '-'
    if i == len(deltas.keys())-1:
        c, ls = 'r', '-.'
    plt.plot(z, D, c=c, ls=ls)
#plt.axhline(1., c='k')
#plt.yscale('log')
plt.xlabel('$z$')
plt.ylabel(r'$\Delta(z)$')
plt.tight_layout()
plt.show()


