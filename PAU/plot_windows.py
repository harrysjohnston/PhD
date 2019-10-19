# coding: utf-8
from functions import *
import h5py
from clone_randoms import *
if len(sys.argv) < 2:
	print '\ngive\n\
1: run ID\n'
	sys.exit()

runid = sys.argv[1]

diag = h5py.File('diagnostics_clonerandoms_%s.h5'%runid, 'r')
windows = h5py.File('windows_%s.h5'%runid, 'r')['windows'][:]
for i in np.random.choice(range(len(windows)), size=30):
    plt.plot(get_z(diag['d_mid'][:]), windows[i], alpha=0.7)
plt.yscale('log')
plt.xlabel('z')
plt.show()
plt.tight_layout()

