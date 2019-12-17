# coding: utf-8
from functions import *
import h5py
from clone_randoms import *
if len(sys.argv) < 2:
	print '\ngive\n\
1: run ID\n\
2: number of windows to plot (opt, default=50)\n'
	sys.exit()
if len(sys.argv) > 2:
	nwin = int(sys.argv[2])
else:
	nwin = 50

runid = sys.argv[1]

plt.figure()
diag = h5py.File('diagnostics_clonerandoms_%s.h5'%runid, 'r')
windows = h5py.File('windows_%s.h5'%runid, 'r')['windows'][:]
chi = (3.*diag['V_mid'][:]/diag['Om'].value)**(1./3.)
for i in np.random.choice(range(len(windows)), size=nwin):
    plt.plot(get_z(chi), windows[i], alpha=0.7)
    #plt.plot(diag['V_mid'][:], windows[i], alpha=0.7)
#plt.yscale('log')
plt.xlabel('z')
#plt.xlabel('V')
plt.show()
plt.tight_layout()
diag.close()

