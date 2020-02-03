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

f, ax1 = plt.subplots()
ax2 = ax1.twiny()
plt.sca(ax1)
diag = h5py.File('diagnostics_clonerandoms_%s.h5'%runid, 'r')
Om = diag['Om'].value
V_mid = diag['V_mid'][:]

windows = h5py.File('windows_%s.h5'%runid, 'r')['windows'][:]
chi = (3.*V_mid/Om)**(1./3.)
for i in np.random.choice(range(len(windows)), size=nwin):
#i = 0
#while i < nwin:
	#j = np.random.choice(range(len(windows)))
	#if not any(windows[j][:3] != 0):
	#	continue
	#plt.plot(get_z(chi), windows[i], alpha=0.7)
	#i += 1
	plt.plot(V_mid, windows[i], alpha=0.7)
#plt.yscale('log')
new_tick_locations = np.arange(1,30,step=4)*1e6
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
tick_function = lambda x: get_z((3.*x / Om)**(1./3.))
ax2.set_xticklabels(['%.3f'%i for i in tick_function(new_tick_locations)])
ax2.set_xlabel("$z$")
plt.xlabel('z')
plt.xlabel('V')
plt.show()
plt.tight_layout()
diag.close()

