# coding: utf-8
import matplotlib
matplotlib.use('Agg')
from functions import *

t = Table.read('PAUS_wcolumns.fits')
kcid = h5py.File('PAUS_uniq.kcorrs', 'r')['ID'][:]
c1 = (t['mag_i'] > 18) & (t['mag_i'] < 22.5)
c2 = (t['mag_y'] > 18) & (t['mag_y'] < 22.5)
c3 = c1 | c2
c4 = np.isin(t['ref_id'], np.unique(kcid))
c5 = t['star_flag'] != 1
cut = c3 & c4 & c5

# make LePhare-based red_sequence
u = t['lp_mu']
g = t['lp_mg']
r = t['lp_mr']
i = t['lp_mi']
m = (2.1-1.8)/(-25--17.2)
c = 2.1 - m * -25
y = m*i + c
t['red_sequence_LePhare'] = u - i > y
t['blue_cloud_LePhare'] = u - i < y

# make Cigale-based red_sequence
t['red_sequence_Cigale_3cluster_normi'] = t['class3_normi'] == 3
t['blue_cloud_Cigale_3cluster_normi'] = t['class3_normi'] == 1
t['red_sequence_Cigale_3cluster'] = t['class3'] == 3
t['blue_cloud_Cigale_3cluster'] = t['class3'] == 1
t['red_sequence_Cigale_2cluster'] = t['class2'] == 1
t['blue_cloud_Cigale_2cluster'] = t['class2'] == 2

# more columns
t['comoving_mpc_bcnz_zb'] = MICEcosmo.comoving_distance(t['bcnz_zb']) * MICEcosmo.h
t['numeric_id'] = t['ref_id']
t['qz'][np.isnan(t['qz'])] = np.inf # so that inverse weight -> 0

t[cut].write('PAUS_cut.fits', overwrite=1)

# LePhare plot
tr = t['red_sequence_LePhare']
tb = t['blue_cloud_LePhare']
plt.figure()
plt.plot(np.sort(i), y[np.argsort(i)], 'C1-', lw=0.7)
#plt.scatter(i, u-i, s=0.01, alpha=0.1)
plt.scatter(i[tr], (u-i)[tr], c='r', s=0.01, alpha=0.1)
plt.scatter(i[tb], (u-i)[tb], c='b', s=0.01, alpha=0.1)
plt.xlabel('$i$')
plt.ylabel('$u-i$')
plt.xlim(-27, -12)
plt.ylim(-0.1, 4.4)
h = [
	new_handle(marker='.',c='r',ls='',label='PAUS W3 red'),
	new_handle(marker='.',c='b',ls='',label='PAUS W3 blue'),
	new_handle(c='C1',ls='-',label='$u-i=1.138-0.038i$')
]
l = ['PAUS W3 red', 'PAUS W3 blue', '$u-i=1.138-0.038i$']
plt.legend(h,l,fontsize=14,title='LePhare')
#plt.show()
plt.tight_layout()
plt.savefig('PAU_CMD_LePhare.png', bbox_inches='tight')

# Cigale plots
f, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(14, 4.8))
u = t['u_abs']
i = t['i_abs']

tr = t['red_sequence_Cigale_2cluster']
tb = t['blue_cloud_Cigale_2cluster']
plt.sca(ax[0])
plt.scatter(i[tr], (u-i)[tr], c='r', s=0.01, alpha=0.1)
plt.scatter(i[tb], (u-i)[tb], c='b', s=0.01, alpha=0.1)
plt.xlabel('$i$')
plt.ylabel('$u-i$')
plt.xlim(-27, -12)
plt.ylim(-0.1, 4.4)
h = [
	new_handle(marker='.',c='r',ls='',label='PAUS W3 red'),
	new_handle(marker='.',c='b',ls='',label='PAUS W3 blue')
]
l = ['PAUS W3 red', 'PAUS W3 blue']
plt.legend(h,l,fontsize=14,title='Cigale 2-cluster')
#plt.show()
plt.tight_layout()
#plt.savefig('PAU_CMD_Cigale_2cluster.png', bbox_inches='tight')

tr = t['red_sequence_Cigale_3cluster_normi']
tb = t['blue_cloud_Cigale_3cluster_normi']
tg = ~tr & ~tb
plt.sca(ax[1])
plt.scatter(i[tr], (u-i)[tr], c='r', s=0.01, alpha=0.1)
plt.scatter(i[tb], (u-i)[tb], c='b', s=0.01, alpha=0.1)
plt.scatter(i[tg], (u-i)[tg], c='g', s=0.01, alpha=0.1)
plt.xlabel('$i$')
#plt.ylabel('$u-i$')
plt.xlim(-27, -12)
plt.ylim(-0.1, 4.4)
h = [
	new_handle(marker='.',c='r',ls='',label='PAUS W3 red'),
	new_handle(marker='.',c='b',ls='',label='PAUS W3 blue'),
	new_handle(marker='.',c='g',ls='',label='PAUS W3 green')
]
l = ['PAUS W3 red', 'PAUS W3 blue', 'PAUS W3 green']
plt.legend(h,l,fontsize=14,title='Cigale 3-cluster')
#plt.show()
plt.tight_layout()
plt.savefig('PAU_CMD_Cigale.png', bbox_inches='tight')

