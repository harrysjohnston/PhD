import matplotlib
matplotlib.use('Agg')
from functions import *

cpath = sys.argv[1]
t = Table.read(cpath)

#w = t['pgm']**2. / ((t['pgm']*0.3)**2. + t['dem']**2.)
#w[t['pgm'] < 0.1] = 0.
inv_w = 0.25**2. + (t['dem']/t['pgm'])**2.
w = 1. / inv_w
w /= w.sum()
qz_w = 1. / t['qz']
qz_w /= qz_w.sum()

t['KSBweight'] = 1. / w
t['KSB_qz_weight'] = 1. / (w * qz_w)

zcut = (t['bcnz_zb'] > 0.1) & (t['bcnz_zb'] < 0.8)

for rc in ['_LePhare','_Cigale_2cluster','_Cigale_3cluster', '_Cigale_3cluster_normi']:
	redcut = t['red_sequence'+rc]
	bluecut = t['blue_cloud'+rc]
	sample_cuts = [zcut, (zcut&redcut), (zcut&bluecut)]
	sample_labels = ['total', 'red', 'blue']
	for lab, cut in zip(sample_labels, sample_cuts):
		lab += rc

		mcorr = np.nan * np.ones(len(t))
		g1 = np.nan * np.ones(len(t))
		g2 = np.nan * np.ones(len(t))

		mcorr = np.sum(t['mum'][cut]*w[cut]) / np.sum(w[cut])
		g1[cut] = (t['e1c'] / t['pgm'] / mcorr)[cut]
		g2[cut] = (t['e2c'] / t['pgm'] / mcorr)[cut]

		t['mcorr_'+lab] = mcorr
		t['g1_'+lab] = g1
		t['g2_'+lab] = g2
		t[lab+'_ok'] = np.isfinite(g1) & np.isfinite(g2) & (mcorr > 0.1) & (t['pgm'] > 0.1)

		plt.figure() ; plt.title(lab)
		plt.scatter(g1, t['e1c'], s=1, alpha=0.1)
		plt.scatter(g2, t['e2c'], s=1, alpha=0.1)
		myhist(nn(g1))
		myhist(nn(g2), ls='--')
		plt.tight_layout()
		plt.savefig('KSB_shapes_%s.png'%lab, bbox_inches='tight')
		plt.close()

		print lab, 'g1 frac not finite:', frac(~np.isfinite(g1))
		print lab, 'g2 frac not finite:', frac(~np.isfinite(g2))

t.write(sys.argv[1], overwrite=1)

