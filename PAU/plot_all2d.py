# coding: utf-8
import matplotlib
matplotlib.use('Agg')
from functions import *
#run_startup()
dat = read_files('.dat', dire=sys.argv[1], asc=1)
dkeys = np.sort(dat.keys())
for k in dkeys:
	try:
		dict3d = popen(join(sys.argv[1], k.replace('.dat','.p')))
	except:
		print k, 'failed'
		continue
	r = dict3d['r'][:,0]
	r1 = np.append(-r[::-1], r)
	noise = dat[k]['noise']
	Pi = dict3d['Pi']
	midPi = midpoints(Pi)
	pairsdict = {}
	if 'wgp' in k:
		pairsdict['DS'] = dict3d['DS_3d']
		pairsdict['RS'] = dict3d['RS_3d']
	elif 'wgg' in k:
		pairsdict['DD'] = dict3d['DD3d']
		pairsdict['RD'] = dict3d['RD3d']
		pairsdict['DR'] = dict3d['DR3d']
		pairsdict['RR'] = dict3d['RR3d']

	normalise = 1
	if normalise:
		if 'wgp' in k:
			pairsdict['DS'] = dict3d['DS_3d']
			pairsdict['RS'] = dict3d['RS_3d'] * dict3d['DSntot'] / dict3d['RSntot']
		elif 'wgg' in k:
			pairsdict['DD'] = dict3d['DD3d']
			pairsdict['RD'] = dict3d['RD3d'] * dict3d['DDntot'] / dict3d['RDntot']
			pairsdict['DR'] = dict3d['DR3d'] * dict3d['DDntot'] / dict3d['DRntot']
			pairsdict['RR'] = dict3d['RR3d'] * dict3d['DDntot'] / dict3d['RRntot']

	construct = 1
	w = dict3d['w3d']
	if construct and 'wgp' not in k:
		w1 = pairsdict['DD'] - pairsdict['DR'] - pairsdict['RD'] + pairsdict['RR']
		w1 /= pairsdict['RR']
		w = [w, w1]
	else:
		w = [w]

	if 'wgp' in k:
		f, ax = plt.subplots(2, 2, figsize=(12,7))
	elif 'wgg' in k:
		f, ax = plt.subplots(3, 2, figsize=(12,10))
	ax = ax.flatten()
	#ax1 = ax[0].twinx()

	plt.sca(ax[0])
	for iw in w:
		wp = np.sum(iw * (Pi[1]-Pi[0]), axis=0)
		if 'wgp' in k:
			fac = r**0.8
			plt.ylabel(r'$r^{0.8}w_{\rm{g+}}(r_{p})$')
		else:
			fac = r
			plt.ylabel(r'$rw_{\rm{gg}}(r_{p})$')
		plt.semilogx(r, wp*fac)
	plt.axhline(0, c='k', ls=':')
	plt.xlabel('$r_{p}$')
	plt.title(k)
	#plt.sca(ax1)
	#ax1.tick_params(right=0)
	#plt.plot(r, noise, 'r:')

	w = w[0]
	plt.sca(ax[1])
	plt.yscale('symlog',linthreshy=10)
	plt.xscale('symlog',linthreshx=10)
	w1 = np.append(w[:,::-1], w, axis=1)
	#w1 = np.log10(abs(w1)) * np.sign(w1)
	#c = plt.imshow(w1, cmap='coolwarm', norm=MidpointNormalize(midpoint=0),
				#aspect='auto', origin='lower', extent=[-np.log10(r.max()),np.log10(r.max()),Pi.min(),Pi.max()])
	c = plt.pcolor(r1, Pi, np.log10(w1), cmap='PuRd')#, norm=MidpointNormalize(midpoint=0))
	cbar = plt.colorbar(c)
	cbar.ax.set_ylabel(r'log$(w)$')
	#cbar.ax.set_ylabel(r'$w \, / \, \sigma$')
	plt.ylabel(r'$\Pi$')

	mean_subtract = 0
	for i, pk in enumerate(np.sort(pairsdict.keys())):
		plt.sca(ax[2+i])
		plt.yscale('symlog',linthreshy=10)
		plt.xscale('symlog',linthreshx=10)
		pairs = pairsdict[pk]
		pairs = np.log10(pairs)
		pairs1 = pairs.copy()
		if mean_subtract:
			pairs -= pairs.mean(axis=0)
			cmap = 'coolwarm'
			cblab = 'log(%s) - mean'%pk
		else:
			cmap = 'viridis'
			cblab = 'log(%s)'%pk
		pairs = np.append(pairs[:,::-1], pairs, axis=1)
		#Pi1 = np.append(Pi[::-1], Pi)
		c = plt.pcolor(r1, Pi, pairs, cmap='coolwarm', norm=MidpointNormalize(midpoint=0))
					#aspect='auto', origin='lower', extent=[-np.log10(r.max()),np.log10(r.max()),Pi.min(),Pi.max()])
		#c = plt.imshow(pairs, cmap=cmap, aspect='auto', origin='lower', norm=MidpointNormalize(midpoint=0),
		#				extent=[-np.log10(r.max()),np.log10(r.max()),Pi.min(),Pi.max()])
						#vmin=-0.3, vmax=0.3)
		cbar = plt.colorbar(c)
		cbar.ax.set_ylabel(cblab)

		int_over_rp = 1
		if int_over_rp:
			plt.axhline(0, c='k', ls=':')
			sum_pairs = np.sum(pairs1, axis=1)
			axy = plt.gca().twiny()
			axy.plot(sum_pairs / sum_pairs.max(), midPi)
			axy.tick_params(labeltop=0)

	plt.tight_layout()
	plt.savefig(join(sys.argv[1], k.replace('.dat','.png')), bbox_inches='tight')
    
plt.show()





