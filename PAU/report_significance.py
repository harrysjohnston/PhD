from functions import *

Njk = int(sys.argv[1])
off_diag = 1
plot = 1
drop_largerp = 1
print '\n'
print 'Njk =', Njk
print 'off-diagonals =', bool(off_diag), 'hard-coded!'
print 'plot =', bool(plot), 'hard-coded!'
print 'exclude largest-rp bin =', bool(drop_largerp), 'hard-coded!'

rep_dict = {}
for dire in sys.argv[2:]:
	rep_dict[dire] = {}
	try:
		wgp_r = ascii.read(join(dire, 'wgp_red.dat'))
		wgp_b = ascii.read(join(dire, 'wgp_blue.dat'))
		wgp_r_cov = np.loadtxt(join(dire, 'wgp_red_wgplus.cov'))
		wgx_r_cov = np.loadtxt(join(dire, 'wgp_red_wgcross.cov'))
		wgp_b_cov = np.loadtxt(join(dire, 'wgp_blue_wgplus.cov'))
		wgx_b_cov = np.loadtxt(join(dire, 'wgp_blue_wgcross.cov'))
	except IOError:
		print dire, 'failed'
		continue

	if drop_largerp:
		cut = wgp_r['rnom'] < 18.
		wgp_r = wgp_r[cut]
		wgp_b = wgp_b[cut]
		wgp_r_cov = wgp_r_cov[cut][:, cut]
		wgp_b_cov = wgp_b_cov[cut][:, cut]
		wgx_r_cov = wgx_r_cov[cut][:, cut]
		wgx_b_cov = wgx_b_cov[cut][:, cut]

	if not off_diag:
		wgp_r_cov = np.diag(wgp_r_cov) * np.identity(len(wgp_r_cov))
		wgp_b_cov = np.diag(wgp_b_cov) * np.identity(len(wgp_b_cov))
		wgx_r_cov = np.diag(wgx_r_cov) * np.identity(len(wgx_r_cov))
		wgx_b_cov = np.diag(wgx_b_cov) * np.identity(len(wgx_b_cov))

	print '\n', dire, ':'
	for signal in ['wgplus', 'wgcross']:
		rep_dict[dire][signal] = {}
		w_r = {'wgplus':wgp_r, 'wgcross':wgp_r}[signal]
		w_r_cov = {'wgplus':wgp_r_cov, 'wgcross':wgx_r_cov}[signal]
		w_b = {'wgplus':wgp_b, 'wgcross':wgp_b}[signal]
		w_b_cov = {'wgplus':wgp_b_cov, 'wgcross':wgx_b_cov}[signal]

		hf_r = hartlap(Njk, len(w_r[signal]))
		chi2_r = dotdot(w_r[signal], hf_r * np.linalg.inv(w_r_cov), w_r[signal])
		pval_r = chi2pval(chi2_r, len(w_r[signal]))
		sigma_r = p2sigma(pval_r)

		hf_b = hartlap(Njk, len(w_b[signal]))
		chi2_b = dotdot(w_b[signal], hf_b * np.linalg.inv(w_b_cov), w_b[signal])
		pval_b = chi2pval(chi2_b, len(w_b[signal]))
		sigma_b = p2sigma(pval_b)

		print '%s,red: p = %.2f | %.2f sigma (HF = %.2f)'%(signal,pval_r,sigma_r,hf_r)
		print '%s,blue: p = %.2f | %.2f sigma (HF = %.2f)'%(signal,pval_b, sigma_b,hf_b)

		rep_dict[dire][signal]['red'] = (chi2_r, pval_r, sigma_r)
		rep_dict[dire][signal]['blue'] = (chi2_b, pval_b, sigma_b)

	print '\n'

if plot:
	f, ax = plt.subplots(2, 3, sharey='row', sharex=True, figsize=(12, 8))
	axis_x = [iter(range(100)),
			  iter(range(100)),
			  iter(range(100))]
	labels = [[],[],[]]
	for dire in np.sort(sys.argv[2:]):
		if 'LePhare' in dire: acol = 0
		elif 'Cigale_3cluster' in dire: acol = 1
		elif 'Cigale_2cluster' in dire:	acol = 2
		x = next(axis_x[acol])
		dire_lab = basename(dire)
		for pattern in ['OUTPUTS_PAUS_','LePhare_','Cigale_2cluster_','Cigale_3cluster_','normi_']:
			dire_lab = dire_lab.replace(pattern, '')
		dire_lab = dire_lab.replace('fibonacci_', 'dyn-$\Pi-$')
		dire_lab = dire_lab.replace('uniform_', 'unif-$\Pi-$')
		dire_lab = dire_lab.replace('unwindowed', 'U')
		dire_lab = dire_lab.replace('windowed', 'W')
		if '_qz50' in dire_lab: dire_lab = 'Qz50 ' + dire_lab.replace('_qz50', '')
		labels[acol].append(dire_lab.replace('_',' '))
		for signal in ['wgplus', 'wgcross']:
			if signal == 'wgplus': arow = 0
			if signal == 'wgcross': arow = 1
			sigma_r = rep_dict[dire][signal]['red'][2]
			sigma_b = rep_dict[dire][signal]['blue'][2]
			ax[arow, acol].plot([x], [sigma_r], 'rs', mfc='none', mec='r')
			ax[arow, acol].plot([x], [sigma_b], 'bs', mfc='none', mec='b')
	ax[1, 0].set_xticks(range(len(labels[0])))
	ax[1, 1].set_xticks(range(len(labels[1])))
	ax[1, 2].set_xticks(range(len(labels[2])))
	ax[1, 0].set_xticklabels(labels[0], rotation=60, ha='right', rotation_mode='anchor', fontsize=9)
	ax[1, 1].set_xticklabels(labels[1], rotation=60, ha='right', rotation_mode='anchor', fontsize=9)
	ax[1, 2].set_xticklabels(labels[2], rotation=60, ha='right', rotation_mode='anchor', fontsize=9)
	ax[0, 0].set_ylabel(r'$w_{\rm{g+}}$ detection $[\sigma]$', fontsize=14)
	ax[1, 0].set_ylabel(r'$w_{\rm{g\times}}$ detection $[\sigma]$', fontsize=14)
	ax[0, 0].set_title('LePhare')
	ax[0, 1].set_title('Cigale 3-cluster')
	ax[0, 2].set_title('Cigale 2-cluster')
	for a in ax.flatten():
		a.tick_params(which='minor', bottom=0, top=0)
		a.set_yticks(range(10))
		a.autoscale()
		a.axhline(2, ls=':', c='k')
		a.axhline(3, ls='--', c='k')
	for i in range(len(labels)):
		for a in ax[:, i]:
			for j in range(len(labels[i])):
				if 'Qz50' in labels[i][j]:
					a.axvspan(j-0.5, j+0.5, color='grey', alpha=0.2, lw=0)
				if 'dyn' in labels[i][j]:
					a.axvspan(j-0.5, j+0.5, facecolor='none', edgecolor='cyan', hatch='\\\\\\\\', alpha=0.6, lw=0)
				if 'zph' in labels[i][j]:
					a.axvspan(j-0.5, j+0.5, facecolor='none', edgecolor='cyan', hatch='////', alpha=0.6, lw=0)
	plt.tight_layout()
	plt.subplots_adjust(hspace=0, wspace=0)
	plt.savefig('summary_plots/PAUS_IA_summary_Njk%s_fig.pdf'%Njk, bbox_inches='tight')
	plt.savefig('summary_plots/PAUS_IA_summary_Njk%s_fig.png'%Njk, bbox_inches='tight')
	plt.show()
	












