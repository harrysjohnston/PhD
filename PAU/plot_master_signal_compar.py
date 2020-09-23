# coding: utf-8
import matplotlib
matplotlib.use('Agg')
from functions import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'-gama',
		type=str,
		nargs='*',
		help='paths to GAMA correlations -- up to 2 paths')
	parser.add_argument(
		'-paus',
		type=str,
		help='path to PAUS correlations')
	parser.add_argument(
		'-wgx',
		type=int,
		default=0,
		help='1 = do plots for wgx (default=0)')
	parser.add_argument(
		'-paper',
		type=int,
		default=0,
		help='1 = do plots the paper (i.e. remove extra labels etc.)')
	args = parser.parse_args()
signal_id = ['wgplus','wgcross'][args.wgx]

f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(9,9))
axins = inset_axes(ax1, loc='lower left', width='75%', height='30%', borderpad=3)
axins.set_xscale('log')
ax1.set_xscale('log')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax1.axhline(0, c='k', lw=0.6)
axins.axhline(0, c='k', lw=0.6)
if 'qz' in args.paus:
	ax1.set_title('PAUS W3 -- best 50%% Qz (%s)'%args.paus.replace('OUTPUTS_PAUS_','').replace('_qz50','').replace('_', '--'))
	if args.paper:
		ax1.set_title(r'PAUS W3 -- best 50% Qz in $0.1<z_{\rm{phot.}}<0.8$')
else:
	ax1.set_title('PAUS W3 -- all galaxies (%s)'%args.paus.replace('OUTPUTS_PAUS_','').replace('_', '--'))
	if args.paper:
		ax1.set_title(r'PAUS W3 -- all galaxies in $0.1<z_{\rm{phot.}}<0.8$')

ddict = {}
ddict['GAMA'] = {}
for gama_dir in args.gama:
	ddict['GAMA'][gama_dir] = read_files('.dat', dire=gama_dir, asc=1)
ddict['PAUS'] = {}
ddict['PAUS'][args.paus] = read_files('.dat', dire=args.paus, asc=1)

rpp = 0.8
s3 = splitN(3, 0.03)

for cat in ['PAUS', 'GAMA']:
	if cat == 'PAUS':
		mark = '.'
		lst = '-'
		lw = 1.3
		elw = 1.5
		cs = 0
		mfc = None
		ms = 2
	elif cat == 'GAMA':	
		mark = 'o'
		lst = ':'
		lw = 0.6
		elw = None
		mfc = 'none'
	for dire in ddict[cat].keys():
		if cat == 'GAMA':
			if 'ztonry' not in dire.lower().replace('_','') or 'zphot' in dire:
				ms = 12
				cs = 4
			else:
				ms = 7
				cs = 1.3
		else:
			ms = None
			cs = 0
		for col in ['red', 'blue']:
			lab = '%s %s'%(cat, col)
			if col == 'red':
				s = 0.98
			elif col == 'blue':
				s = 1.02
			if cat == 'GAMA':
				if 'ztonry' not in dire.lower().replace('_','') or 'zphot' in dire:
					s *= 1.01
					lab += r' $z_{\rm{phot.}}$'
				else:
					s *= 0.99
			for sig in ['wgp', 'wgg']:
				try: dat = ddict[cat][dire][sig+'_'+col+'.dat']
				except:
					print '\n', cat, dire, sig+'_'+col+'.dat', 'failed!'
					continue
				if cat == 'GAMA': dat = dat[:-1]
				if sig == 'wgp':
					try:
						r, w, e = dat['rnom'], dat[signal_id], dat['%s_jackknife_err'%signal_id]
					except KeyError:
						r, w = dat['rnom'], dat[signal_id]
						e = np.zeros_like(r)
					ax1.errorbar(r*s, r**rpp*w, r**rpp*e, fmt=mark, mfc=mfc, ms=ms, ls=lst, lw=lw, elinewidth=elw, c=col, capsize=cs, label=lab)
					axins.errorbar(r*s, r**rpp*w, r**rpp*e, fmt=mark, mfc=mfc, ms=ms, lw=lw, elinewidth=elw, c=col, capsize=cs, label=lab)
				if sig == 'wgg':
					try:
						r, w, e = dat['rnom'], dat['wgg'], dat['wgg_jackknife_err']
					except KeyError:
						r, w = dat['rnom'], dat['wgg']
						e = np.zeros_like(r)
					ax2.errorbar(r*s, w, e, fmt=mark, mfc=mfc, ms=ms, ls=lst, lw=lw, elinewidth=elw, c=col, capsize=cs, label=lab)

x1, x2, y1, y2 = 0.12, 5.3, -0.16, 0.15
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticks([0.2, 1, 3])
axins.set_xticklabels(['0.2','1.0','3.0'], fontdict={'fontsize':11})
#axins.set_yticklabels([])
mark_inset(ax1, axins, loc1=1, loc2=3, fc='none', ec='k', alpha=0.4, lw=0.5)
h, l = ax1.get_legend_handles_labels()
ax2.legend(h, l, ncol=3, loc='best', fontsize=13, frameon=0)
ax1.set_ylim(-1.8, 1.4)
ax2.set_ylim(0.9, None)

maxrp = 18.
xlim = ax1.get_xlim()
ax1.axvspan(maxrp, xlim[1], facecolor='none', lw=0, edgecolor='grey', hatch='////', alpha=0.5)
ax2.axvspan(maxrp, xlim[1], facecolor='none', lw=0, edgecolor='grey', hatch='////', alpha=0.5)
ax1.set_xlim(*xlim)

if not args.wgx:
	ax1.set_ylabel(r'$r_{p}^{0.8}w_{\rm{g+}}\,[h^{-1}\rm{Mpc}]^{1.8}$')
else:
	ax1.set_ylabel(r'$r_{p}^{0.8}w_{\rm{g\times}}\,[h^{-1}\rm{Mpc}]^{1.8}$')
ax2.set_ylabel(r'$w_{\rm{gg}}\,[h^{-1}\rm{Mpc}]$')
ax2.set_xlabel(r'$r_{p}\,\,[h^{-1}\rm{Mpc}]$')
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.show()
outname = 'PAUvsGAMA_figure_%s_%s'%(basename(args.paus).replace('OUTPUTS_','').strip('/'), args.gama[1].replace('OUTPUTS_','').strip('/'))

if not args.wgx:
	plt.savefig('master_figures/'+outname+'.png', bbox_inches='tight')
	plt.savefig('master_figures/'+outname+'.pdf', bbox_inches='tight')




