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

f, ax = plt.subplots(2, 2, sharex=True, figsize=(10, 9))
ax1 = ax[0]
ax2 = ax[1]
#axins = inset_axes(ax1, loc='lower left', width='75%', height='30%', borderpad=3)
#axins.set_xscale('log')
ax1[0].set_xscale('log')
ax2[0].set_yscale('log')
ax2[1].set_yscale('log')
for a in ax1:
	a.axhline(0, c='k', lw=0.6)
#axins.axhline(0, c='k', lw=0.6)
kw = {'ha':'center','va':'center','fontsize':14}
if 'qz' in args.paus:
	if args.paper:
		f.text(0.5, 0.98, r'PAUS W3 -- best 50% Qz in $0.1<z_{\rm{phot.}}<0.8$', **kw)
	else:
		f.text(0.5, 0.98, 'PAUS W3 -- best 50%% Qz (%s)'%args.paus.replace('OUTPUTS_PAUS_','').replace('_qz50','').replace('_', '--'), **kw)
else:
	if args.paper:
		f.text(0.5, 0.98, r'PAUS W3 -- all galaxies in $0.1<z_{\rm{phot.}}<0.8$', **kw)
	else:
		f.text(0.5, 0.98, 'PAUS W3 -- all galaxies (%s)'%args.paus.replace('OUTPUTS_PAUS_','').replace('_', '--'), **kw)

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
		mark = 'p'
		lst = '-'
		lw = 1.3
		elw = 1.5
		cs = 4
		mfc = None
		ms = 6
	elif cat == 'GAMA':	
		mark = 'o'
		lst = ':'
		lw = 0.6
		elw = None
		mfc = 'none'
	for dire in ddict[cat].keys():
		if cat == 'GAMA':
			if 'ztonry' not in dire.lower().replace('_','') or 'zphot' in dire:
				mark = 'v'
				mfc = None
				ms = 7
				cs = 2
			else:
				ms = 7
				cs = 1.3
		for col in ['red', 'blue']:
			lab = '%s %s'%(cat, col)
			if col == 'red':
				s = 0.98
				acol = 0
			elif col == 'blue':
				s = 1.02
				acol = 1
			if cat == 'GAMA':
				if 'ztonry' not in dire.lower().replace('_','') or 'zphot' in dire:
					s *= 1.01
					lab += r' $z_{\rm{phot.}}$'
					if col == 'red': coli = 'salmon'
					if col == 'blue': coli = 'deepskyblue'
				else:
					s *= 0.99
					if col == 'red': coli = 'firebrick'
					if col == 'blue': coli = 'darkblue'
			else:
				coli = col
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
					ax1[acol].errorbar(r*s, r**rpp*w, r**rpp*e, fmt=mark, mfc=mfc, ms=ms, ls=lst, lw=lw, elinewidth=elw, c=coli, capsize=cs, label=lab)
					#axins.errorbar(r*s, r**rpp*w, r**rpp*e, fmt=mark, mfc=mfc, ms=ms, lw=lw, elinewidth=elw, c=col, capsize=cs, label=lab)
				if sig == 'wgg':
					try:
						r, w, e = dat['rnom'], dat['wgg'], dat['wgg_jackknife_err']
					except KeyError:
						r, w = dat['rnom'], dat['wgg']
						e = np.zeros_like(r)
					ax2[acol].errorbar(r*s, w, e, fmt=mark, mfc=mfc, ms=ms, ls=lst, lw=lw, elinewidth=elw, c=coli, capsize=cs, label=lab)

x1, x2, y1, y2 = 0.12, 5.3, -0.16, 0.15
#axins.set_xlim(x1, x2)
#axins.set_ylim(y1, y2)
#axins.set_xticks([0.2, 1, 3])
#axins.set_xticklabels(['0.2','1.0','3.0'], fontdict={'fontsize':11})
##axins.set_yticklabels([])
#mark_inset(ax1, #axins, loc1=1, loc2=3, fc='none', ec='k', alpha=0.4, lw=0.5)
h1, l1 = ax1[0].get_legend_handles_labels()
h2, l2 = ax1[1].get_legend_handles_labels()
h = h1 + h2
l = l1 + l2
ax2[0].legend(h, l, ncol=1, loc='best', fontsize=13, frameon=0)
ax1[0].set_ylim(-0.9, 1.4)
ax1[1].set_ylim(-0.6, 0.5)
for a in ax2:
	a.set_ylim(0.9, 1.2e3)

maxrp = 18.
xlim = ax1[0].get_xlim()
for a in ax.flatten():
	a.axvspan(maxrp, xlim[1], facecolor='none', lw=0, edgecolor='grey', hatch='////', alpha=0.5)
ax1[0].set_xlim(*xlim)

if not args.wgx:
	ax1[0].set_ylabel(r'$r_{p}^{0.8}w_{\rm{g+}}\,[h^{-1}\rm{Mpc}]^{1.8}$')
else:
	ax1[0].set_ylabel(r'$r_{p}^{0.8}w_{\rm{g\times}}\,[h^{-1}\rm{Mpc}]^{1.8}$')
ax2[0].set_ylabel(r'$w_{\rm{gg}}\,[h^{-1}\rm{Mpc}]$')
ax2[0].set_xlabel(r'$r_{p}\,\,[h^{-1}\rm{Mpc}]$')
ax2[1].set_xlabel(r'$r_{p}\,\,[h^{-1}\rm{Mpc}]$')
plt.tight_layout()
plt.subplots_adjust(hspace=0, top=0.96)
plt.show()
outname = 'PAUvsGAMA_figure_%s_%s'%(basename(args.paus).replace('OUTPUTS_','').strip('/'), args.gama[1].replace('OUTPUTS_','').strip('/'))

if not args.wgx:
	plt.savefig('master_figures/'+outname+'.png', bbox_inches='tight')
	plt.savefig('master_figures/'+outname+'.pdf', bbox_inches='tight')




