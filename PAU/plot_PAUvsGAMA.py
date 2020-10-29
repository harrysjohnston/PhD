# coding: utf-8
from functions import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

if len(sys.argv)<3:
	print '\n\
Give:\n\
1: "qz" or "tot"\n\
2: "004" or "009"\n'
	sys.exit()

f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(9,9))
axins = inset_axes(ax1, loc='lower left', width='75%', height='30%', borderpad=3)
axins.set_xscale('log')
ax1.set_xscale('log')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax1.axhline(0, c='k', lw=0.6)
axins.axhline(0, c='k', lw=0.6)
pz = sys.argv[2]
gama_ia = read_files('vs', dire='../Wcorr_pz%sfullGAMA_rb_c0.66/to_plot/'%pz)
gama_ia[r'GR ~ $z_{\rm{phot.}}$'] = gama_ia['highZ_vs_highZ_Red']
gama_ia[r'GB ~ $z_{\rm{phot.}}$'] = gama_ia['highZ_vs_highZ_Blue']
del gama_ia['highZ_vs_highZ_Red'], gama_ia['highZ_vs_highZ_Blue']
gama_ia = merge_two_dicts(gama_ia, read_files('vs', dire='../Wcorr_fullGAMA_rb_c0.66/to_plot/'))
gama_ia['GAMA red'] = gama_ia['highZ_vs_highZ_Red']
gama_ia['GAMA blue'] = gama_ia['highZ_vs_highZ_Blue']
del gama_ia['highZ_vs_highZ_Red'], gama_ia['highZ_vs_highZ_Blue']
for k in gama_ia.keys():
	gama_ia[k] = gama_ia[k][:-1]
#gama_ia['GAMA total'] = np.loadtxt('../HH_datfiles/wgp_z2_r')
#gama_ia['GAMA red'][:, 3] = np.loadtxt('../Vanilla+SDSSdatfiles/wgp_z2_r',usecols=3)
#gama_ia['GAMA blue'][:, 3] = np.loadtxt('../Vanilla+SDSSdatfiles/wgp_z2_b',usecols=3)
gama_gg = read_files('wgg*all.dat', idnot='pair*ic', dire='playground/GAMA')
for k in gama_gg.keys():
	gama_gg[k] = gama_gg[k][:-1]

#gama_gg = merge_two_dicts(gama_gg, read_files('wgg*f12*.dat', idnot='pair*ic*pz', dire='playground/GAMA'))
#gama_gg['wgg_ps_f12_rr.dat'][:, 2] = ascii.read('../OUTPUTS_GAMA_f0oldrand/wgg_red_unwindowed.dat')['wgg_jackknife_err']
#gama_gg['wgg_ps_f12_bb.dat'][:, 2] = ascii.read('../OUTPUTS_GAMA_f0oldrand/wgg_blue_unwindowed.dat')['wgg_jackknife_err']
#gama_gg['wgg_ps.dat'][:, 2] = ascii.read('../OUTPUTS_GAMA_f0oldrand/wgg_unwindowed.dat')['wgg_jackknife_err']

for k in gama_gg.keys():
	if pz == '004' and '009' in k:
		del gama_gg[k]
	elif pz == '009' and '004' in k:
		del gama_gg[k]
if pz == '004': gama_zphot_label = r'$\mathcal{N}\,\left(z_{\rm{spec.}}, 0.004\right)$'
if pz == '009': gama_zphot_label = r'$\mathcal{N}\,\left(z_{\rm{spec.}}, 0.009\right)$'

if sys.argv[1] == 'qz':
	ax1.set_title('PAUS W3 -- best 50% Qz')
	pau_ia = read_files('wgp_*dat', dire='OUTPUTS_PAUS_KSBqzselected_fg1fg2_PW1/', asc=1)
	pau_gg = read_files('wgg_ps', idnot='pair*ic', dire='playground/qzselect')
elif sys.argv[1] == 'tot': 
	ax1.set_title('PAUS W3 -- all galaxies')
	pau_ia = read_files('wgp_*dat', dire='OUTPUTS_PAUS_KSBtotal_fg1fg2_PW1/', asc=1)
	pau_gg = read_files('wgg_ps', idnot='pair*ic', dire='playground')

#pau_gg['wgg_ps_f0_rr.dat'][:, 2] = ascii.read('OUTPUTS_PAUS_KSBtotal_fg1fg2_AS/wgg_red_unwindowed.dat')['wgg_jackknife_err']
#pau_gg['wgg_ps_f0_bb.dat'][:, 2] = ascii.read('OUTPUTS_PAUS_KSBtotal_fg1fg2_AS/wgg_blue_unwindowed.dat')['wgg_jackknife_err']
#if len(pau_gg['wgg_ps_f0.dat']) != len(pau_gg['wgg_ps_f0_rr.dat']):
#	r = pau_gg['wgg_ps_f0_rr.dat'][:, 0]
#	w = np.interp(r, pau_gg['wgg_ps_f0.dat'][:, 0], pau_gg['wgg_ps_f0.dat'][:, 1])
#	e = ascii.read('OUTPUTS_PAUS_KSBtotal_fg1fg2_AS/wgg_unwindowed.dat')['wgg_jackknife_err']
#	pau_gg['wgg_ps_f0.dat'] = np.column_stack((r, w, e))
#else:
#	pau_gg['wgg_ps_f0.dat'][:, 2] = ascii.read('OUTPUTS_PAUS_KSBtotal_fg1fg2_AS/wgg_unwindowed.dat')['wgg_jackknife_err']

no_total = 1
rpp = 0.8
s3 = splitN(3, 0.03)
for k in np.sort(pau_ia.keys()):
	if no_total:
		if 'wgp_un' in k: continue
	r, w, e = pau_ia[k]['rnom'], pau_ia[k]['wgplus'], pau_ia[k]['wgplus_jackknife_err']
	if 'red' in k:
		c = 'r'
		lab = 'PAUS W3 red'
		s = s3[0]
	elif 'blue' in k:
		c = 'b'
		lab = 'PAUS W3 blue'
		s = s3[2]
	else:
		c = 'goldenrod'
		lab = 'PAUS W3 total'
		s = s3[1]
	ax1.errorbar(r*s, r**rpp*w, r**rpp*e, fmt='.-', lw=1.3, elinewidth=1.5, c=c, label=lab, capsize=0)
	axins.errorbar(r*s, r**rpp*w, r**rpp*e, fmt='.', lw=1.3, elinewidth=1.5, c=c, label=lab, capsize=0)
for k in np.sort(gama_ia.keys()):
	if no_total:
		if 'total' in k: continue
	r, w, e = gama_ia[k].T[[0,1,3]]
	if 'red' in k or 'GR' in k:
		c = 'r'
		s = s3[0]
	elif 'blue' in k or 'GB' in k:
		c = 'b'
		s = s3[2]
	else:
		c = 'goldenrod'
		s = s3[1]
	if '~' in k or 'pz' in k:
		ms = 12
		cs = 4
		s *= 1.01
	else:
		ms = 7
		cs = 1.3
		s *= 0.99
	lab = k
	if 'GB' in k or 'GR' in k: lab = gama_zphot_label
	ax1.errorbar(r*s, r**rpp*w, r**rpp*e, fmt='o:', lw=0.6, ms=ms, c=c, mfc='none', label=lab, capsize=cs)
	axins.errorbar(r*s, r**rpp*w, r**rpp*e, fmt='o', lw=0.6, ms=ms, c=c, mfc='none', label=lab, capsize=cs)
for k in np.sort(pau_gg.keys()):
	if no_total:
		if 'rr' not in k and 'bb' not in k: continue
	r, w, e = pau_gg[k].T
	if 'rr' in k:
		c = 'r'
		s = s3[0]
	elif 'bb' in k:
		c = 'b'
		s = s3[2]
	else:
		c = 'goldenrod'
		s = s3[1]
	ax2.errorbar(r*s, w, e, fmt='.-', lw=1.3, elinewidth=1.5, c=c, capsize=0)
for k in np.sort(gama_gg.keys()):
	if no_total:
		if ('rr' not in k and 'bb' not in k and 'red' not in k and 'blue' not in k): continue
	try:
		r, w, e, jkerr = gama_gg[k].T
	except:
		r, w, e = gama_gg[k].T
		jkerr = e
		print k, 'no jackknife'
	if 'rr' in k or 'red' in k or 'GR' in k:
		c = 'r'
		s = s3[0]
	elif 'bb' in k or 'blue' in k or 'GB' in k:
		c = 'b'
		s = s3[2]
	else:
		c = 'goldenrod'
		s = s3[1]
	fmt = 'o:'
	if '~' in k or 'pz' in k:
		ms = 12
		cs = 4
		s *= 1.01
	else:
		ms = 7
		cs = 1.3
		s *= 0.99
	ax2.errorbar(r*s, w, jkerr, fmt=fmt, lw=0.6, ms=ms, c=c, mfc='none', capsize=cs)

x1, x2, y1, y2 = 0.12, 5.3, -0.16, 0.15
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticks([0.2, 1, 3])
axins.set_xticklabels(['0.2','1.0','3.0'], fontdict={'fontsize':11})
#axins.set_yticklabels([])
mark_inset(ax1, axins, loc1=1, loc2=3, fc='none', ec='k', alpha=0.4, lw=0.5)
h, l = ax1.get_legend_handles_labels()
ax2.legend(h, l, ncol=3, loc='best', fontsize=12)
ax1.set_ylim(None, 1.4)
#ax2.set_ylim(0.3, 1e3)
ax1.set_ylabel(r'$r_{p}^{0.8}w_{\rm{g+}}\,[h^{-1}\rm{Mpc}]^{1.8}$')
ax2.set_ylabel(r'$w_{\rm{gg}}\,[h^{-1}\rm{Mpc}]$')
ax2.set_xlabel(r'$r_{p}\,\,[h^{-1}\rm{Mpc}]$')
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.show()
plt.savefig('PAUvsGAMA_figure_%s%s.pdf'%(sys.argv[1], sys.argv[2]), bbox_inches='tight')




