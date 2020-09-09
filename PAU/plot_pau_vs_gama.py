# coding: utf-8
from functions import *

rpp = 0.8

f, (ax, ax1) = plt.subplots(2, figsize=(7,7))
flist = glob('OUTPUTS_PAUS_KSB%s/wgp*dat'%sys.argv[1])# + \
        #glob('../OUTPUTS_GAMA%s/wgp*dat'%sys.argv[1])
flist = np.sort(flist)
flist = [fl for fl in flist if 'ps' not in fl and 'wgp_un' not in fl]

gama_r = np.loadtxt('../Wcorr_fullGAMA_rb_c0.66/to_plot/highZ_vs_highZ_Red')
gama_b = np.loadtxt('../Wcorr_fullGAMA_rb_c0.66/to_plot/highZ_vs_highZ_Blue')
r, w, wx, e, xe = gama_b[:, [0, 1, 4, 3, 6]].T
w *= 0.934
e *= 0.934
ax.errorbar(r, r**rpp*w, r**rpp*e, fmt='P:', c='b', ms=8, mfc='none', capsize=4, label='GAMA/blue', alpha=0.6)
ax1.errorbar(r, r**rpp*wx, r**rpp*xe, fmt='X:', c='b', ms=8, mfc='none', capsize=4, alpha=0.6)
r, w, wx, e, xe = gama_r[:, [0, 1, 4, 3, 6]].T
w *= 0.934
e *= 0.934
ax.errorbar(r, r**rpp*w, r**rpp*e, fmt='P:', c='r', ms=8, mfc='none', capsize=4, label='GAMA/red', alpha=0.6)
ax1.errorbar(r, r**rpp*wx, r**rpp*xe, fmt='X:', c='r', ms=8, mfc='none', capsize=4, alpha=0.6)

paudat = {}
for i, datf in enumerate(flist):
	lab = datf
	for p in ['wgp_','OUTPUTS_','_KSBtotal','.dat','OUTPUTS_','../','_fg1fg2','_PW1','_unwindowed']:
		lab = lab.replace(p,'')
		
	if 'red' in datf: c = 'r'
	elif 'blue' in datf: c = 'b'
	else: c = 'goldenrod'
	if 'GAMA' in datf: ls = '--'
	else: ls = '-'
	si = splitN(len(flist), 0.03)[i]
	
	dat = ascii.read(datf)
	r = dat['rnom'].data
	wgp, wgx = dat['wgplus'].data, dat['wgcross'].data
	if 'GAMA' in datf:
		wgp *= 0.934
		wgx *= 0.934
	try:
		wgperr, wgxerr = dat['wgplus_jackknife_err'].data, dat['wgcross_jackknife_err'].data
		#wgpjm, wgxjm = dat['wgplus_jackknife_mean'].data, dat['wgcross_jackknife_mean'].data
	except:
		wgperr = wgxerr = np.zeros_like(wgp)

	eb = ax.errorbar(r*si, r**rpp*wgp, r**rpp*wgperr, fmt='+-', ms=8, c=c, ls=ls, capsize=2, label=lab, alpha=0.8)
	ebx = ax1.errorbar(r*si, r**rpp*wgx, r**rpp*wgxerr, fmt='x-', ms=6, ls=eb[0].get_linestyle(), c=eb[0].get_color(), label=None, capsize=2, alpha=0.8)
#    ebx[-1][0].set_linestyle(':')
	#ax.plot(r*si, r**rpp*wgpjm, 's', c=eb[0].get_color(), markerfacecolor='none', ms=10, label=None)
	paudat[lab] = dat

for a in [ax,ax1]:
    a.set_xscale('log')
    a.axhline(0, c='k', ls=':', lw=1, zorder=10)

ax.legend(loc='best')
ax.set_ylabel(r'$r_{p}^{%s}w_{\rm{g+}}(r_{p})\quad[h^{-1}\rm{Mpc}]^{%s}$'%(rpp, 1+rpp),fontsize=14)
ax1.set_ylabel(r'$r_{p}^{%s}w_{\rm{gx}}(r_{p})\quad[h^{-1}\rm{Mpc}]^{%s}$'%(rpp, 1+rpp), fontsize=14)
ax1.set_xlabel(r'$r_{p}\,[h^{-1}\rm{Mpc}]$')
plt.tight_layout()
plt.show()



