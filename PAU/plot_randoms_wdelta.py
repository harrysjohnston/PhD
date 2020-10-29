# coding: utf-8
from functions import *
from clone_randoms import *
if 'gama' in sys.argv[1] or 'phot_rand' in sys.argv[1] or 'spec_rand' in sys.argv[1] or 'zphot_' in sys.argv[1]:
	cat = fopen('../SMLambdarApMatchedPhotom.fits')
	zph = 'zphot_2'
	zcol = 'Z_TONRY'
	lab = 'GAMA'
	xlim = -0.02, 0.6
else:
	cat = fopen('PAUS_KSB.fits')
	zph = 'bcnz_zb'
	zcol = 'ZBEST'
	lab = 'PAUS W3'
	xlim = -0.05, 1.25
if len(sys.argv) > 2:
	MODE = 'multi'
	diagfiles = sys.argv[1:]
else:
	MODE = 'single'
	diagfiles = [sys.argv[1]]

f, ax = plt.subplots(2, sharex=True, figsize=(7,9))
cmaps = iter([plt.cm.Reds, plt.cm.Blues, plt.cm.Greens, plt.cm.Oranges, plt.cm.Purples])
hline = ax[1].axhline(1, c='k', ls='-')
cm = ['C%s'%i for i in range(len(diagfiles))]
#ax[1].set_yscale('log')

h, l = [], []
for j, df in enumerate(diagfiles):
	diag = h5py.File(df,'r')
	nr = diag['n_r'][:]
	dmid = diag['d_mid'][:]
	dmid1 = np.linspace(dmid.min(), dmid.max(), 50)
	delt = diag['Delta_d'][:]
	#cm = next(cmaps)(np.linspace(0.4,0.9,len(nr)))

	if 'zph' in df and 'Unwindowed' not in df:
		df_label = 'zph-windowed'
	elif 'zph' in df and 'Unwindowed' in df:
		df_label = 'zph-unwindowed'
	elif 'Unwindowed' in df:
		df_label = 'unwindowed'
	else:
		df_label = 'windowed'
	if 'phot_rand' in df:
		df_label = '$n(z_{phot})$'
	if 'spec_rand' in df:
		df_label = '$n(z_{spec})$'
	if 'zphot_1' in df:
		df_label = '$n(z_{phot},1)$'
	if 'zphot_2' in df:
		df_label = '$n(z_{phot},2)$'

	if df == diagfiles[0]:
		hist = ax[0].hist(cat[zph], bins=len(dmid)//5, color=hline.get_c(), histtype='step', lw=1, density=1, label=lab+r' $z_{\rm{phot.}}$')
		hist2 = ax[0].hist(cat[zcol][(cat[zcol]>0)&np.isfinite(cat[zcol])], bins=len(dmid)//13, color=hline.get_c(), histtype='step', ls='--', lw=1, density=1, label=lab+r' $z_{\rm{spec.}}$')

	#for i in range(len(nr)):
	#	try:
	#		norm = np.trapz(nr[i], x=get_z(dmid))
	#		#if i in range(len(nr))[::2]:
	#		ax[0].plot(get_z(dmid), nr[i]/norm, c=cm[i], alpha=0.8, lw=0.9)
	#	except:
	#		pass
	#	#if i > 12:
	#		#z_int = np.linspace(hist[1].min(), hist[1].max(), len(hist[1])-1)
	#		#delt_int = np.interp(z_int, get_z(dmid), delt[i])
	#		#ax[1].plot(z_int, delt_int, c=cm[i], alpha=0.5)
	norm = np.trapz(nr[-1], x=get_z(dmid))
	#if i in range(len(nr))[::2]:
	#if 'master' not in df:
		#ax[1].plot(get_z(dmid), smooth(delt[-1], 5), c=cm[j], alpha=0.7, lw=1.6)
	#else:
	ax[1].plot(get_z(dmid1), np.interp(dmid1, dmid, delt[-1]), c=cm[j], lw=1.6)
	h1, = ax[0].plot(get_z(dmid), nr[-1]/norm, c=cm[j], lw=1.6)

	#h, l = ax[0].get_legend_handles_labels()
	#h1 = tuple([new_handle(ls='-',c=c,alpha=0.7) for c in cm])
	#h1 = new_handle(ls='-',c=cm[j],alpha=0.7)
	h += [h1]
	#l += ['%s\niterations'%df_label]
	l += [df_label]

	diag.close()

h = [hist[-1][0], hist2[-1][0]] + h
l = [lab+r' $z_{\rm{phot.}}$',lab+r' $z_{\rm{spec.}}$'] + l
ax[0].legend(h, l, loc='best', fontsize=11, frameon=0)
ax[1].set_xlabel('$z$')
ax[1].set_ylim(0, 2.5)
ax[1].set_xlim(xlim)
ax[0].set_ylabel('PDF')
ax[1].set_ylabel(r'$\Delta(z) \, / \, N_{\rm{clone}}$')
plt.tight_layout()
plt.show()

plt.savefig('%s_randoms_wdelta.pdf'%lab.replace(' ',''), bbox_inches='tight')


