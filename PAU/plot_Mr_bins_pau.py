# coding: utf-8
from functions import *
f, (a1,a2,a3) = plt.subplots(3, sharex=True, figsize=(8,9))
flist = np.sort(glob('OUTPUTS_PAUS_KSB*Mr*/wgp*dat'))
spl = splitN(3, 0.02)
rpp = 0.8
Mrdict = dict(zip(flist, ['$M_{r} > -19.69$', '$M_{r} > -20.90$', '$M_{r} > -20.01$',
    '$M_{r} < -21.00$', '$M_{r} < -21.78$', '$M_{r} < -21.26$',
    '$-21.00 > M_{r} > -19.69$', '$-21.78 > M_{r} > -20.90$', '$-21.26 > M_{r} > -20.01$']))
flist1 = [i for i in flist if 'high' in i]+[i for i in flist if 'mid' in i]+[i for i in flist if 'low' in i]
for df in flist1:
    lab = Mrdict[df]
    if 'wgp_un' in df:
        a = a1
        c = 'goldenrod'
    elif 'wgp_red' in df:
        a = a2
        c = 'red'
    elif 'wgp_blue' in df:
        a = a3
        c = 'blue'
    if 'lowMr' in df:
        lst = '-'
        s = spl[0]
    elif 'midMr' in df:
        lst = '--'
        s = spl[1]
    elif 'highMr' in df:
        lst = ':'
        s = spl[2]
    dat = ascii.read(df)
    r, w, e = dat['rnom'], dat['wgplus'], dat['wgplus_jackknife_err']
    eb = a.errorbar(r*s, r**rpp*w, r**rpp*e, ls=lst, c=c, capsize=2, label=lab)
    eb[-1][0].set_linestyle(lst)
for a in (a1,a2,a3):
    a.legend(loc='best', ncol=3)
    a.set_xscale('log')
    a.set_ylabel(r'$r_{p}^{%s}w_{\rm{g+}}(r_{p})$'%rpp)
    a.axhline(0, c='k', ls=':', lw=0.8)
a3.set_xlabel(r'$r_{p}$')
plt.tight_layout()
