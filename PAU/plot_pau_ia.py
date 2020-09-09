# coding: utf-8
from functions import *
for x in ['total_fg1fg2_AS', 'total_fg1fg2_PW1','total_fg1_AS', 'total_fg1_PW1', 'total_fg2_AS', 'total_fg2_PW1', 'total_f0_AS', 'total_f0_PW1']:
    ipython.magic('run plot_pau_vs_gama.py %s'%x)
    x1 = x.replace('total','').replace('g','').replace('f1f','f1').replace('_AS','').replace('_PW1','')
    rr = np.loadtxt('playground/wgp_ps%s_rr.dat'%x1)
    ax.errorbar(rr[:,0], rr[:,0]**0.8*rr[:,1], rr[:,0]**0.8*rr[:,3], fmt='r:', capsize=3)
    ax1.errorbar(rr[:,0], rr[:,0]**0.8*rr[:,2], rr[:,0]**0.8*rr[:,3], fmt='r:', capsize=3)
    bb = np.loadtxt('playground/wgp_ps%s_bb.dat'%x1)
    ax.errorbar(bb[:,0], bb[:,0]**0.8*bb[:,1], bb[:,0]**0.8*rr[:,3], fmt='b:', capsize=3)
    ax1.errorbar(bb[:,0], bb[:,0]**0.8*bb[:,2], bb[:,0]**0.8*rr[:,3], fmt='b:', capsize=3)
    un = np.loadtxt('playground/wgp_ps%s.dat'%x1)
    ax.errorbar(un[:,0], un[:,0]**0.8*un[:,1], un[:,0]**0.8*un[:,3], fmt=':', c='goldenrod', capsize=3)
    ax1.errorbar(un[:,0], un[:,0]**0.8*un[:,2], un[:,0]**0.8*un[:,3], fmt=':', c='goldenrod', capsize=3)

plt.show()
