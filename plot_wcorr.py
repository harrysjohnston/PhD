#!/user/bin/env python
from __future__ import print_function, division
from astropy.io import fits
import numpy as np
from astropy.io import ascii
from os.path import join, isdir, basename, normpath
from os import listdir, mkdir
import os
import argparse
import csv
from astropy import cosmology
from astropy.cosmology import Planck13
import matplotlib
import matplotlib.pyplot as plt
import sys
import scipy.integrate as scint
import scipy.stats as stat
from scipy.stats import chi2

def plot(path):
	listDir = listdir(path)
	cut = [i.startswith('largePi') for i in listDir]
	cut = np.invert(cut)
	listDir = listDir[cut]
	listDir.sort()
	listDir = np.array(listDir)[[1,0,3,2]]
	dataArr = [np.loadtxt(join(path,i)) for i in listDir]
	# [r_p, wgplus, wgcross, wgerr] for each of 4 datasets
	r_p = dataArr[0][:,0]

	f, axarr = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(15,10))
	f.subplots_adjust(hspace=0, wspace=0)
	axarr[0,0].errorbar(r_p, dataArr[0][1], yerr=dataArr[0][3],
						elinewidth=2, color='r', capsize=0,
						label='w(g+)')
	axarr[0,0].errorbar(r_p, dataArr[0][2], yerr=dataArr[0][3],
						elinewidth=2, color='g', capsize=0,
						label='w(gx)', alpha=0.5)
	axarr[0,1].errorbar(r_p, dataArr[1][1], yerr=dataArr[1][3],
						elinewidth=2, color='b', capsize=0,
						label='w(g+)')
	axarr[0,1].errorbar(r_p, dataArr[1][2], yerr=dataArr[1][3],
						elinewidth=2, color='g', capsize=0,
						label='w(gx)', alpha=0.5)
	axarr[1,0].errorbar(r_p, dataArr[2][1], yerr=dataArr[2][3],
						elinewidth=2, color='r', capsize=0,
						label='w(g+)')
	axarr[1,0].errorbar(r_p, dataArr[2][2], yerr=dataArr[2][3],
						elinewidth=2, color='g', capsize=0,
						label='w(gx)', alpha=0.5)
	axarr[1,1].errorbar(r_p, dataArr[3][1], yerr=dataArr[3][3],
						elinewidth=2, color='b', capsize=0,
						label='w(g+)')
	axarr[1,1].errorbar(r_p, dataArr[3][2], yerr=dataArr[3][3],
						elinewidth=2, color='g', capsize=0,
						label='w(gx)', alpha=0.5)
	arr_ind = [(0,0), (0,1), (1,0), (1,1)]
	for i, ind in enumerate(arr_ind):
		a = axarr[ind]
		a.set_xscale('log')
		a.set_xlim(0.25,70)
		a.set_ylim(-0.5,0.4)
		a.plot(x, [0]*len(x), lw=2, ls='--', color='c')
		# a.set_xlabel('Comoving transverse separation (Mpc/h)')
		# a.set_ylabel('Correlations')
		a.set_title('%s'%listDir[i], fontsize=12)
		a.legend(loc='upper right')
		a.grid()
	plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
	plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
	# plt.setp([a.get_xlabels() for a in axarr[:, 1]], visible=False)
	axarr[1,0].set_xlabel('Comoving transverse separation (Mpc/h)')
	axarr[1,0].set_ylabel('Correlations')
	# DO CUTS!!
	ZC = np.loadtxt('%s/../ZC_cuts'%path, delimiter=',')
	axarr[0,0].set_title('Cuts: z%s, c%s'%(ZC[0],ZC[1]))
	plt.show()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
	parser.add_argument(
		'path',
		help='full path of directory containing random-subtracted correlation data')
	args = parser.parse_args()

	plot(args.path)






