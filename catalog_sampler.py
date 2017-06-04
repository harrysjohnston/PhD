#!/user/bin/env python
from __future__ import print_function, division
from astropy.io import fits
import numpy as np
from astropy.io import ascii
from os.path import join, isdir, basename, normpath
from os import listdir, mkdir
import copy
import os
import argparse
import csv
from astropy import cosmology
from astropy.cosmology import Planck13
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import scipy.integrate as scint
import scipy.stats as stat
from scipy.stats import chi2
try:
	import healpy as hp
except ImportError:
	print('healpy import error')
import gc
import astropy.stats as astat
import jackknife3d as jk3d

class RealCatalogue:

	def __init__(self, path, DEI, mc, SDSS): # ADD MORE self.SPECS HERE; FEWER ARGS FOR FNS!
		"""""
		read-in catalogue

		"""""
		self.path = path
		KSBheads = ['RA_1_1','DEC_1_1','Z_TONRY','e1c','e2c','RankBCG_1','logmstar','pgm','absmag_g_1','absmag_i_1','col3']
		DEIheads = ['RA_GAMA','DEC_GAMA','Z_TONRY','e1','e2','RankBCG','logmstar','pgm','absmag_g_1','absmag_i_1','MASK']
		SDSSheads = ['ra','dec','z','e1','e2','sigma_gamma','rf_g-r']
		if DEI:
			self.headers = DEIheads
			self.DEI = 1
			self.SDSS = 0
		elif SDSS:
			self.headers = SDSSheads
			self.DEI = 0
			self.SDSS = 1
		else:
			self.headers = KSBheads
			self.DEI = 0
			self.SDSS = 0
		hdulist = fits.open(path)
		data = hdulist[1].data
		print('SELECTING R_MAG < %s'%mc)
		data = data[data['%s'%(['absmag_r_1', 'M_r'][self.SDSS])]<mc]
		self.data = data
		del data
		gc.collect()
		self.columns = hdulist[1].columns.names
		# ascii (.asc) file IDs
		self.labels = ['highZ_Red', 'highZ_Blue', 'lowZ_Red', 'lowZ_Blue',
						'highZ_Red_UnMasked', 'highZ_Blue_UnMasked', 'lowZ_Red_UnMasked', 'lowZ_Blue_UnMasked',
						'highZ', 'lowZ']
		ext_labels = ['%s_largePi'%i for i in self.labels[:4]]
		ext_labels = self.labels[:4]+ext_labels
		ext_labels.sort()
		self.ext_labels = ext_labels
		# bodies of wcorr output file IDs
		self.wcorrLabels = ['highZ_vs_highZ_Red', 'highZ_vs_highZ_Blue', 'lowZ_vs_lowZ_Red', 'lowZ_vs_lowZ_Blue']
		# Npatch = {3.0:96,4.0:54,9.0:24}
		# self.Npatch = Npatch[psize]

		# MEASURE WG+ SIGNALS

	def cut_data(self, pgm_, z_, colour_, lmstar_, LRG, LRG1, LRG2, LRG3, LRG4, LRG5, LRG6, BCGdens, BCGshap, *bitmask_):
		"""""
		cut catalogue according to bitmasks, 
		PGM, & into subsamples

		"""""
		for head in self.headers:
			assert head in self.columns, "'%s' not in columns, see headers: %s"%(head,self.columns)

		if z_!=None:
			self.zstr = '%.f'%z_
		else:
			self.zstr = None
		if colour_!=None:
			self.cstr = '%.f'%colour_
		elif LRG:
			self.cstr = 'LRGs'
		else:
			self.cstr = None

		self.pre_count = len(self.data)
		z = self.data[self.headers[2]]
		self.pre_z = z
		colour = self.data[self.headers[-3]] - self.data[self.headers[-2]]
		total_bitmasks = self.data[self.headers[-1]]
		if not self.SDSS:
			logmstar = self.data[self.headers[6]]+np.log10(self.data['fluxscale'])

		if self.SDSS:
			colour = self.data[self.headers[-1]]
			total_bitmasks = np.zeros_like(colour)
			logmstar = np.ones_like(colour)

		if LRG:
			g_s,r_s,i_s = self.data['dered_g'],self.data['dered_r'],self.data['dered_i']
			ext_r = self.data['extinction_r']
			rpetro = self.data['petroMag_r']
			cpar = 0.7*(g_s-r_s)+1.2*(r_s-i_s-0.18)
			cperp = (r_s-i_s)-(g_s-r_s)/4.-0.18
			dperp = (r_s-i_s)-(g_s-r_s)/8.

			#cutI
			# print('CUTTING by r_s NOT r_petro-ext_r...!')
			rscut = (16.<=(rpetro-ext_r))&((rpetro-ext_r)<=19.2)
			# rscut = (16.<=(r_s))&((r_s)<=19.2)
			rscpar = (r_s)<(13.1+(cpar/0.3))
			abscperp = abs(cperp)<0.2
			#cut2
			rscut2 = (16.<=(r_s))&((r_s)<=19.5)
			cperpcut = cperp>(0.45-(g_s-r_s/6.))
			gsrscut = (g_s-r_s)>(1.3+0.25*(r_s-i_s))
			#LOWZcut
			rscutL = (16.<=(r_s))&((r_s)<=19.6)
			rscparL = (r_s)<(13.5+(cpar/0.3))

			cut1 = abscperp
			if LRG1:
				cut1 &= rscut
			if LRG2:
				cut1 &= rscpar
			cut2 = cperpcut
			if LRG3:
				cut2 &= rscut2
			if LRG4:
				cut2 &= gsrscut
			LOWZcut = abscperp
			if LRG5:
				LOWZcut &= rscutL
			if LRG6:
				LOWZcut &= rscparL
			# print('selecting LOWZ exact...!!')
			LRGcut = LOWZcut|cut1|cut2

		if self.DEI|self.SDSS:
			pgm = np.ones_like(colour)
		else:
			pgm = self.data['pgm']
		pgm_cut = np.array((pgm > pgm_))
		zeroPgm_cut = np.array((pgm != 0))
		print('pgm cut [unique]: \t', np.unique(pgm_cut))

		# define colour, redshift, bitmask & BCG cuts
		if colour_ != None:
			red_cut = np.array((colour > colour_)) # larger (B-V) <-> 'redder' colour
			blue_cut = ~red_cut
			print('c cut [unique]: \t', colour_, np.unique(red_cut))
		else:
			red_cut = np.array([True]*len(self.data))
			blue_cut = red_cut
			print('Red catalog == Blue catalog')
		if z_ != None:
			z_cut = np.array((z > z_)) # HIGH-Z
			z_cut_r = ~z_cut # LOW-Z
		else:
			z_cut = np.array([True]*len(self.data))
			z_cut_r = z_cut
			print('highZ catalog == lowZ catalog')
		print('z cut [unique]: \t', z_, np.unique(z_cut))
		if lmstar_ != None:
			lmstar_cut = np.array((logmstar>lmstar_[0])&(logmstar<lmstar_[1]))
		elif LRG:
			lmstar_cut = LRGcut
			print('cutting for (~%s) LRGs...!'%sum(LRGcut))
		else:
			lmstar_cut = np.array([True]*len(self.data))
		print('lmstar cut [unique]: \t', lmstar_, np.unique(lmstar_cut))

		bitmask_cut = [True]*len(total_bitmasks)

		print("CUTTING MASK!=0...")
		bitmask_cut = np.where(total_bitmasks==0,True,False)

		if bitmask_[0] != None:
			bitmask_ = bitmask_[0]
			for i in range(len(bitmask_)):
				# construct bitmask cut
				bitmask_cut &= np.where(bitmask_[i] & total_bitmasks == bitmask_[i], False, True)

		assert len(bitmask_cut) == len(total_bitmasks), "bitmask testing broken"
		bitmask_cut = np.array(bitmask_cut)
		print('bitmask cut [unique]: \t', np.unique(bitmask_cut))

		BCGcut = np.where(self.data[self.headers[5]]==1,True,False)
		BCG_dc = [True]*len(self.data)
		BCG_sc = BCG_dc
		BCGargs = [0,0]
		if BCGdens:
			BCG_dc &= BCGcut
			BCGargs[1] = 1
		if BCGshap:
			BCG_sc &= BCGcut
			BCGargs[0] = 1

		# apply cuts
		self.highz_R = self.data[(z_cut&red_cut&bitmask_cut&BCG_sc&lmstar_cut)]
		self.highz_B = self.data[(z_cut&blue_cut&bitmask_cut&BCG_sc&lmstar_cut)]
		self.lowz_R = self.data[(z_cut_r&red_cut&bitmask_cut&BCG_sc&lmstar_cut)]
		self.lowz_B = self.data[(z_cut_r&blue_cut&bitmask_cut&BCG_sc&lmstar_cut)]

		self.highz_R_UnM = self.data[(z_cut&red_cut&BCG_sc&lmstar_cut)]
		self.highz_B_UnM = self.data[(z_cut&blue_cut&BCG_sc&lmstar_cut)]
		self.lowz_R_UnM = self.data[(z_cut_r&red_cut&BCG_sc&lmstar_cut)]
		self.lowz_B_UnM = self.data[(z_cut_r&blue_cut&BCG_sc&lmstar_cut)]

		self.highz = self.data[(z_cut&BCG_dc&lmstar_cut)]
		self.lowz = self.data[(z_cut_r&BCG_dc&lmstar_cut)]

		self.samplecounts = []

		# save cuts for later use
		self.Rmags = []
		for sample in [self.highz_R,self.highz_B,self.lowz_R,self.lowz_B]:
			self.Rmags.append(np.mean(sample['%s'% (['absmag_r_1', 'M_r'] [self.SDSS]) ]))
		self.zcut = z_cut
		self.zcut_r = z_cut_r
		self.redcut = red_cut
		self.bluecut = blue_cut
		self.lmstarcut = lmstar_cut
		self.bitmaskcut = bitmask_cut
		self.pgmcut = pgm_cut
		self.BCG_dc = BCG_dc
		self.BCG_sc = BCG_sc
		self.BCGargs = BCGargs

		self.hzr_cut = (self.zcut&self.redcut&self.bitmaskcut&self.BCG_sc&self.lmstarcut)
		self.hzb_cut = (self.zcut&self.bluecut&self.bitmaskcut&self.BCG_sc&self.lmstarcut)
		self.lzr_cut = (self.zcut_r&self.redcut&self.bitmaskcut&self.BCG_sc&self.lmstarcut)
		self.lzb_cut = (self.zcut_r&self.bluecut&self.bitmaskcut&self.BCG_sc&self.lmstarcut)
		self.hzrU_cut = (self.zcut&self.redcut&self.BCG_sc&self.lmstarcut)
		self.hzbU_cut = (self.zcut&self.bluecut&self.BCG_sc&self.lmstarcut)
		self.lzrU_cut = (self.zcut_r&self.redcut&self.BCG_sc&self.lmstarcut)
		self.lzbU_cut = (self.zcut_r&self.bluecut&self.BCG_sc&self.lmstarcut)
		self.hzU_cut = (self.zcut&self.BCG_sc&self.lmstarcut)
		self.lzU_cut = (self.zcut_r&self.BCG_sc&self.lmstarcut)
		self.samplecuts = np.array([self.hzr_cut,self.hzb_cut,self.lzr_cut,self.lzb_cut,
							self.hzrU_cut,self.hzbU_cut,self.lzrU_cut,self.lzbU_cut,
							self.hzU_cut,self.lzU_cut])

		return None

	def cut_columns(self, subsample, h):
		"""""
		take subsample data 
		& isolate columns for wcorr

		"""""
		table = subsample
		RA = np.deg2rad(table[self.headers[0]])
		DEC = np.deg2rad(table[self.headers[1]])
		Z = table[self.headers[2]]
		if self.DEI|self.SDSS:
			pgm = np.ones_like(Z)
		else:
			pgm = table['pgm']
		e1 = table[self.headers[3]]/pgm
		e2 = table[self.headers[4]]/pgm
		# e2 *= -1 # for RA increasing leftward, c.f. x-axis increasing rightward ???
		e_weight = np.where(pgm<0.1,0,pgm)
		if self.DEI:
			e_weight = np.where(table['flag_DEIMOS']=='0000',1,0)
		if self.SDSS:
			e_weight = np.where( ((abs(e1)>9.9)|(abs(e2)>9.9)) ,0,1)
		if len(e1)>1000:
			e1m,e2m = e1[np.array(e_weight,dtype=bool)],e2[np.array(e_weight,dtype=bool)]
			e1m,e2m = np.mean(e1m),np.mean(e2m)
			e1 -= e1m
			e2 -= e2m

		e1,e2,e_weight = map(lambda x: np.nan_to_num(x), [e1,e2,e_weight])

		# random re-shuffle test - density-shape corr should now ~ 0
		# e12 = list(zip(e1,e2))
		# np.random.shuffle(e12)
		# e12 = np.array(e12)
		# e1 = e12[:,0]
		# e2 = e12[:,1]
		#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   

		if 'chi_comov_TONRY' in table.columns.names:
			comov = table['chi_comov_TONRY']
		else:
			comov = Planck13.comoving_distance(Z)
		comov *= h
		new_table = np.column_stack((RA,DEC,comov,e1,e2,e_weight))

		self.samplecounts.append(len(new_table))

		return new_table,Z

	def save_tables(self, new_table, outfile_root_, label, z_cut, c_cut, notes):
		"""""
		save subsample tables to ascii

		"""""
		if notes != None:
			outfile_root_ += '_%s'%notes
		if (z_cut != None) & (c_cut != None):
			outfile_root = outfile_root_ + "_z" + str(z_cut) + "_c" + str(c_cut)
		elif z_cut != None:
			outfile_root = outfile_root_ + "_z" + str(z_cut)
		elif c_cut != None:
			outfile_root = outfile_root_ + "_c" + str(c_cut)
		else:
			outfile_root = outfile_root_ + "_allz"

		self.new_root = outfile_root
		if not isdir(outfile_root):
			mkdir(outfile_root)

		ascii.write(new_table, join(outfile_root, label + ".asc"), names=['#RA/rad', '#DEC/rad', '#comov_dist/Mpc/h', '#e1', '#e2', '#e_weight'])
		sample_no = "%s # objects:\t%s"%(label,len(new_table))
		return sample_no

	def save_swotfiles(self, subsample, label):
		"""""
		save samples in SWOT format (ra[deg],dec[deg],z)

		"""""
		RA = subsample[self.headers[0]]
		DEC = subsample[self.headers[1]]
		Z = subsample[self.headers[2]]
		newtable = np.column_stack((RA,DEC,Z))

		ascii.write(newtable,join(self.new_root,'swot_%s.asc'%label),names=['#ra[deg]','dec[deg]','z'],delimiter='\t')

		return Z # for random downsampling

	def save_swotpatches(self, patch, label, pnum):
		RA = patch[self.headers[0]]
		DEC = patch[self.headers[1]]
		Z = patch[self.headers[2]]
		newpatch = np.column_stack((RA,DEC,Z))

		sw_pDir = join(self.new_root,'swot_%s'%label)
		if not isdir(sw_pDir):
			mkdir(sw_pDir)

		ascii.write(newpatch,join(sw_pDir,label+'%spatch.asc'%str(pnum).zfill(2)),names=['#ra[deg]','dec[deg]','z'],delimiter='\t')

	def make_combos(self, densColours):
		# construct sets of filenames, counts, & IDs for wcorr-calls

		wcorr_ind = [[8,8,9,9], [4,5,6,7]][densColours]

		self.wcorr_combos = [
		[self.labels[wcorr_ind[0]]+'.asc', self.samplecounts[wcorr_ind[0]], self.labels[0]+'.asc', self.samplecounts[0], catalog.wcorrLabels[0]],
		[self.labels[wcorr_ind[1]]+'.asc', self.samplecounts[wcorr_ind[1]], self.labels[1]+'.asc', self.samplecounts[1], catalog.wcorrLabels[1]],
		[self.labels[wcorr_ind[2]]+'.asc', self.samplecounts[wcorr_ind[2]], self.labels[2]+'.asc', self.samplecounts[2], catalog.wcorrLabels[2]],
		[self.labels[wcorr_ind[3]]+'.asc', self.samplecounts[wcorr_ind[3]], self.labels[3]+'.asc', self.samplecounts[3], catalog.wcorrLabels[3]]
		]
		print('CHECK THIS wcorr call combinations:\n',self.wcorr_combos)
		self.samplecounts = self.samplecounts[:10]
		[print('# objects %s: \t'%self.labels[i], v) for i, v in enumerate(self.samplecounts)]
		# print('total shapes: \t%s'%np.sum(self.samplecounts[:4]))
		# print('total density: \t%s'%np.sum(self.samplecounts[4:]))

	def prep_wcorr(self, files_path, wcorr_combos, rp_bins, rp_lims, los_bins, los_lim, nproc, large_pi, out_sh):
		shell_script = [
		'#!/bin/tcsh',
		'#PBS -q compute',
		'#PBS -N %s'%out_sh,
		'#PBS -l nodes=1:ppn=%s'%nproc,
		'#PBS -l walltime=120:00:00',
		'#PBS -l mem=50gb',
		'#PBS -o %s'%join(files_path,out_sh+'.out'),
		'#PBS -e %s'%join(files_path,out_sh+'.err'),
		'',
		'module load dev_tools/nov2014/python-anaconda',
		'',
		'cd $PBS_O_WORKDIR',
		'',
		'setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/share/splinter/hj/PhD/CosmoFisherForecast/bjutils/lib/',
		'',
		'date']

		for combo in wcorr_combos: 
			# write 4x wcorr-calls to .sh script, & another 4x if largePi testing
			shell_script.append('')
			outfile = combo[4]
			shell_script.append('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s %s 0 0' %   (files_path, combo[0], combo[1], combo[2], combo[3], rp_bins, rp_lims[0], rp_lims[1], los_bins, los_lim, outfile, nproc)
			)
			shell_script.append('')
			if large_pi:
				outfile += '_largePi'
				shell_script.append('')
				shell_script.append('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s %s %s 0' %  (files_path, combo[0], combo[1], combo[2], combo[3], rp_bins, rp_lims[0], rp_lims[1], los_bins, los_lim, outfile, nproc, large_pi)
				)
			shell_script.append('')

			shell_script.append('date')

		File = join(files_path, '%s.sh'%out_sh)
		Write = open(File, 'w')
		Text = '\n'.join(shell_script)
		Write.write(str(Text))
		Write.close()

	def plot_wcorr(self, files_path, wcorrIDs, BT, JK, largePi):
		wcorrOutputs = []
		rand_wcorrOutputs = []
		for item in wcorrIDs:
			# construct filenames.dat
			wcorrOutputs.append('%s'%join(files_path, ('wcorr_' + item + '.dat')))
			rand_wcorrOutputs.append('%s'%join(files_path, ('wcorr_rand_' + item + '.dat')))
		realData = [] # clean this up by directly building the to_plot data file,
		randData = [] # starting from np.empty([bins=6,cols=8]) - poss fewer cols
		wgplus = [] #   if no BT or JK
		wgcross = []
		wgerr = []
		Pproperrs = []
		Xproperrs = []
		Pproperrs2 = []
		Xproperrs2 = []
		rand_wgplus = []
		rand_wgcross = []
		rand_wgerr = []
		easyPlotDir = join(files_path, 'to_plot')
		if not isdir(easyPlotDir):
			mkdir(easyPlotDir)
		flist = listdir(files_path)
		zcheck = np.array(['wcorr_'+i+'.dat' in flist for i in wcorrIDs])
		wcorrOutputs,rand_wcorrOutputs = np.array(wcorrOutputs),np.array(rand_wcorrOutputs)
		wcorrOutputs,rand_wcorrOutputs = wcorrOutputs[zcheck],rand_wcorrOutputs[zcheck]
		new_labels = [j for i,j in enumerate(self.labels[:4]) if 'wcorr_'+self.wcorrLabels[i]+'.dat' in listdir(files_path)]
		for i, path in enumerate(wcorrOutputs):
			realData.append(np.loadtxt(path))
			randData.append(np.loadtxt(rand_wcorrOutputs[i]))
			realData[i][:,3] -= randData[i][:,3]
			realData[i][:,4] -= randData[i][:,4] # subtracting randoms from +/x
			realErr = realData[i][:,6]
			randErr = randData[i][:,6]

			Pproperr,Xproperr = np.ones_like(realErr),np.ones_like(realErr)
			if BT:
				BTerrs = np.loadtxt(join(files_path,'BTerrs_%s'%new_labels[i]))
				if largePi:
					BTerrs = np.loadtxt(join(files_path,'BTerrs_%s_largePi'%new_labels[i]))
				Perr = BTerrs[:,0]
				Xerr = BTerrs[:,1]
				Pproperr = np.sqrt((Perr**2) + (randErr**2))
				Xproperr = np.sqrt((Xerr**2) + (randErr**2)) # propagate errors
				Pproperrs.append(Pproperr)
				Xproperrs.append(Xproperr)

			Pproperr2,Xproperr2 = np.ones_like(realErr),np.ones_like(realErr)
			if JK:
				try:
					JKerrs = np.loadtxt(join(files_path,'JKerrs_%s'%new_labels[i]))
					if largePi:
						JKerrs = np.loadtxt(join(files_path,'JKerrs_%s_largePi'%new_labels[i]))
				except IOError:
					print('IOError THROWN for JKerrs_%s.. (largePi=%s)!!'%(new_labels[i],largePi))
					JKerrs = np.ones_like(np.array([realErr,randErr]).T)
				Perr2 = JKerrs[:,0]
				Xerr2 = JKerrs[:,1]
				Pproperr2 = np.sqrt((Perr2**2) + (randErr**2))
				Xproperr2 = np.sqrt((Xerr2**2) + (randErr**2)) # propagate errors
				Pproperrs2.append(Pproperr2)
				Xproperrs2.append(Xproperr2)

			propgErrs = np.sqrt((realErr**2) + (randErr**2)) # propagate errors
			wgplus.append(realData[i][:,3])
			wgcross.append(realData[i][:,4])
			wgerr.append(propgErrs)
			rand_wgplus.append(randData[i][:,3])
			rand_wgcross.append(randData[i][:,4])
			rand_wgerr.append(randData[i][:,6])
			# save reduced data to ascii for easy plotting

			reducedData = np.column_stack((realData[0][:,0], realData[i][:,3], Pproperr, Pproperr2, realData[i][:,4], Xproperr, Xproperr2, propgErrs)) # = [r_p, wgplus, BTPerr, JKPerr, wgcross, BTXerr, JKXerr, analyticErrs]
			ascii.write(reducedData, join(easyPlotDir, basename(normpath(path))[6:-4]), delimiter='\t', names=['r_p', 'wg+', 'BT+err', 'JK+err', 'wgx', 'BTxerr', 'JKxerr', 'analyticerrs'])

		return wcorrOutputs

		# COMPUTE SAMPLE COVARIANCE & ERRORS

	def chi2(self, path2data, expec, dof, covartype):
		filesList = np.array(listdir(path2data))
		wcorrData = np.array(['_vs_' in x for x in filesList])
		covarData = np.array(['%scovar'%covartype in x for x in filesList])
		wcorrList = filesList[wcorrData]
		covarList = filesList[covarData]
		wcorrList.sort() # hzB, hzBlPi, hzR, hzRlPi, lzB, lzBlPi....
		covarList.sort() # all Plus, then all Xross

		dataArr = np.array([np.loadtxt(join(path2data, i),skiprows=1) for i in wcorrList])
		dataArr = np.array([[i[:,1],i[:,4]] for i in dataArr]) # +, x signals
		covarArr = np.array([np.loadtxt(join(path2data, i),skiprows=1) for i in covarList])
		halflen = len(covarArr)/2 # 1st half will be + covars, 2nd will be x
		covarSigma = []
		chiFunc = lambda x: chi2.pdf(x, dof)
		# compute chi2, p-values, significance for a) plus and b) cross signals
		for j,ARR in enumerate([covarArr[:halflen],covarArr[halflen:]]):
			print('plus (0), cross (1): %s'%j) # THIS FOR-LOOP BREAKS FOR ALL-Z RUNS
			for i,cov in enumerate(ARR):
				# print('%s'%self.ext_labels[i])
				cov = np.mat(cov)
				# print('covariance:\n',cov)
				sig = np.mat(dataArr[i][j])
				if (~cov.any(axis=0)).any():
					zero_ind = int(np.where(~cov.any(axis=0))[1])
					cov = np.delete(cov,zero_ind,axis=0)
					cov = np.delete(cov,zero_ind,axis=1)
					sig = np.delete(sig,zero_ind)
					print('SINGULAR MATRIX - low-z lack of large-r_p sampling? - discarding bin %s'%(zero_ind+1))
				# print('signal:\n',sig)
				N_d = len(sig)
				N_s = self.Npatch
				invCov = np.linalg.inv(cov) * (N_s-1)/(N_s-N_d-2) # HARTLAP FACTOR
				chi = (sig*invCov)*sig.T
				fchi = float(chi)
				p_val = scint.quad(chiFunc, fchi, np.inf)[0]
				# print('chi2: ',fchi,', p =',p_val)
				xs = abs(stat.norm.interval((1-p_val), loc=0, scale=1)[0])
				covarSigma.append([fchi,p_val,xs])
		covarSigma = np.array(covarSigma)
		
		chi2Stats = np.column_stack((covarList,covarSigma[:,0],covarSigma[:,1],covarSigma[:,2]))
		ascii.write(chi2Stats, join(path2data,'%schi2'%covartype), delimiter='\t', names=['dataset','chi^2','p-val','x-sigma'])

		return None

	def patch_data(self, patchSize, *bitmask_):
		# find survey area
		nside = 2048
		fullSky = 41252.96 # square degrees
		npix = hp.nside2npix(nside)
		pixar = hp.nside2pixarea(nside,degrees=True)
		ra = self.data[self.headers[0]]
		dec = self.data[self.headers[1]]
		theta = np.deg2rad(90.-dec)
		phi = np.deg2rad(ra)
		pixIDs = hp.ang2pix(nside,theta,phi,nest=False)
		GKmap = np.bincount(pixIDs,minlength=npix)
		GKpix = np.where(GKmap!=0,1,0)

		# find masked pixel IDs
		if not self.SDSS:
			kidsBitmap = hp.read_map('/share/splinter/hj/PhD/KiDS_counts_N2048.fits', dtype=int)
			bitmask_cut = [True]*len(kidsBitmap)

			print("PATCHES: CUTTING MASK!=0")
			bitmask_cut = np.where(kidsBitmap==0,True,False)

			if bitmask_[0] != None:
				bitmask_ = bitmask_[0]
				for i in range(len(bitmask_)):
					# construct bitmask cut
					bitmask_cut &= np.where(bitmask_[i] & kidsBitmap == bitmask_[i], False, True)
			lostpixIDs = [j for j,b in enumerate(bitmask_cut) if b==False]
			lostfromGK = np.where(GKpix[lostpixIDs]!=0,1,0) #1=pixel-lost,0=not-lost
			lostmap = np.bincount(lostpixIDs,minlength=npix)
			# hp.write_map(join(self.new_root,'lostpixmap.fits'),lostmap)
			print('Lost npix, fraction of area: %s, %.3f'%(sum(lostfromGK),sum(lostfromGK)/sum(GKpix)))

			if lostpixIDs != []:
				thetaPhis = hp.pix2ang(nside,lostpixIDs)
				# lost pixel coords;
				lostpixra,lostpixdec = np.rad2deg(thetaPhis[1]),(90.-np.rad2deg(thetaPhis[0]))
			else:
				lostpixra,lostpixdec = (np.array([1e3]),np.array([1e3]))
			del kidsBitmap,lostmap
			gc.collect()
		else:
			lostpixra,lostpixdec = (np.array([1e3]),np.array([1e3]))

		# divide catalog.data into patches...
		raHist = np.histogram(ra,bins=100)
		decHist = np.histogram(dec,bins=100)
		raCounts, raEdges = raHist
		decCounts, decEdges = decHist
		Counts, Edges = [raCounts, decCounts], [raEdges, decEdges]
		Lowers = [[],[]]
		Uppers = [[],[]]
		# find mins/maxs in ra/dec of surveyed areas
		for i, count in enumerate(Counts):
			while len(count)!=0:
				nonz = np.nonzero(count)[0] # index of 1st populated bin
				Lowers[i].append(Edges[i][nonz[0]]) # corresponding min coord
				if 0 in count:
					leadz = list(count).index(0) # index of 1st empty bin
					Uppers[i].append(Edges[i][leadz]) # corresp. max coord
					# shorten counts/edges, & repeat...
					count = count[leadz:]
					Edges[i] = Edges[i][leadz:]
					nonz = np.nonzero(count)[0]
					count = count[nonz[0]:]
					Edges[i] = Edges[i][nonz[0]:]
				else: # ...until absolute maxs in ra/dec
					Uppers[i].append(Edges[i][-1])
					count=[] # kill while-loop

		raLs,raUs,decLs,decUs = map(lambda x: np.array(x),[Lowers[0],Uppers[0],Lowers[1],Uppers[1]])
		print('raLs,raUs,decLs,decUs: ',raLs,raUs,decLs,decUs)
		# ranges in ra/dec
		deltaR = raUs-raLs
		deltaD = decUs-decLs
		print('deltaR, deltaD', deltaR, deltaD)
		RDratio = deltaR/deltaD
		# use ratio to determine integer nos. of 'columns/rows' of patches
		dLen = np.array([1]*len(RDratio))
		rLen = dLen*np.round(RDratio)
		rPatchside = deltaR/rLen
		dPatchside = deltaD/dLen
		patchAr = rPatchside*dPatchside

		while any(abs(patchAr-patchSize)>0.5):
			initDiff = abs(patchAr-patchSize)
			dLen += 1
			rLen = dLen*np.round(RDratio)
			rPatchside = deltaR/rLen
			dPatchside = deltaD/dLen
			patchAr = rPatchside*dPatchside
			print('try patchArs :',patchAr)
			newDiff = abs(patchAr-patchSize)
			if newDiff[0]>initDiff[0]:
				print('MISSED PATCH-AREA TARGET, reverting to closest..')
				dLen -= 1
				rLen = dLen*np.round(RDratio)
				rPatchside = deltaR/rLen
				dPatchside = deltaD/dLen
				patchAr = rPatchside*dPatchside
				print('patchArs :',patchAr)
				break
		[print('region %s: %.2f deg^2'%(j+1,patchAr[j])) for j in range(len(patchAr))] 
		print('ra sides (#,deg): ',rLen,rPatchside)
		print('dec sides (#,deg): ',dLen,dPatchside)
		patch_Ars = np.array([[i]*(dLen[j])*(rLen[j]) for j,i in enumerate(patchAr)]).flatten()

		# contsruct patch edges = 'ra/decPatches'
		raLims = np.column_stack((raLs,raUs))
		decLims = np.column_stack((decLs,decUs))

		raPatches = [np.linspace(raLims[i][0],raLims[i][1],num=rLen[i]+1) for i in range(len(raLims))]
		decPatches = [np.linspace(decLims[i][0],decLims[i][1],num=dLen[i]+1) for i in range(len(decLims))]
		raPatches,decPatches = map(lambda x: np.array(x),[raPatches,decPatches])

		# create column/row cuts from patch edges
		raCuts = []
		decCuts = []
		pixraCuts = []
		pixdecCuts = []
		for j in range(len(raPatches)):
			raCuts.append(np.array(
				[np.where((ra>=raPatches[j][i])&(ra<=raPatches[j][i+1]),True,False) for i in range(len(raPatches[j])-1)]
				)
			)
			# and define from lostpixra
			pixraCuts.append(np.array(
				[np.where((lostpixra>=raPatches[j][i])&(lostpixra<=raPatches[j][i+1]),True,False) for i in range(len(raPatches[j])-1)]
				)
			)
		for j in range(len(decPatches)):
			decCuts.append(np.array(
				[np.where((dec>=decPatches[j][i])&(dec<=decPatches[j][i+1]),True,False) for i in range(len(decPatches[j])-1)]
				)
			)
			# and define from lostpixdec
			pixdecCuts.append(np.array(
				[np.where((lostpixdec>=decPatches[j][i])&(lostpixdec<=decPatches[j][i+1]),True,False) for i in range(len(decPatches[j])-1)]
				)
			)

		raCuts,decCuts,pixraCuts,pixdecCuts = map(lambda x: np.array(x),[raCuts,decCuts,pixraCuts,pixdecCuts])
		raCuts,decCuts,pixraCuts,pixdecCuts = map(lambda x: x.reshape(-1,x.shape[-1]),[raCuts,decCuts,pixraCuts,pixdecCuts]) # flatten 3 cut-arrays (from 3 survey regions) into 1

		# combine into patch-cuts
		patchCuts = []
		pixpatchCuts = []
		for j in decCuts:
			[patchCuts.append(i&j) for i in raCuts]
		for j in pixdecCuts:
			[pixpatchCuts.append(i&j) for i in pixraCuts]
		patchCuts = np.array(patchCuts)
		pixpatchCuts = np.array(pixpatchCuts)
		print('pixpatchCuts shape (patches,pixels): (%s, %s)'%pixpatchCuts.shape)
		assert patchCuts.shape[0] == raCuts.shape[0]*decCuts.shape[0],'patch-cuts broken'

		# combine patch & z/colour cuts, coerce into ndarrays for slicing
		highzR_pcuts = [(self.zcut&self.redcut&self.bitmaskcut&self.BCG_sc&self.lmstarcut&pc) for pc in patchCuts]
		highzB_pcuts = [(self.zcut&self.bluecut&self.bitmaskcut&self.BCG_sc&self.lmstarcut&pc) for pc in patchCuts]
		lowzR_pcuts = [(self.zcut_r&self.redcut&self.bitmaskcut&self.BCG_sc&self.lmstarcut&pc) for pc in patchCuts]
		lowzB_pcuts = [(self.zcut_r&self.bluecut&self.bitmaskcut&self.BCG_sc&self.lmstarcut&pc) for pc in patchCuts]

		highzR_UnM_pcuts = [(self.zcut&self.redcut&self.BCG_sc&self.lmstarcut&pc) for pc in patchCuts]
		highzB_UnM_pcuts = [(self.zcut&self.bluecut&self.BCG_sc&self.lmstarcut&pc) for pc in patchCuts]
		lowzR_UnM_pcuts = [(self.zcut_r&self.redcut&self.BCG_sc&self.lmstarcut&pc) for pc in patchCuts]
		lowzB_UnM_pcuts = [(self.zcut_r&self.bluecut&self.BCG_sc&self.lmstarcut&pc) for pc in patchCuts]

		highz_pcuts = [(self.zcut&self.BCG_dc&self.lmstarcut&pc) for pc in patchCuts]
		lowz_pcuts = [(self.zcut_r&self.BCG_dc&self.lmstarcut&pc) for pc in patchCuts]

		# all lists of patch-cuts -> arrays
		highzR_pcuts,highzB_pcuts,lowzR_pcuts,lowzB_pcuts, highzR_UnM_pcuts,highzB_UnM_pcuts,lowzR_UnM_pcuts,lowzB_UnM_pcuts, highz_pcuts,lowz_pcuts = map(lambda x: np.array(x),[highzR_pcuts,highzB_pcuts,lowzR_pcuts,lowzB_pcuts, highzR_UnM_pcuts,highzB_UnM_pcuts,lowzR_UnM_pcuts,lowzB_UnM_pcuts, highz_pcuts,lowz_pcuts])

		# cut data into 10xpatch-arrays, w/ each element is a fits-table
		hizR_patches,hizB_patches,lozR_patches,lozB_patches, hizR_UnM_patches,hizB_UnM_patches,lozR_UnM_patches,lozB_UnM_patches, hiz_patches,loz_patches = map(lambda x: np.array([self.data[pc] for pc in x]),[highzR_pcuts,highzB_pcuts,lowzR_pcuts,lowzB_pcuts, highzR_UnM_pcuts,highzB_UnM_pcuts,lowzR_UnM_pcuts,lowzB_UnM_pcuts, highz_pcuts,lowz_pcuts])

		del highzR_pcuts,highzB_pcuts,lowzR_pcuts,lowzB_pcuts,highzR_UnM_pcuts,highzB_UnM_pcuts,lowzR_UnM_pcuts,lowzB_UnM_pcuts,highz_pcuts,lowz_pcuts
		gc.collect()

		self.patchedData = np.array([hizR_patches,hizB_patches,lozR_patches,lozB_patches, hizR_UnM_patches,hizB_UnM_patches,lozR_UnM_patches,lozB_UnM_patches, hiz_patches,loz_patches])

		# count lost pixels within each patch for weighting
		npixLost = np.array([np.count_nonzero(i) for i in pixpatchCuts])
		pArLost = npixLost*pixar
		weights = 1-(pArLost/patch_Ars)
		print('patch weights (check none << 0) :\n',weights)
		weights = np.where(weights>0,weights,0)
		self.patchWeights = weights
		# print('patch weights: ',self.patchWeights)

		return self.patchedData, self.patchWeights

	def save_patches(self, patch, outfile_root, label, p_num, largePi):
		if largePi:
			label += '_largePi'
		patchDir = join(outfile_root,label)
		if not isdir(patchDir):
			mkdir(patchDir)
		patchName = join(patchDir,label+'%spatch.asc'%str(p_num).zfill(2))
		ascii.write(patch, patchName, delimiter='\t', names=['#RA/rad', '#DEC/rad', '#comov_dist/Mpc/h', '#e1', '#e2', '#e_weight'])
		return patchDir

	def wcorr_patches(self, patchDir, rp_bins, rp_lims, los_bins, los_lim, nproc, largePi):
		patches = [patch for patch in listdir(patchDir) if ('wcorr' not in patch)&('_' in patch)]
		patches.sort()
		label = basename(normpath(patchDir))
		if 'highZ' in label:
			dDir = join(patchDir,'..','highZ')
		if 'lowZ' in label:
			dDir = join(patchDir,'..','lowZ')
		dpatches = listdir(dDir)
		dpatches.sort()
		os.system('cp %s/* %s'%(dDir,patchDir))
		print('correlating patch (BCG=%s) densities with %s (BCG=%s) shapes...'%(self.BCGargs[0],label,self.BCGargs[1]))
		# print('is BCGs: shape=%s, dens=%s'%(self.BCGargs[0],self.BCGargs[1]))
		for i,p in enumerate(patches):
			pCount = len(np.loadtxt(join(patchDir,p)))
			dCount = len(np.loadtxt(join(patchDir,dpatches[i])))
			# print("patch %s, density popn %s, shapes popn %s"%((i+1),dCount,pCount))
			os.system('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s %s 0 0'%(patchDir,dpatches[i],dCount,p,pCount,rp_bins,rp_lims[0],rp_lims[1],los_bins,los_lim,p[:-9],nproc))
			if largePi:
				os.system('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s_largePi %s 1 0'%(patchDir,dpatches[i],dCount,p,pCount,rp_bins,rp_lims[0],rp_lims[1],los_bins,los_lim,p[:-9],nproc))
			del pCount,dCount
			gc.collect()

	def bootstrap_signals(self, patchDir, patchWeights, largePi):
		pwcorrs = [x for x in listdir(patchDir) if ('.dat' in x)&('Pi' not in x)]
		if largePi:
			pwcorrs = [x for x in listdir(patchDir) if 'Pi' in x]
		psigs = []
		for i,x in enumerate(pwcorrs):
			path = join(patchDir,x)
			data = np.array(np.loadtxt(path))
			psigs.append([data[:,3],data[:,4],[patchWeights[i]]*len(data[:,3])])
			# [+, x, pweight]
		psigs = np.array(psigs)
		BTsignals = astat.bootstrap(psigs,bootnum=100000)
		# construct errors for each bin in r_p, for + & x corrs
		wgplus = BTsignals[:,:,0,:]
		wgcross = BTsignals[:,:,1,:]
		pws = BTsignals[:,:,2,:]
		# shape = (BTs, patch-signal/weight, rp bins)
		del BTsignals,psigs
		gc.collect()

		# calculate mean over patches (weighted)
		Pmeans = np.average(wgplus,axis=1,weights=pws)
		Xmeans = np.average(wgcross,axis=1,weights=pws)
		Pmeds = np.median(wgplus,axis=1)
		Xmeds = np.median(wgcross,axis=1)
		# shape = (BTs, mean/med-signal-in-rp-bin)
		del wgplus,wgcross
		gc.collect()

		# calculate covariance matrix & corr-coeffs (for +)
		Cp,Cx = np.cov(Pmeans,rowvar=0),np.cov(Xmeans,rowvar=0)
		Rp,Rx = np.corrcoef(Pmeans,rowvar=0),np.corrcoef(Xmeans,rowvar=0)

		# calculate stdev over BT-realis'ns
		Pstds = np.sqrt(np.diag(Cp))
		Xstds = np.sqrt(np.diag(Cx))
		# shape = (stdev-on-means-in-rp-bins)

		BTstds = np.column_stack((Pstds,Xstds))
		label = basename(normpath(patchDir))
		BTerrs_out = join(patchDir,'..','BTerrs_%s'%label)
		if largePi:
			BTerrs_out = join(patchDir,'..','BTerrs_%s_largePi'%label)
		ascii.write(BTstds, BTerrs_out, delimiter='\t', names=['#w(g+)err','#w(gx)err'])
		cov_combos = [[Cp,'P'],[Cx,'X']]

		toplotDir = join(patchDir,'../to_plot')
		if not isdir(toplotDir):
			mkdir(toplotDir)
		if not largePi:
			corrName = join(toplotDir,'BTcorrcoeff_%s'%label)
			ascii.write(Rp, corrName, delimiter='\t')

			for covs in cov_combos:
				covName = join(toplotDir,'BTcovar%s_%s'%(covs[1],label))
				ascii.write(covs[0], covName, delimiter='\t')
		else:
			for covs in cov_combos:
				covName = join(toplotDir,'BTcovar%s_%s_largePi'%(covs[1],label))
				ascii.write(covs[0], covName, delimiter='\t')

		# BTanalysis = np.column_stack((Pmeans,Xmeans,Pmeds,Xmeds))
		# BTanalysis_root = join(toplotDir,'BTanalysis_%s'%label)
		# ascii.write(BTanalysis,BTanalysis_root,delimiter='\t',
		# 	names=['#wg+_pmean_rp%s'%i for i in range(Pmeans.shape[-1])]+['#wgx_pmean_rp%s'%i for i in range(Xmeans.shape[-1])]+['#wg+_pmed_rp%s'%i for i in range(Pmeds.shape[-1])]+['#wgx_pmed_rp%s'%i for i in range(Xmeds.shape[-1])])
		# del BTanalysis
		# gc.collect()
		return None

	def jackknife_patches(self, patchDir, patchWeights):
		# resample patches & save JK samples
		patches = [x for x in listdir(patchDir) if ('patch' in x)&('_' in x)]
		patch_cats = np.array([np.loadtxt(join(patchDir,i)) for i in patches])
		patch_cats = np.array([i for i in patch_cats if (i.shape!=(0,))&(len(i.shape)!=1)])
		JKdir = join(patchDir,'JKsamples')
		if not isdir(JKdir):
			mkdir(JKdir)
		jkweights = np.empty([len(patch_cats)])
		for i,pc in enumerate(patch_cats):
			del_one = np.delete(patch_cats,i,axis=0)
			new_cat = np.concatenate(del_one)
			cat_name = join(JKdir,'JKsample%s.asc'%(str(i).zfill(2)))
			ascii.write(new_cat, cat_name, delimiter='\t', names=['#RA/rad', '#DEC/rad', '#comov_dist/Mpc/h', '#e1', '#e2', '#e_weight'])
			del_one_weights = np.delete(patchWeights,i)
			jkweights[i] = np.mean(del_one_weights)
		self.jkweights = jkweights
		# print('CHECK THIS jackknife sample weights:\n',self.jkweights) 
		# JK SAMPLE WEIGHTS DIFFER ON ORDER <=0.01 FOR PSIZE 9, 4...

		return jkweights

	def wcorr_jackknife(self, patchDir, rp_bins, rp_lims, los_bins, los_lim, nproc, largePi, densColours):
		# wcorr JK samples
		JKdir = join(patchDir,'JKsamples')
		JKsamples = [x for x in listdir(JKdir) if ('.asc' in x)&('JKsample' in x)]
		JKsamples.sort()
		shapdens_dict = {'highZ_Red':'highZ_Red_UnMasked.asc',
				'highZ_Blue':'highZ_Blue_UnMasked.asc',
				'lowZ_Red':'lowZ_Red_UnMasked.asc',
				'lowZ_Blue':'lowZ_Blue_UnMasked.asc',
				'highZ_Red_UnMasked':'highZ_Red_UnMasked.asc',
                                'highZ_Blue_UnMasked':'highZ_Blue_UnMasked.asc',
                                'lowZ_Red_UnMasked':'lowZ_Red_UnMasked.asc',
                                'lowZ_Blue_UnMasked':'lowZ_Blue_UnMasked.asc',
				'highZ_Red_largePi':'highZ_Red_UnMasked.asc',
                                'highZ_Blue_largePi':'highZ_Blue_UnMasked.asc',
                                'lowZ_Red_largePi':'lowZ_Red_UnMasked.asc',
                                'lowZ_Blue_largePi':'lowZ_Blue_UnMasked.asc',
				'highZ_Red_UnMasked_largePi':'highZ_Red_UnMasked.asc',
                                'highZ_Blue_UnMasked_largePi':'highZ_Blue_UnMasked.asc',
                                'lowZ_Red_UnMasked_largePi':'lowZ_Red_UnMasked.asc',
                                'lowZ_Blue_UnMasked_largePi':'lowZ_Blue_UnMasked.asc'}
		if densColours:
			#dlabel = basename(normpath(patchDir))+'_UnMasked.asc'
			dlabel = shapdens_dict[patchDir]
			dens_sample = join(patchDir,'..',dlabel)
		else:
			if 'highZ' in patchDir:
				dens_sample = join(patchDir,'..','highZ.asc')
			if 'lowZ' in patchDir:
				dens_sample = join(patchDir,'..','lowZ.asc')
			dlabel = basename(normpath(dens_sample))
		#slabel = basename(normpath(patchDir))
		slabel = shapdens_dict[patchDir][:-13]
		dCount = len(np.loadtxt(dens_sample))
		print("correlating (BCG=%s) density sample with jackknife_i %s (BCG=%s) shapes (largePi=%s)..."%(self.BCGargs[0],slabel,self.BCGargs[1],largePi))
		os.system('cp %s %s'%(dens_sample,JKdir))
		for i,jk in enumerate(JKsamples):
			jkCount = len(np.loadtxt(join(JKdir,jk)))
			# print("JK%s, shapes popn %s"%((i+1),jkCount))
			if largePi:
				os.system('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s_largePi %s 1 0'%(JKdir,dlabel,dCount,jk,jkCount,rp_bins,rp_lims[0],rp_lims[1],los_bins,los_lim,jk[:-4],nproc))
			else:
				os.system('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s %s 0 0'%(JKdir,dlabel,dCount,jk,jkCount,rp_bins,rp_lims[0],rp_lims[1],los_bins,los_lim,jk[:-4],nproc))
		return None

	def pearson_r(self,covar_matrix):
		c = covar_matrix
		d = np.diag(c)
		stddev = np.sqrt(d.real)
		c /= stddev[:, None]
		c /= stddev[None, :]
		# Clip real and imaginary parts to [-1, 1]
		np.clip(c.real, -1, 1, out=c.real)
		if np.iscomplexobj(c):
		    np.clip(c.imag, -1, 1, out=c.imag)
		return c

	def jackknife(self, patchDir, jkweights, error_scaling, largePi):
		# read-in JK signals
		JKdir = join(patchDir,'JKsamples')
		jkwcorrs = [x for x in listdir(JKdir) if ('Pi' not in x)&('.dat' in x)]
		# N patches <-> N signal files
		if largePi:
			jkwcorrs = [x for x in listdir(JKdir) if ('Pi' in x)&('.dat' in x)]
		jkwcorrs.sort()
		jksignals = np.array([np.loadtxt(join(JKdir,i)) for i in jkwcorrs])

		wgp,wgx = jksignals[:,:,3],jksignals[:,:,4]
		Nobs,Nvar = wgp.shape

		# compute jackknife covariance & pearson-r corrcoeffs
		Cp,Cx = np.cov(wgp,rowvar=0, aweights=jkweights),np.cov(wgx,rowvar=0, aweights=jkweights)

		Cp,Cx = Cp*((Nobs-1)**2)/Nobs, Cx*((Nobs-1)**2)/Nobs # jackknife normalisation
		Cp, Cx = Cp*(error_scaling**(-2)) , Cx*(error_scaling**(-2)) # jackknife area-scaling for lost (edge) patches
		Cp_,Cx_ = copy.copy(Cp),copy.copy(Cx)
		Rp,Rx = self.pearson_r(Cp),self.pearson_r(Cx)

		# compute JK stdev on signals
		Pstds,Xstds = np.sqrt(np.diag(Cp_)),np.sqrt(np.diag(Cx_))
		del jksignals
		gc.collect()

		JKstds = np.column_stack((Pstds,Xstds))
		label = basename(normpath(patchDir))
		if largePi&('largePi' not in label):
			label += '_largePi'
		JKerrs_out = join(patchDir,'..','JKerrs_%s'%label)
		ascii.write(JKstds, JKerrs_out, delimiter='\t', names=['#w(g+)err','#w(gx)err'])
		cov_combos = [[Cp_,'P'],[Cx_,'X']]

		toplotDir = join(patchDir,'../to_plot')
		if not isdir(toplotDir):
			mkdir(toplotDir)
		if not largePi:
			corrName = join(toplotDir,'JKcorrcoeff_%s'%label)
			ascii.write(Rp, corrName, delimiter='\t')

		for covs in cov_combos:
			covName = join(toplotDir,'JKcovar%s_%s'%(covs[1],label))
			ascii.write(covs[0], covName, delimiter='\t')
		return None

	def map_test(self, catalogs):
		npix = hp.nside2npix(128)
		hmap = [0]*npix
		for cat in catalogs:
			ra = cat['RA_1_1']
			dec = cat['DEC_1_1']
			theta = np.deg2rad(90.-dec)
			phi = np.deg2rad(ra)
			pix = hp.ang2pix(128,theta,phi)
			hmap += np.bincount(pix, minlength=npix)
		hp.write_map('patchTest_map.fits',hmap)
		# plt.savefig('/share/splinter/hj/PhD/patchTest.pdf')
		return None

class RandomCatalogue(RealCatalogue):

	def __init__(self, path, densColours, sdss):
		"""""
		read-in catalogue

		"""""
		self.path = path
		self.sdss = sdss
		if not sdss:
			hdulist = fits.open(path)
			self.data = hdulist[1].data
			self.columns = hdulist[1].columns.names
		else:
			self.data = np.loadtxt(path).T
		self.headers = [ ['RA','DEC','Z'], [0,1,2] ][sdss]
		self.labels = [['rand_highZ','rand_lowZ'],['rand_highZ_Red','rand_highZ_Blue','rand_lowZ_Red','rand_lowZ_Blue']][densColours]
		self.samples = []
		self.samplecounts = []

	def cut_data(self, d_sample_Zs):
		"""""
		cut catalogue into redshift subsamples

		"""""
		if not self.sdss:
			assert 'RA' in self.columns, "'RA' not in columns, see column headers: "+ str(self.columns)
			assert 'DEC' in self.columns, "'DEC' not in columns, see column headers: "+ str(self.columns)
			assert 'Z' in self.columns, "'Z' not in columns, see column headers: "+ str(self.columns)
			assert 'RAND_NUM' in self.columns, "'RAND_NUM' not in columns, see column headers: "+ str(self.columns)

		# compute shapes' n(z) & target for rands
		nbin = 10
		print('# z-bins for RANDOM downsampling = %s'%nbin)
		sh_zhist = np.histogram(d_sample_Zs,bins=nbin,range=(0.,0.5))
		sh_nz = sh_zhist[0]
		sh_Nz = sh_nz/len(d_sample_Zs)

		# compute rand n(z) & downsample each bin
		rnd_z = self.data[self.headers[2]]
		rnd_zhist = np.histogram(rnd_z,bins=nbin,range=(0.,0.5))
		zbins = rnd_zhist[1]
		rnd_nz = rnd_zhist[0]
		f_redc = 11*(sh_nz/rnd_nz)

		nztune = np.zeros_like(rnd_z,dtype=bool)
		if 'RAND_NUM' not in self.headers:
			rand_num = np.random.random(len(rnd_z))
		else:
			rand_num = self.data['RAND_NUM']
		for i,nzbin in enumerate(rnd_nz):
			bincut = (zbins[i]<rnd_z)&(rnd_z<=zbins[i+1])
			nzcut = (f_redc[i]<rand_num)&(rand_num<=f_redc[i]*2)
			nzbincut = bincut&nzcut
			nztune = np.where(nzbincut,True,nztune)

		if self.sdss:
			new_dat = (self.data.T[nztune]).T
		else:
			new_dat = self.data[nztune]
		self.samples.append(new_dat)
		self.samplecounts.append(len(new_dat))
		newz = new_dat[self.headers[2]]
		new_nz = np.histogram(newz,bins=nbin,range=(0.,0.5))[0]
		new_Nz = new_nz/len(newz)
		print('real/random N(z) (should be ~1): ',sh_Nz/new_Nz)

	def cut_columns(self, subsample, h): 
		"""""
		take subsample data 
		& isolate columns for wcorr

		"""""

		table = subsample
		RA = np.deg2rad(table[self.headers[0]])
		DEC = np.deg2rad(table[self.headers[1]])
		Z = table[self.headers[2]]
		e1 = 2*np.random.random(len(table))
		e2 = e1
		# e2 *= -1 # for RA increasing leftward, c.f. x-axis increasing rightward
		e_weight = np.ones_like(e1)
		comov = Planck13.comoving_distance(Z)
		comov *= h
		new_table = np.column_stack((RA,DEC,comov,e1,e2,e_weight))
		return new_table

if __name__ == "__main__":
	# Command-line args...
	parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
	parser.add_argument(
	'Catalog',
	help='full path of REAL catalogue to be sampled into ascii table(s)')
	parser.add_argument(
	'Path',
	help='full path of destination directory for subsample ascii catalogues, where further directories will be created and appended with "_z_<z_cut>_c_<colour_cut>" if applicable. ***If using -plotNow, give full path to wcorr output directory***')
	parser.add_argument(
	'-Random',
	help="optional; path to RANDOM catalogue to be correspondingly sampled",
	default=None)
	parser.add_argument(
	'-pgm_cut',
	type=np.float32,
	help='shear polarizability cut, defaults to 0.1',
	default=0.1)
	parser.add_argument(
	'-zCut',
	type=np.float32,
	help='lowZ vs. highZ redshift threshold, between 0 - 0.5. Defaults to 0.22, set to 0 for no redshift cut',
	default=None)
	parser.add_argument(
	'-cCut',
	type=np.float32,
	help='red vs. blue colour threshold, rest-frame g-i (typically>1.0) for KiDSxGAMA, or g-r (typ.>0.63) for SDSS Main. Defaults to None',
	default=None)
	parser.add_argument(
	'-lmstarCut',
	nargs=2,
	type=np.float32,
	help='stellar mass window cut [min, max], logmstar i.e. 10^(x), defaults to None',
	default=None)
	parser.add_argument(
	'-LRGs',
	type=int,
	help='employ BOSS-LOWZ-style observer-mag cuts for LRG selection (1), or not (0), defaults to 0',
	default=0)
	parser.add_argument(
	'-bitmaskCut',
	nargs='*',
	type=int,
	help='list of bitmask IDs (powers of 2) to exclude from catalogue, eg. text-file w/ ROWS; 4, 16, 4096...etc...',
	default=None)
	parser.add_argument(
	'-H',
	type=np.float32,
	help='reduced Planck constant, defaults to 0.7',
	default=0.7)
	parser.add_argument(
	'-rpBins',
	type=int,
	help='specify no. of (log-spaced) bins in comoving transverse separation r_p (Mpc/h), for measurement of density-shape correlations. Defaults to 8',
	default=8)
	parser.add_argument(
	'-rpLims',
	nargs=2,
	type=np.float32,
	help='specify upper & lower (2 args, space-separated) limit in comoving transverse separation r_p (Mpc/h), for measurement of density-shape correlations. Defaults to 0.1 - 60 Mpc/h',
	default=[0.1, 60])
	parser.add_argument(
	'-losBins',
	type=int,
	help='specify no. of bins (positive branch) in comoving line-of-sight separation \Pi (Mpc/h), for measurement of density-shape correlations. Defaults to 15; i.e. 30 bins in total',
	default=15)
	parser.add_argument(
	'-losLim',
	type=np.float32,
	help='specify cut-off in comoving line-of-sight separation \Pi (Mpc/h), for measurement of density-shape correlations. Defaults to 60; s/t range is -60 to 60 Mpc/h',
	default=60)
	parser.add_argument(
	'-nproc',
	type=int,
	help='no. processors to be used in correlation measurement, default = 12',
	default=12)
	parser.add_argument(
	'-largePi',
	type=int,
	choices=[0,1],
	help='specify regular (0) or regular + large-Pi systematics tests (1), defaults to 1',
	default=0)
	parser.add_argument(
	'-wcorr',
	type=int,
	choices=[0,1],
	help='initiate wcorr density-shape correlation measurements (1), or not (0), defaults to 0',
	default=0)
	parser.add_argument(
	'-notes',
	help='notes on any changed wcorr parameters, for appendage to directory name',
	default=None)
	parser.add_argument(
	'-plot',
	help='reduce data & save for easy plotting, after correlations completed (1), or not (0), defaults to 1',
	choices=[0,1],
	type=int,
	default=1)
	parser.add_argument(
	'-plotNow',
	help='bypass sampling functions and reduce/save data for plotting (1), or not (0). Sampling and correlations must be ALREADY COMPLETED i.e. "raw" data already existing. Give arg="Path" as the path to the .dat files (Catalog arg must still be path of readable .fits catalog). Defaults to 0',
	choices=[0,1],
	type=int,
	default=0)
	parser.add_argument(
	'-chiSqu',
	help='bypass sampling functions and calc chi^2 stats from reduced correlation data (1), or not (0). Give arg="Path" as the path to the .dat files (Catalog arg must still be path of readable .fits catalog). Defaults to 0',
	choices=[0,1],
	type=int,
	default=0)
	parser.add_argument(
	'-expec',
	help='model values for chi^2 statistics. Defaults to zeros at all points i.e. chi2 testing for signal detection',
	default=0)
	parser.add_argument(
	'-bootstrap',
	type=int,
	choices=[0,1],
	help='perform bootstrap error determination (1) or not (0), defaults to 1',
	default=1)
	parser.add_argument(
	'-jackknife',
	type=int,
	choices=[0,1],
	help='perform jackknife error determination (1) or not (0), defaults to 1',
	default=1)
	parser.add_argument(
	'-patchSize',
	help='****fix me - will break if not 3, 4, or 9****\npreferred mean patch size (deg^2) for sample covariance determinations, defaults to 9sqdeg - USE THIS ARG IF DOING CHI2',
	type=np.float32,
	default=9)
	parser.add_argument(
	'-BCGdens',
	help='select only BCGs for real density samples (1) or not (0), defaults to 0',
	type=int,
	default=0)
	parser.add_argument(
	'-BCGshap',
	help='select only BCGs for shapes samples (1) or not (0), defaults to 0',
	type=int,
	default=0)
	parser.add_argument(
	'-DEIMOS',
	help='DEIMOS shapes (1), or KSB shapes (0), defaults to 1',
	type=int,
	default=1)
	parser.add_argument(
	'-SDSS',
	help='SDSS Main sample (1), or KiDSxGAMA (0), defaults to 0',
	type=int,
	default=0)
	parser.add_argument(
	'-jk3d',
	help='slice jackknife patches in redshift (1), or not (0), defaults to 1',
	type=int,
	default=1)
	parser.add_argument(
	'-cubeZdepth',
	help='specify target redshift increment of jackknife cubes, default=0.06',
	type=np.float32,
	default=0.06)
	parser.add_argument(
	'-rmagCut',
	help='R-band magnitude above which to exclude faint sources, defaults to 0',
	type=np.float32,
	default=0)
	parser.add_argument(
	'-densColours',
	help='use redshift&colour-cut samples as position samples for correlations (1), or just redshift-cut samples (0), defaults to 0',
	type=int,
	default=0)
	parser.add_argument(
	'-LRG1',
	type=int,
	choices=[0,1],
	default=1,
	help='if selecting LRGs, choose (0) to drop colour/magnitude cut 1 of 6, or (1) to keep the cut. Default=1 for all -LRGX')
	parser.add_argument(
	'-LRG2',
	type=int,
	choices=[0,1],
	default=1,
	help='(0) drop cut 2 of 6, or (1) to keep the cut')
	parser.add_argument(
	'-LRG3',
	type=int,
	choices=[0,1],
	default=1,
	help='(0) drop cut 3 of 6, or (1) to keep the cut')
	parser.add_argument(
	'-LRG4',
	type=int,
	choices=[0,1],
	default=1,
	help='(0) drop cut 4 of 6, or (1) to keep the cut')
	parser.add_argument(
	'-LRG5',
	type=int,
	choices=[0,1],
	default=1,
	help='(0) drop cut 5 of 6, or (1) to keep the cut')
	parser.add_argument(
	'-LRG6',
	type=int,
	choices=[0,1],
	default=1,
	help='(0) drop cut 6 of 6, or (1) to keep the cut')
	args = parser.parse_args()

	catalog = RealCatalogue(args.Catalog, args.DEIMOS, args.rmagCut, args.SDSS)

	if args.plotNow:
		# reduce & save data files, returning filename-list
		print('PLOTTING')
		ldir = listdir(args.Path)
		ia_datfiles = [i for i in catalog.wcorrLabels if 'wcorr_'+i+'.dat' in ldir]
		wcorrOuts = catalog.plot_wcorr(args.Path, ia_datfiles, args.bootstrap, args.jackknife, 0)
		largePi_outs = [basename(normpath(out[:-4] + '_largePi.dat')) for out in wcorrOuts]
		# check for largePi .dat files
		isIn = [i in ldir for i in largePi_outs]
		uniq = np.unique(isIn)
		if uniq.all() == True:
			IDs = [outs[6:-4] for outs in largePi_outs]
			a = catalog.plot_wcorr(args.Path, IDs, args.bootstrap, args.jackknife, 1)

		if args.chiSqu:
			print('CALC CHI^2')
			# calculate chi^2 statistics & save to ascii
			if args.bootstrap:
				catalog.chi2(join(args.Path,'to_plot'), args.expec, args.rpBins, 'BT')
				print('Bootstrap chi2 done')
			if args.jackknife:
				catalog.chi2(join(args.Path,'to_plot'), args.expec, args.rpBins, 'JK')
				print('Jackknife chi2 done')
			else:
				print('No sample covariances estimated -> no chi^2')
			sys.exit()
		sys.exit()

	if args.chiSqu:
		print('CALC CHI^2')
		# calculate chi^2 statistics & save to csv
		if args.bootstrap:
			catalog.chi2(join(args.Path,'to_plot'), args.expec, args.rpBins, 'BT')
			print('Bootstrap chi2 done')
		if args.jackknife:
			catalog.chi2(join(args.Path,'to_plot'), args.expec, args.rpBins, 'JK')
			print('Jackknife chi2 done')
		else:
			print('No sample covariances estimated -> no chi^2')
		sys.exit()

	catalog.cut_data(args.pgm_cut, args.zCut, args.cCut, args.lmstarCut, args.LRGs, args.LRG1, args.LRG2, args.LRG3, args.LRG4, args.LRG5, args.LRG6, args.BCGdens, args.BCGshap, args.bitmaskCut)
	samples = [catalog.highz_R,catalog.highz_B,catalog.lowz_R,catalog.lowz_B, 
				catalog.highz_R_UnM,catalog.highz_B_UnM,catalog.lowz_R_UnM,catalog.lowz_B_UnM, 
				catalog.highz,catalog.lowz]
	cuts = 'z-cut: %s\t colour-cut (g-i): %s'%(args.zCut,args.cCut)
	outfile_root = join(args.Path,'Wcorr')

	print('CUTTING/SAVING SAMPLES...')
	sample_zs = []
	swot_z = []
	for i, sample in enumerate(samples):
		new_table,sample_z = catalog.cut_columns(sample, args.H)
		sample_zs.append(sample_z)
		sample_num = catalog.save_tables(new_table, outfile_root, catalog.labels[i], args.zCut, args.cCut, args.notes)
		np.savetxt(join(catalog.new_root,catalog.labels[i]+'_galZs.txt'),sample_z)
		if i>3:
			swot_z.append(catalog.save_swotfiles(sample,catalog.labels[i]))
	print('SAVING TO DIRECTORY: %s'%catalog.new_root)
	if args.densColours:
		sample_zs = sample_zs[4:-2]
		swot_z = swot_z[:-2]
	else:
		sample_zs = sample_zs[8:]
		swot_z = swot_z[-2:]

	if args.bootstrap or args.jackknife:
		print('COMPUTING SAMPLE COVARIANCES...')

		# jkData.shape = (10 subsamples, N patches/cubes)
		# jkData, jkWeights = catalog.patch_data(args.patchSize, args.bitmaskCut)
		jkData, jkWeights, error_scaling = jk3d.resample_data(catalog.data, catalog.samplecuts, patchside=args.patchSize, do_sdss=args.SDSS, do_3d=args.jk3d, cube_zdepth=args.cubeZdepth, largePi=0, bitmaskCut=args.bitmaskCut)
		Njkregions = len(jkData[0])

		for i,sam in enumerate(jkData):
			p_pops = np.array([len(x) for x in sam])
			popd_sam = sam[p_pops!=0]
			for j,p in enumerate(popd_sam):
				new_p,patch_z = catalog.cut_columns(p, args.H)
				pDir = catalog.save_patches(new_p, catalog.new_root, catalog.labels[i], j, 0) # save_patches returns str(patchDir)
				if i>3:
					catalog.save_swotpatches(p, catalog.labels[i], j)
		for lab in catalog.labels[:4]:
			pDir = join(catalog.new_root,lab)
			if ((args.zCut==None)&(lab.startswith('low')))|((args.cCut==None)&('Blue' in lab)):
				print('no z/colour-cut; skipping %s..'%lab)
			else:
				if args.jackknife:
					jksample_weights = catalog.jackknife_patches(pDir, jkWeights)
					catalog.wcorr_jackknife(pDir, args.rpBins, args.rpLims, args.losBins, args.losLim, args.nproc, 0, args.densColours)
					catalog.jackknife(pDir, jksample_weights, error_scaling, 0)
				# if args.largePi:
				# 	if args.jackknife:
				# 		jksample_weights = catalog.jackknife_patches(pDir, jkWeights)
				# 		catalog.wcorr_jackknife(pDir, args.rpBins, args.rpLims, args.losBins, args.losLim, args.nproc, 1, args.densColours)
				# 		catalog.jackknife(pDir, jksample_weights, error_scaling, 1)
		if args.largePi:
			jkData, jkWeights, error_scaling = jk3d.resample_data(catalog.data, catalog.samplecuts, patchside=args.patchSize, do_sdss=args.SDSS, do_3d=args.jk3d, cube_zdepth=args.cubeZdepth, largePi=1, bitmaskCut=args.bitmaskCut)
			Njkregions = len(jkData[0])

		for i,sam in enumerate(jkData):
			for j,p in enumerate(sam):
				new_p,patch_z = catalog.cut_columns(p, args.H)
				pDir = catalog.save_patches(new_p, catalog.new_root, catalog.labels[i], j, 1)

			if args.jackknife:
				jksample_weights = catalog.jackknife_patches(pDir, jkWeights)
				catalog.wcorr_jackknife(pDir, args.rpBins, args.rpLims, args.losBins, args.losLim, args.nproc, 1, args.densColours)
				catalog.jackknife(pDir, jksample_weights, error_scaling, 1)

	catalog.make_combos(args.densColours)

	if (args.zCut==None)&(args.cCut==None):
		adjusted_combos = [catalog.wcorr_combos[0]]
	elif (args.zCut==None)&(args.cCut!=None):
		adjusted_combos = catalog.wcorr_combos[:2]
	elif (args.zCut!=None)&(args.cCut==None):
		adjusted_combos = np.array(catalog.wcorr_combos)[np.array([0,2])]
	elif (args.zCut!=None)&(args.cCut!=None):
		adjusted_combos = catalog.wcorr_combos
	print('CHECK THIS z/colour cut-adjusted combos:\n',adjusted_combos)

	catalog.prep_wcorr(catalog.new_root,adjusted_combos,args.rpBins,args.rpLims,args.losBins,args.losLim,args.nproc,args.largePi,'real_wcorr')

	if args.wcorr:
		print('QSUBBING REAL_WCORR..')
		list_dir = np.array(listdir(catalog.new_root))
		shells = np.array([i.endswith('.sh') for i in list_dir])
		r_shells = np.array([i.startswith('real') for i in list_dir])
		list_dir = list_dir[(shells&r_shells)]
		[os.system('qsub '+ join(catalog.new_root, shell)) for shell in list_dir]

	if args.Random != None:
		catalog2 = RandomCatalogue(args.Random, args.densColours, args.SDSS)
		for sam_z in sample_zs:
			catalog2.cut_data(sam_z) # generates randoms samples & counts

		[print('# objects %s: \t'%catalog2.labels[i],v) for i,v in enumerate(catalog2.samplecounts)]

		cuts = 'z-cut: %s'%args.zCut

		for i, sample in enumerate(catalog2.samples):
			label = catalog2.labels[i]
			if ((args.zCut==None)&('lowZ' in label))|((args.cCut==None)&('Blue' in label)):
				print('no z/colour-cut; skipping %s...'%label)
			else:
				print('CUTTING/SAVING RANDOMS (%s)...'%label)
				new_table = catalog2.cut_columns(sample, args.H)
				sample_num = catalog2.save_tables(new_table,outfile_root,label,args.zCut,args.cCut,args.notes)
				rand_sw_zs = catalog2.save_swotfiles(sample,label)

		rand_ind = [[0,0,1,1], [0,1,2,3]][args.densColours]
		rand_combos = [
		[catalog2.labels[rand_ind[0]]+'.asc', catalog2.samplecounts[rand_ind[0]], catalog.labels[0]+'.asc', catalog.samplecounts[0], 'rand_'+catalog.wcorrLabels[0]],
		[catalog2.labels[rand_ind[1]]+'.asc', catalog2.samplecounts[rand_ind[1]], catalog.labels[1]+'.asc', catalog.samplecounts[1], 'rand_'+catalog.wcorrLabels[1]],
		[catalog2.labels[rand_ind[2]]+'.asc', catalog2.samplecounts[rand_ind[2]], catalog.labels[2]+'.asc', catalog.samplecounts[2], 'rand_'+catalog.wcorrLabels[2]],
		[catalog2.labels[rand_ind[3]]+'.asc', catalog2.samplecounts[rand_ind[3]], catalog.labels[3]+'.asc', catalog.samplecounts[3], 'rand_'+catalog.wcorrLabels[3]]
		]
		if (args.zCut==None)&(args.cCut==None):
			rand_combos = [rand_combos[0]]
		elif (args.zCut!=None)&(args.cCut==None):
			rand_combos = np.array(rand_combos)[np.array([0,2])]
		elif (args.zCut==None)&(args.cCut!=None):
			rand_combos = rand_combos[:2]
		print('CHECK THIS rand combos:\n',rand_combos)

		catalog2.prep_wcorr(catalog.new_root, rand_combos, args.rpBins, args.rpLims, args.losBins, args.losLim, args.nproc, args.largePi, 'rand_wcorr')

		if args.plot:
			with open(join(catalog.new_root, 'rand_wcorr.sh'), 'a') as script:
				script.write(
				'\npython /share/splinter/hj/PhD/catalog_sampler.py %s %s -patchSize %s -plotNow 1 -chiSqu 1 -bootstrap %s -jackknife %s'%(args.Catalog,catalog.new_root,args.patchSize,args.bootstrap,args.jackknife)
				)
				script.write('\n')

		if args.wcorr:
			print('QSUBBING RANDOM_WCORR..')
			list_dir = np.array(listdir(catalog.new_root))
			shells = np.array([i.endswith('.sh') for i in list_dir])
			r_shells = np.array([i.startswith('rand') for i in list_dir])
			list_dir = list_dir[(shells&r_shells)]
			[os.system('qsub '+ join(catalog.new_root, shell)) for shell in list_dir]

	print('FINISHED!')
	with open(join(catalog.new_root,'C-lineArgs_SampleProps.txt'),'a') as script:
		script.write(str(args))
		script.write('\n')
		[script.write('%s: \t%d, mean R_mag: %.4f\n'%(catalog.labels[i],catalog.samplecounts[i],catalog.Rmags[i])) for i in range(len(catalog.labels[:4]))]
		[script.write('%s: \t%d\n'%(catalog.labels[i+4],catalog.samplecounts[i+4])) for i in range(len(catalog.labels[4:]))]
		[script.write('%s: \t%d\n'%(catalog2.labels[i],catalog2.samplecounts[i])) for i in range(len(catalog2.labels))]
	os.system('cp %s %s'%(join(catalog.new_root,'C-lineArgs_SampleProps.txt'),join(catalog.new_root,'to_plot')))
	np.savetxt(join(catalog.new_root,'Rmags.txt'),catalog.Rmags,header='mean absmag_r; hzr,hzb,lzr,lzb')
	os.system('cp %s %s'%(join(catalog.new_root,'Rmags.txt'),join(catalog.new_root,'to_plot')))








