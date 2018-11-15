#!/user/bin/env python
from __future__ import print_function, division
from astropy.io import fits
import numpy as np
from astropy.io import ascii
from os.path import join, isdir, basename, normpath, dirname
from os import listdir, mkdir
import copy
import os
import argparse
import csv
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM as FLCDM
MICEcosmo = FLCDM(Om0=0.25, H0=70, Ob0=0.044)
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
import downsampler as ds
ds_jkfunc = ds.make_jks
import pickle

class RealCatalogue:

	def __init__(self, path, DEI, mc, SDSS, other=0, radians=0, cols=None, largePi=0, MICEdensity=0, SHIFT=0, PLOTNOW=0): # ADD MORE self.SPECS HERE; FEWER ARGS FOR FNS!
		"""""
		read-in catalogue

		"""""
		self.path = path
		self.largePi = largePi
		self.MICEdensity = MICEdensity
		GAMAheads = ['ra', 'dec', 'z', 'e1', 'e2', 'RankBCG', 'logmstar', 'pgm', 'absmag_g', 'absmag_r', 'mask']
		SDSSheads = ['ra', 'dec', 'z', 'e1', 'e2', 'sigma_gamma', 'rf_g-r', 'absmag_r']
		if DEI:					################# GET RID OF PGM #################
			self.headers = dict( zip( GAMAheads, ['RA_GAMA', 'DEC_GAMA', 'Z_TONRY', 'e1_r', 'e2_r', 'RankBCG', 'logmstar', 'pgm', 'absmag_g', 'absmag_r', 'MASK_r'] ) )
			self.DEI = 1
			self.SDSS = 0
		if SDSS:
			self.headers = dict( zip( SDSSheads, SDSSheads ) )
			self.DEI = DEI = 0
			self.SDSS = 1
		if not DEI | SDSS:
			self.headers = dict( zip( GAMAheads, ['RA_1_1', 'DEC_1_1', 'Z_TONRY', 'e1c', 'e2c', 'RankBCG_1', 'logmstar', 'pgm', 'absmag_g_1', 'absmag_r_1', 'col3'] ) )
			self.DEI = 0
			self.SDSS = 0
		self.other = other

		if cols!=None:
			if DEI:
				self.headers = dict( zip( GAMAheads, cols ) )
			if SDSS:
				self.headers = dict( zip( SDSSheads, cols ) )
			print('specifying column headers:')
			for k in (GAMAheads, SDSSheads)[SDSS]:
				print('%s: %s'%(k, self.headers[k]))

		hdulist = fits.open(path)
		data = hdulist[1].data

		if (not SDSS) & SHIFT:
			g12_ra = data[ self.headers['ra'] ]
			shifted_dec = data[ self.headers['dec'] ]
			G12 = (g12_ra > 170) & (g12_ra < 190)
			shifted_dec = np.where(G12, shifted_dec + 1, shifted_dec)
			data[ self.headers['dec'] ] = shifted_dec

		print('SELECTING R_MAG < %s'%mc)
		if DEI:
			fs = data['fluxscale']
		else:
			fs = np.ones(len(data))

		if not PLOTNOW:
			Mr = data[self.headers['absmag_r']]
			Mr = np.where(fs >= 1., Mr - 2.5*np.log10(fs), Mr)
			data = data[Mr < mc]

			if not self.other:
				print('cutting 0.02 < z < 0.5 randoms cannot extend this near, and GAMA too sparse after z=0.5')
				data = data[ (0.02 < data[self.headers['z']]) & (data[self.headers['z']] <= 0.5) ]

			if radians:
				data[self.headers['ra']] = data[self.headers['ra']] * (180./np.pi)
				data[self.headers['dec']] = data[self.headers['dec']] * (180./np.pi)

		self.data = data
		del data
		gc.collect()

		self.columns = hdulist[1].columns.names
		# ascii (.asc) file IDs
		self.labels = ['highZ_Red', 'highZ_Blue', 'lowZ_Red', 'lowZ_Blue',
				'highZ_Red_UnMasked', 'highZ_Blue_UnMasked', 'lowZ_Red_UnMasked', 'lowZ_Blue_UnMasked',
				'highZ', 'lowZ']
		ext_labels = ['%s_largePi'%i for i in self.labels[:4]]
		ext_labels = self.labels[:4] + ext_labels
		ext_labels.sort()
		self.ext_labels = ext_labels
		# bodies of wcorr output file IDs
		self.wcorrLabels = ['highZ_vs_highZ_Red', 'highZ_vs_highZ_Blue', 'lowZ_vs_lowZ_Red', 'lowZ_vs_lowZ_Blue']
		# Npatch = {3.0:96,4.0:54,9.0:24}
		# self.Npatch = Npatch[psize]

		# MEASURE WG+ SIGNALS

	def cut_data(self, pgm_=None, z_=None, colour_=None, lmstar_=None, LRG=0, LRG1=1, LRG2=1, LRG3=1, LRG4=1, LRG5=1, LRG6=1, BCGdens=0, BCGshap=0, bitmask_=None):
		"""""
		cut catalogue according to bitmasks, 
		PGM, & into subsamples

		"""""
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
		z = self.data[self.headers['z']]
		self.pre_z = z
		if self.DEI:
			colour = self.data[self.headers['absmag_g']] - self.data[self.headers['absmag_r']]
			#if ('gminusi' in self.data.columns.names):
				#print('LATEST GAMA; colour = "gminusi"')
				#colour = self.data['gminusi']
		if not self.SDSS:
			total_bitmasks = self.data[self.headers['mask']]
			logmstar = self.data[self.headers['logmstar']]+np.log10(self.data['fluxscale'])

		if self.SDSS:
			colour = self.data[self.headers['rf_g-r']]
			#colour = np.where(np.isnan(colour1), self.data['obs_u-r'], colour1)
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

		if pgm_==None:
			pgm_ = 0
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
			#if self.SDSS:
				#red_cut = np.where(np.isnan(colour1), False, red_cut) # where no rf_g-r colour, DISCARD obs frame galaxies
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
		allowed_bits = lambda x: (x==0) | (x==16) | (x==32) | (x==48) # no masking | secondary | tertiary halos | both
		bitmask_cut = np.where(allowed_bits(total_bitmasks), True, False)

		if bitmask_ != None:
			bitmask_ = bitmask_[0]
			for i in range(len(bitmask_)):
				# construct bitmask cut
				bitmask_cut &= np.where(bitmask_[i] & total_bitmasks == bitmask_[i], False, True)

		assert len(bitmask_cut) == len(total_bitmasks), "bitmask testing broken"
		bitmask_cut = np.array(bitmask_cut)
		print('bitmask cut [unique]: \t', np.unique(bitmask_cut))

		if (not self.SDSS) & ('RankBCG' in self.data.columns.names):
			BCGcut = self.data[self.headers['RankBCG']] == 1
			SATcut = ~BCGcut
			BCG_dc = np.ones_like(BCGcut, dtype=bool)
			BCG_sc = np.ones_like(BCGcut, dtype=bool)
			BCGargs = (BCGshap, BCGdens)
			if BCGdens == 1:
				BCG_dc &= BCGcut
			elif BCGdens == 2:
				BCG_dc &= SATcut
			if BCGshap == 1:
				BCG_sc &= BCGcut
			elif BCGshap == 2:
				BCG_sc &= SATcut
		else:
			BCG_dc = BCG_sc = np.ones_like(bitmask_cut, dtype=bool)
			BCGargs = (0, 0)

		# apply cuts
		self.highz_R = self.data[(z_cut & red_cut & bitmask_cut & BCG_sc & lmstar_cut)]
		self.highz_B = self.data[(z_cut & blue_cut & bitmask_cut & BCG_sc & lmstar_cut)]
		self.lowz_R = self.data[(z_cut_r & red_cut & bitmask_cut & BCG_sc & lmstar_cut)]
		self.lowz_B = self.data[(z_cut_r & blue_cut & bitmask_cut & BCG_sc & lmstar_cut)]

		if self.MICEdensity:
			fs = self.data['fluxscale']
			Mr = self.data[self.headers['absmag_r']]
			Mr = np.where(fs >= 1., Mr - 2.5*np.log10(fs), Mr)
			MICE_limit = Mr < -19.67
			BCG_dc &= MICE_limit
		self.highz_R_UnM = self.data[(z_cut & red_cut & BCG_dc & lmstar_cut)]
		self.highz_B_UnM = self.data[(z_cut & blue_cut & BCG_dc & lmstar_cut)]
		self.lowz_R_UnM = self.data[(z_cut_r & red_cut & BCG_dc & lmstar_cut)]
		self.lowz_B_UnM = self.data[(z_cut_r & blue_cut & BCG_dc & lmstar_cut)]

		self.highz = self.data[(z_cut & BCG_dc & lmstar_cut)]
		self.lowz = self.data[(z_cut_r & BCG_dc & lmstar_cut)]

		self.samples = [self.highz_R, self.highz_B, self.lowz_R, self.lowz_B, 
						self.highz_R_UnM, self.highz_B_UnM, self.lowz_R_UnM, self.lowz_B_UnM, 
						self.highz, self.lowz]
		self.keys = ['z2_r', 'z2_b', 'z1_r', 'z1_b',
					'z2_r_unm', 'z2_b_unm', 'z1_r_unm', 'z1_b_unm',
					'z2', 'z1']
		self.samplecounts = []

		# save cuts & sample properties
		self.Rmags = []
		Sprops = {}
		for key, sample in zip(self.keys, self.samples):
            # save mean diff to pivot magnitude
			if self.DEI:
				fs = sample['fluxscale']
			else:
				fs = np.ones(len(sample))
			Mr = sample[self.headers['absmag_r']]
			Mr = np.where(fs >= 1., Mr - 2.5*np.log10(fs), Mr)
			if key in self.keys[:4]:
				self.Rmags.append(np.mean(Mr + 22.))

			Sprops['L_%s'%key] = 10**(-0.4 * np.mean(Mr + 22.))
			Sprops['Z_%s'%key] = np.mean(sample[self.headers['z']])
			if colour_ is None:
				if self.DEI:
					colour = sample[self.headers['absmag_g']] - sample[self.headers['absmag_r']]
				if self.SDSS:
					colour = sample[self.headers['rf_g-r']]
				Sprops['RF_%s'%key] = (colour > 0.66).sum() / float(len(colour))

		self.Sprops = Sprops

		self.zcut = z_cut
		self.zcut_r = z_cut_r
		self.redcut = red_cut
		self.bluecut = blue_cut
		self.lmstarcut = lmstar_cut
		self.bitmaskcut = bitmask_cut
		self.pgmcut = pgm_cut
		self.BCG_dc = BCG_dc # contains MICE limit if appl.
		self.BCG_sc = BCG_sc
		self.BCGargs = BCGargs

		self.hzr_cut = (self.zcut & self.redcut & self.bitmaskcut & self.BCG_sc & self.lmstarcut)
		self.hzb_cut = (self.zcut & self.bluecut & self.bitmaskcut & self.BCG_sc & self.lmstarcut)
		self.lzr_cut = (self.zcut_r & self.redcut & self.bitmaskcut & self.BCG_sc & self.lmstarcut)
		self.lzb_cut = (self.zcut_r & self.bluecut & self.bitmaskcut & self.BCG_sc & self.lmstarcut)
		self.hzrU_cut = (self.zcut & self.redcut & self.BCG_dc & self.lmstarcut)
		self.hzbU_cut = (self.zcut & self.bluecut & self.BCG_dc & self.lmstarcut)
		self.lzrU_cut = (self.zcut_r & self.redcut & self.BCG_dc & self.lmstarcut)
		self.lzbU_cut = (self.zcut_r & self.bluecut & self.BCG_dc & self.lmstarcut)
		self.hzU_cut = (self.zcut & self.BCG_dc & self.lmstarcut)
		self.lzU_cut = (self.zcut_r & self.BCG_dc & self.lmstarcut)
		self.samplecuts = np.array([self.hzr_cut,self.hzb_cut,self.lzr_cut,self.lzb_cut,
							self.hzrU_cut,self.hzbU_cut,self.lzrU_cut,self.lzbU_cut,
							self.hzU_cut,self.lzU_cut])

		return None

	def cut_columns(self, subsample, h, flipe1, flipe2, Kneighbour, R0cut, shapes=0, mbias=(0., 0.)):
		"""""
		take subsample data
		& isolate columns for wcorr

		"""""
		table = subsample
		RA = np.deg2rad(table[self.headers['ra']])
		DEC = np.deg2rad(table[self.headers['dec']])
		Z = table[self.headers['z']]
		if self.DEI|self.SDSS:
			pgm = np.ones_like(Z)
		else:
			pgm = table['pgm']
		e_weight = np.where(pgm<0.1,0,pgm)
		e1 = table[self.headers['e1']]/pgm
		e2 = table[self.headers['e2']]/pgm

		if flipe1:
			#print('FLIPPING e1 !!!!!!')
			e1 *= -1
		if flipe2:
			#print('FLIPPING e2 !!!!!!')
			e2 *= -1

		if self.DEI:
			#flag_cols = [i for i in table.columns.names if i.startswith('flag_DEIMOS')]
			flag_cols = [i for i in table.columns.names if i == 'flag_DEIMOS_r']
			shape_cut = np.ones_like(e1)
			for fc in flag_cols:
				shape_cut *= np.where(table[fc]=='0000', 1, 0)

			if Kneighbour != 0.:
				neighbour_separation = table['Closest_neighbour']
				pair_radii = table['IsoRadius'] + table['Neighbour_IsoRadius']
				shape_cut *= np.where(neighbour_separation > Kneighbour*pair_radii, 1, 0)

			if R0cut != None:
				assert len(R0cut) == 2, "give upper and lower limit for R0_r"
				R = table['R0_r']
				shape_cut *= np.where( (R>R0cut[0]) & (R<R0cut[1]), 1, 0)

			m1, m2 = mbias
			e1 = e1 * (1 + m1)
			e2 = e2 * (1 + m2)

		if self.SDSS:
			shape_cut = np.where( ((abs(e1)>9.9)|(abs(e2)>9.9)), 0, 1)

		e1,e2 = map(lambda x: np.nan_to_num(x), [e1,e2])

		# do not cut unmasked samples
		if not shapes: shape_cut = np.ones_like(shape_cut)
		shape_cut = np.array(shape_cut, dtype=bool)

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
			comov = MICEcosmo.comoving_distance(Z)
		comov *= h
		new_table = np.column_stack((RA,DEC,comov,e1,e2,e_weight))[shape_cut]

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

		ascii.write(new_table, join(outfile_root, label + ".asc"), names=['# ra[rad]', 'dec[rad]', 'chi[Mpc/h]', 'e1', 'e2', 'e_weight'], delimiter='\t', overwrite=1)
		sample_no = "%s # objects:\t%s"%(label,len(new_table))
		return sample_no

	def save_swotfiles(self, subsample, label):
		"""""
		save samples in SWOT format (ra[deg],dec[deg],z)

		"""""
		try:
			RA = subsample[self.headers['ra']]
			DEC = subsample[self.headers['dec']]
			Z = subsample[self.headers['z']]
			e1 = subsample[self.headers['e1']]
			e2 = subsample[self.headers['e2']]
		except IndexError:
			RA = subsample.T[0]
			DEC = subsample.T[1]
			Z = subsample.T[2]
			e1 = np.ones_like(Z)
			e2 = np.ones_like(Z)
		newtable = np.column_stack((RA,DEC,Z,e1,e2))

		ascii.write(newtable, join(self.new_root, 'swot_%s.asc'%label), names=['# ra[deg]', 'dec[deg]', 'z', 'e1', 'e2'], delimiter='\t', overwrite=1)

		return Z # for random downsampling

	def save_swotpatches(self, patch, label, pnum):
		RA = patch[self.headers['ra']]
		DEC = patch[self.headers['dec']]
		Z = patch[self.headers['z']]
		newpatch = np.column_stack((RA,DEC,Z))

		sw_pDir = join(self.new_root,'swot_%s'%label)
		if not isdir(sw_pDir):
			mkdir(sw_pDir)

		ascii.write(newpatch,join(sw_pDir,label+'%spatch.asc'%str(pnum).zfill(3)),names=['# ra[deg]','dec[deg]','z'],delimiter='\t', overwrite=1)

	def make_combos(self, densColours):
		# construct sets of filenames, counts, & IDs for wcorr-calls

		wcorr_ind = [[8,8,9,9], [4,5,6,7]][densColours]

		self.wcorr_combos = [
		[self.labels[wcorr_ind[0]]+'.asc', self.samplecounts[wcorr_ind[0]], self.labels[0]+'.asc', self.samplecounts[0], self.wcorrLabels[0]],
		[self.labels[wcorr_ind[1]]+'.asc', self.samplecounts[wcorr_ind[1]], self.labels[1]+'.asc', self.samplecounts[1], self.wcorrLabels[1]],
		[self.labels[wcorr_ind[2]]+'.asc', self.samplecounts[wcorr_ind[2]], self.labels[2]+'.asc', self.samplecounts[2], self.wcorrLabels[2]],
		[self.labels[wcorr_ind[3]]+'.asc', self.samplecounts[wcorr_ind[3]], self.labels[3]+'.asc', self.samplecounts[3], self.wcorrLabels[3]]
		]
		print('CHECK THIS wcorr call combinations:')
		for wcc in self.wcorr_combos:
			print(wcc)
		self.samplecounts = self.samplecounts[:10]
		[print('# objects %s: \t'%self.labels[i], v) for i, v in enumerate(self.samplecounts)]

	def prep_wcorr(self, files_path, wcorr_combos, rp_bins, rp_lims, los_bins, los_lim, large_pi, nproc, play, out_sh):
		processor_string = 'nodes=1:ppn=16'
		nprocessors = 16
		if play:
			print('OVERRIDING wcorr processors!')
			processor_string = 'nodes=1:ppn=%s'%nproc
			nprocessors = nproc
		shell_script = [
		'#!/bin/tcsh',
		'#PBS -q compute',
		'#PBS -N %s'%out_sh,
		'#PBS -l %s'%processor_string,
		'#PBS -l walltime=24:00:00',
		'#PBS -l mem=50gb',
		'#PBS -o %s'%join(files_path,out_sh+'.out'),
		'#PBS -e %s'%join(files_path,out_sh+'.err'),
		'',
		'date',
		'',
		'source ~/.login',
		'']

		for combo in wcorr_combos: 
			# write 4x wcorr-calls to .sh script, & another 4x if largePi testing
			shell_script.append('')
			outfile = combo[4]
			shell_script.append( ( '/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s %s 0 0' %
									(files_path, combo[0], combo[1], combo[2], combo[3], rp_bins, rp_lims[0], rp_lims[1], los_bins, los_lim, outfile, nprocessors) )
			)
			shell_script.append('')
			if large_pi:
				outfile += '_largePi'
				shell_script.append('')
				shell_script.append( ( '/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s %s 1 0' %
									(files_path, combo[0], combo[1], combo[2], combo[3], rp_bins, rp_lims[0], rp_lims[1], los_bins, los_lim, outfile, nprocessors) )
				)
			shell_script.append('')

			shell_script.append('date')

		File = join(files_path, '%s.sh'%out_sh)
		Write = open(File, 'w')
		Text = '\n'.join(shell_script)
		Write.write(str(Text))
		Write.close()
		os.system('chmod +x %s'%File)

	def plot_treecorr(self, path, largePi=0):
		# needs to create output file in same format as below
		outd = join(path, 'to_plot')
		if not isdir(outd):
			mkdir(outd)

		for label in self.labels:
			wlabel = label.split('_')[0] + '_vs_' + label
			wcorrf = join(path, 'wcorr_' + wlabel + '.dat')
			jkf = join(path, 'JKerrs_' + label)
			if largePi:
				jkf += '_largePi'
				wlabel += '_largePi'
			outf = join(outd, wlabel)

			try:
				rp, wp, wx, werr = np.loadtxt(wcorrf, usecols=(0, 3, 4, 6)).T
				ones = np.ones_like(rp, dtype=float)
				try:
					jkp, jkx = np.loadtxt(jkf).T
				except:
					jkp = jkx = ones
				data = np.column_stack((rp, wp, ones, jkp, wx, ones, jkx, werr))
				ascii.write(data, outf, delimiter='\t', overwrite=1,
							names=['r_p', 'wg+', 'BT+err', 'JK+err', 'wgx', 'BTxerr', 'JKxerr', 'analyticerrs'])
			except IOError:
				continue

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
					if largePi:
						JKerrs = np.loadtxt(join(files_path,'JKerrs_%s_largePi'%new_labels[i]))
					else:
						JKerrs = np.loadtxt(join(files_path,'JKerrs_%s'%new_labels[i]))
				except IOError:
					print('IOError THROWN for JKerrs_%s.. (largePi=%s)!!'%(new_labels[i],largePi))
					JKerrs = np.ones_like(np.array([realErr,randErr]).T)
				Perr2 = JKerrs[:,0]
				Xerr2 = JKerrs[:,1]
				#Pproperr2 = np.sqrt((Perr2**2) + (randErr**2))
				#Xproperr2 = np.sqrt((Xerr2**2) + (randErr**2)) # propagate errors
				Pproperr2 = Perr2.copy()
				Xproperr2 = Xerr2.copy()
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
			ascii.write(reducedData, join(easyPlotDir, basename(normpath(path))[6:-4]), delimiter='\t', names=['r_p', 'wg+', 'BT+err', 'JK+err', 'wgx', 'BTxerr', 'JKxerr', 'analyticerrs'], overwrite=1)

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
		ascii.write(chi2Stats, join(path2data,'%schi2'%covartype), delimiter='\t', names=['dataset','chi^2','p-val','x-sigma'], overwrite=1)

		return None

	def save_patches(self, patch, outfile_root, label, p_num, largePi):
		if largePi:
			label += '_largePi'
		patchDir = join(outfile_root,label)
		if not isdir(patchDir):
			mkdir(patchDir)
		patchName = join(patchDir,label+'%spatch.asc'%str(p_num).zfill(3))
		ascii.write(patch, patchName, delimiter='\t', names=['# ra[rad]', 'dec[rad]', 'chi[Mpc/h]', 'e1', 'e2', 'e_weight'], overwrite=1)
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
		ascii.write(BTstds, BTerrs_out, delimiter='\t', names=['# w(g+)err', 'w(gx)err'], overwrite=1)
		cov_combos = [[Cp,'P'],[Cx,'X']]

		toplotDir = join(patchDir,'../to_plot')
		if not isdir(toplotDir):
			mkdir(toplotDir)
		if not largePi:
			corrName = join(toplotDir,'BTcorrcoeff_%s'%label)
			ascii.write(Rp, corrName, delimiter='\t', overwrite=1)

			for covs in cov_combos:
				covName = join(toplotDir,'BTcovar%s_%s'%(covs[1],label))
				ascii.write(covs[0], covName, delimiter='\t', overwrite=1)
		else:
			for covs in cov_combos:
				covName = join(toplotDir,'BTcovar%s_%s_largePi'%(covs[1],label))
				ascii.write(covs[0], covName, delimiter='\t', overwrite=1)

		return None

	def jackknife_patches(self, patchDir):
		# resample patches & save JK samples
		patches = [x for x in listdir(patchDir) if ('patch' in x)&('_' in x)]
		patch_cats = np.array([np.loadtxt(join(patchDir,i)) for i in patches])

		JKdir = join(patchDir,'JKsamples')
		if not isdir(JKdir):
			mkdir(JKdir)
		for i,pc in enumerate(patch_cats):
			del_one = np.delete(patch_cats,i,axis=0)
			new_cat = np.concatenate(del_one)
			cat_name = join(JKdir,'JKsample%s.asc'%(str(i).zfill(3)))
			ascii.write(new_cat, cat_name, delimiter='\t', names=['# ra[rad]', 'dec[rad]', 'chi[Mpc/h]', 'e1', 'e2', 'e_weight'], overwrite=1)

		return None

	def run_treecorr(self, densf, shapesf, drandf, srandf, config, outfile, **kwargs):
		# config & kwargs constructed from command line args for the main script
		import treecorr_3DCF
		tcw = treecorr_3DCF.compute_w
		config1 = config.copy()

		r, wgp, wgx, err = tcw([densf, shapesf], [drandf, srandf], config1, **kwargs)
		# mimic BJ code output
		one_col = np.ones_like(r, dtype=float)
		out_arr = np.column_stack((r, one_col, one_col, wgp, wgx, one_col, err, one_col)) # need to include shot noise as second-to-last column
		np.savetxt(outfile, out_arr)

	def wcorr_jackknife(self, patchDir, rp_bins, rp_lims, los_bins, los_lim, nproc, largePi, densColours, treecorr=None, **kwargs):
		# wcorr JK samples - this function gets called only for shapes samples
		JKdir = join(patchDir,'JKsamples')
		JKsamples = [x for x in listdir(JKdir) if x.startswith('JKsample') & x.endswith('.asc')]
		JKsamples.sort()
		pdir_key = basename(normpath(patchDir))
		#print('CHECK THESE:\n\n\n\n')
		if densColours:
			dens_dir = join(dirname(patchDir), pdir_key + '_UnMasked')
			if largePi:
				dens_dir = join(dirname(patchDir), pdir_key.replace('_largePi', '_UnMasked_largePi'))
		else:
			if 'highZ' in patchDir:
				dens_dir = join(dirname(patchDir), 'highZ')
			if 'lowZ' in patchDir:
				dens_dir = join(dirname(patchDir), 'lowZ')
			if largePi:
				dens_dir += '_largePi'
		#print('pdir_key:\t', pdir_key)
		#print('dens_dir:\t', dens_dir)
		#print('\n\n\n\n')

		print("correlating (BCG=%s) density sample with jackknife_i %s (BCG=%s) shapes (largePi=%s)..."%(self.BCGargs[0],pdir_key,self.BCGargs[1],largePi))
		for i,jk in enumerate(JKsamples):
			dpath = join(dens_dir, 'JKsamples', jk)
			randjk = 'rand_'+jk
			rdpath = join(dens_dir, 'JKsamples', randjk)
			os.system('cp %s %s_d'%(dpath, join(JKdir, jk)))
			os.system('cp %s %s'%(rdpath, join(JKdir, randjk)))
			dCount = np.loadtxt(dpath).shape[0]
			rdCount = np.loadtxt(rdpath).shape[0]
			jkCount = np.loadtxt(join(JKdir, jk)).shape[0]

			if not treecorr:
				if largePi:
					os.system('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s_d %s %s %s %s %s %s %s %s %s_largePi %s 1 0'%(JKdir,jk,dCount,jk,jkCount,rp_bins,rp_lims[0],rp_lims[1],los_bins,los_lim,jk[:-4],nproc))
					os.system('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s_largePi %s 1 0'%(JKdir,randjk,rdCount,jk,jkCount,rp_bins,rp_lims[0],rp_lims[1],los_bins,los_lim,'rand_'+jk[:-4],nproc))
					real_out, rand_out = join(JKdir, 'wcorr_'+jk[:-4]+'_largePi.dat'), join(JKdir, 'wcorr_rand_'+jk[:-4]+'_largePi.dat')

				else:
					os.system('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s_d %s %s %s %s %s %s %s %s %s %s 0 0'%(JKdir,jk,dCount,jk,jkCount,rp_bins,rp_lims[0],rp_lims[1],los_bins,los_lim,jk[:-4],nproc))
					os.system('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s %s 0 0'%(JKdir,randjk,rdCount,jk,jkCount,rp_bins,rp_lims[0],rp_lims[1],los_bins,los_lim,'rand_'+jk[:-4],nproc))
					real_out, rand_out = join(JKdir, 'wcorr_'+jk[:-4]+'.dat'), join(JKdir, 'wcorr_rand_'+jk[:-4]+'.dat')

				realcorr, randcorr = np.loadtxt(real_out), np.loadtxt(rand_out)
				realcorr[:, 3:5] -= randcorr[:, 3:5]
				np.savetxt(real_out, realcorr)

			elif treecorr:
				# densf, shapesf, drandf, srandf, config, outfile,
				# estim='PW1', np=16, **kwargs(nbins_rpar=30, random_oversampling=10., verbosity=1, load_RRs, save_RRs)
				outf = (join(JKdir,'wcorr_'+jk[:-4]+'.dat'),
						join(JKdir,'wcorr_'+jk[:-4]+'_largePi.dat')) [largePi]
				print('%i / %i' % (i+1, len(JKsamples)))
				self.run_treecorr(dpath, join(JKdir,jk), rdpath, rdpath, kwargs['tc_config'], outf, **kwargs['tc3dcf_kwargs'])

			# clean up - if needing to analyse JKs, use arg 'makejk_only'- bypasses this function
			if not treecorr:
				os.system('rm %s'%rand_out) # random wcorr
			os.system('rm %s_d'%join(JKdir, jk)) # copied density JK sample
			os.system('rm %s'%join(JKdir, randjk)) # copied random JK sample
			os.system('rm %s'%join(JKdir, jk)) # shapes JK sample
			if densColours | ('_Blue' in patchDir):
				os.system('rm %s'%dpath) # density JK sample
				os.system('rm %s'%rdpath) # random sample

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
		print('=======================\t', '/: '.join(patchDir.split('/')[-2:]), ' N surviving JK regions: ', len(jksignals), '....!!!!') 

		wgp, wgx = jksignals[:,:,3], jksignals[:,:,4]
		Nobs, Nvar = wgp.shape

		# compute jackknife covariance & pearson-r corrcoeffs
		Cp,Cx = np.cov(wgp,rowvar=0, aweights=jkweights),np.cov(wgx,rowvar=0, aweights=jkweights)

		# jackknife normalisation: numpy.cov calculates a bootstrap-style normalisation
		# according to [N - 1]^-1 (equiv.to.)= [sum(w)-(sum(w**2)/sum(w))]^-1
		# convert to jackknife-style by *= (N - 1)**2 / N, therefore;
		Norm_jk_nom = jkweights.sum() - ( (jkweights**2).sum() / jkweights.sum() )
		Norm_jk_denom = Norm_jk_nom + 1 
		Norm_jack = Norm_jk_nom**2 / Norm_jk_denom
		# ===>
		Cp, Cx = Norm_jack * Cp, Norm_jack * Cx

		Cp, Cx = Cp*(error_scaling**(-2)) , Cx*(error_scaling**(-2)) # jackknife area-scaling for lost (edge) patches
		Cp_,Cx_ = copy.copy(Cp),copy.copy(Cx)
		Rp,Rx = self.pearson_r(Cp),self.pearson_r(Cx)

		# compute JK stdev on signals
		Pstds,Xstds = np.sqrt(np.diag(Cp_)),np.sqrt(np.diag(Cx_))
		del jksignals
		gc.collect()

		JKstds = np.column_stack((Pstds,Xstds))
		label = basename(normpath(patchDir))
		if largePi & ('largePi' not in label):
			label += '_largePi'
		JKerrs_out = join(patchDir,'..','JKerrs_%s'%label)
		ascii.write(JKstds, JKerrs_out, delimiter='\t', names=['# w(g+)err','w(gx)err'], overwrite=1)
		cov_combos = [[Cp_,'P'],[Cx_,'X']]

		toplotDir = join(patchDir,'../to_plot')
		if not isdir(toplotDir):
			mkdir(toplotDir)
		if not largePi:
			corrName = join(toplotDir,'JKcorrcoeff_%s'%label)
			#ascii.write(Rp, corrName, delimiter='\t', overwrite=1)

		for covs in cov_combos:
			covName = join(toplotDir,'JKcovar%s_%s'%(covs[1],label))
			ascii.write(covs[0], covName, delimiter='\t', overwrite=1)
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

	def __init__(self, path, densColours, sdss, SHIFT=0):
		"""""
		read-in catalogue

		"""""
		self.path = path
		self.sdss = sdss
		hdulist = fits.open(path)
		self.data = hdulist[1].data
		self.columns = hdulist[1].columns.names
		self.headers = dict( zip( ['ra', 'dec', 'z'], [ ['RA','DEC','Z'], ['ra','dec','z'] ][sdss] ) )
		self.labels = [['rand_highZ','rand_lowZ'],['rand_highZ_Red','rand_highZ_Blue','rand_lowZ_Red','rand_lowZ_Blue']][densColours]
		self.samples = []
		self.samplecounts = []

		if (not sdss) & SHIFT:
			random_ra = self.data[ self.headers['ra'] ].copy()
			shifted_dec = self.data[ self.headers['dec'] ].copy()
			G12 = (random_ra > 170) & (random_ra < 190)
			shifted_dec = np.where(G12, shifted_dec + 1, shifted_dec)
			self.data[ self.headers['dec'] ] = shifted_dec

	def cut_columns(self, table, h): 
		"""""
		take subsample data 
		& isolate columns for wcorr

		"""""

		RA = np.deg2rad(table.T[0])
		DEC = np.deg2rad(table.T[1])
		Z = table.T[2]
		e1 = e2 = e_weight = np.ones_like(Z)
		comov = MICEcosmo.comoving_distance(Z)
		comov *= h
		new_table = np.column_stack((RA,DEC,comov,e1,e2,e_weight))
		return new_table

class MyArgumentParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        print(arg_line)
        if arg_line.startswith('#'):
            print('##',arg_line)
            return ''
        else:
            return arg_line.split()

if __name__ == "__main__":
	# Command-line args...
	parser = MyArgumentParser(fromfile_prefix_chars='@')
	parser.add_argument(
	'-Catalog',
	default='MUST_GIVE_CATALOG_PATH',
	help='full path of REAL catalogue to be sampled into ascii table(s)')
	parser.add_argument(
	'-Path',
	default='/share/splinter/hj/PhD/',
	help='full path of destination directory for subsample ascii catalogues, where further directories will be created and appended with "_z_<z_cut>_c_<colour_cut>" if applicable. Default= PhD directory on splinter share. ***If using -plotNow, give full path to wcorr output directory***')
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
	help='lowZ vs. highZ redshift threshold, between 0 - 0.6. Defaults to None',
	default=None)
	parser.add_argument(
	'-cCut',
	type=np.float32,
	help='red vs. blue colour threshold, rest-frame g-r (>0.66) for KiDSxGAMA, or rf_g-r (>0.66) for SDSS Main. Defaults to None',
	default=None)
	parser.add_argument(
	'-Kneighbour',
	type=np.float32,
	help='float between 0 - 1, filters neighbouring pairs whose separation < K*sum(radii) for shapes correlations (DEIMOS). Default=0',
	default=0)
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
	help='reduced Planck constant, defaults to 0.70',
	default=0.70)
	parser.add_argument(
	'-rpBins',
	type=int,
	help='specify no. of (log-spaced) bins in comoving transverse separation r_p (Mpc/h), for measurement of density-shape correlations. Defaults to 11',
	default=11)
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
	help='specify regular (0) or regular + large-Pi systematics tests (1), defaults to 0',
	default=0)
	parser.add_argument(
	'-largePiOnly',
	type=int,
	choices=[0,1],
	help='skip jackknife of normal samples, and go straight to (if) largePi (1), or perform regular JK first (0). Default=0',
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
	help='perform bootstrap error determination (1) or not (0), defaults to 0',
	default=0)
	parser.add_argument(
	'-jackknife',
	type=int,
	choices=[0,1],
	help='perform jackknife error determination (1) or not (0), defaults to 1',
	default=1)
	parser.add_argument(
	'-patchSize',
	help='target angular scale [deg] for jackknife patches, can give 1x float, or 2x (ra, dec), default=4',
	nargs='*',
	type=np.float32,
	default=4)
	parser.add_argument(
	'-BCGdens',
	help='1 = take BCGs only for density // 0 = take all galaxies for density // 2 = take satellites only for density',
	type=int,
	default=0)
	parser.add_argument(
	'-BCGshap',
	help='1 = take BCGs only for shapes // 0 = take all galaxies for shapes // 2 = take satellites only for shapes',
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
	'-flipe1',
	help='flip the sign on e1',
	type=int,
	default=0)
	parser.add_argument(
	'-flipe2',
	help='flip the sign on e2',
	type=int,
	default=0)
	parser.add_argument(
	'-jk3d',
	help='slice jackknife patches in redshift (1), or not (0), defaults to 1',
	type=int,
	default=1)
	parser.add_argument(
	'-occ_thresh',
	help='specify patch-occupation threshold for SDSS survey-edge filtering; float between 0 - 1, default=0.99',
	type=np.float32,
	default=0.67)
	parser.add_argument(
	'-cubeZdepth',
	help='specify target comoving Mpc/h increment of jackknife cubes, default=150',
	type=np.float32,
	default=150.)
	parser.add_argument(
	'-rmagCut',
	help='R-band magnitude above which to exclude faint sources, defaults to 0',
	type=np.float32,
	default=0)
	parser.add_argument(
	'-densColours',
	help='use redshift&colour-cut samples as position samples for correlations (1), or just redshift-cut samples (0), defaults to 1',
	type=int,
	default=1)
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
	parser.add_argument(
	'-makejk_only',
	type=int,
	default=0,
	help='use existing jk correlations to construct covariances (1), or take new correlations (0), default=0')
	parser.add_argument(
	'-cols',
	nargs='*',
	type=str,
	default=None,
	help='specify column headers for GAMA: (ra dec z e1 e2 RankBCG logmstar pgm absmag_g absmag_r mask), or SDSS: (ra dec z e1 e2 sigma_gamma rf_g-r), in order. Default=None; default column headers')
	parser.add_argument(
	'-R0cut',
	nargs=2,
	type=np.float32,
	help='define resolution (r-band) window for non-zero correlation weighting, default=None=all')
	parser.add_argument(
	'-mbias',
	nargs=2,
	type=np.float32,
	default=(-0.0040, -0.0035),
	help='m1, m2 multiplicative ellipticity biases, default=(-0.0040, -0.0035) - C.Georgiou etal. Table 1')
	parser.add_argument(
	'-MICEdens',
	type=int,
	default=0,
	help='1=limit density samples to M_r < -18.9, aligning with MICE sims resolution from which clustering errors and hence sample biases are derived. Default=0')
	parser.add_argument(
	'-make_shear',
	type=int,
	default=0,
	help='1=convert ellipticity correlations/covariances to shear (post-processing), calculating responsivities from shape sample ellipticity distributions. Default=0')
	parser.add_argument(
	'-SHIFT',
	type=int,
	default=0,
	help='1=SHIFT ALL DECLINATIONS IN G12 +1deg - for jackknife testing. Default=0')
	parser.add_argument(
	'-play',
	type=int,
	default=1,
	help='1 = carry out sample wg+ correlations in shell (default), 0 = spawn separate jobs')
	parser.add_argument(
	'-other',
	type=int,
	default=0,
	help='1 = neither GAMA/SDSS; make no redshift boundary cuts before sampling')
	parser.add_argument(
	'-radians',
	type=int,
	default=0,
	help='1 = RA/DEC columns are given in radians; convert to degrees on initialisation')
	parser.add_argument(
	'-unit_weights',
	type=int,
	default=0,
	help='1 = force jackknife patch weights to unity')
	parser.add_argument(
	'-treecorr',
	type=str,
	help='use TreeCorr to measure wg+ -- give path to a configuration file with sections for (i) TreeCorr config, and (ii) treecorr_3DCF script kwargs')
	args = parser.parse_args()
	SHIFT = args.SHIFT

	if args.Catalog.startswith('MUST'):
		print(args.Catalog.split('_'))
		sys.exit()

	print('=======================\tREADING CATALOG: %s'%args.Catalog)

	catalog = RealCatalogue(args.Catalog, args.DEIMOS, args.rmagCut, args.SDSS, other=args.other, radians=args.radians, cols=args.cols, largePi=args.largePi, MICEdensity=args.MICEdens, SHIFT=SHIFT, PLOTNOW=args.plotNow)

	if args.plotNow & (args.treecorr is not None):
		print('PLOTTING TREECORR')
		catalog.plot_treecorr(args.Path, largePi=args.largePi)
		sys.exit()

	elif args.plotNow:
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

	catalog.cut_data(pgm_=args.pgm_cut, z_=args.zCut, colour_=args.cCut, lmstar_=args.lmstarCut, LRG=args.LRGs, LRG1=args.LRG1, LRG2=args.LRG2, LRG3=args.LRG3, LRG4=args.LRG4, LRG5=args.LRG5, LRG6=args.LRG6, BCGdens=args.BCGdens, BCGshap=args.BCGshap, bitmask_=args.bitmaskCut)

	samples = [catalog.highz_R,catalog.highz_B,catalog.lowz_R,catalog.lowz_B, 
				catalog.highz_R_UnM,catalog.highz_B_UnM,catalog.lowz_R_UnM,catalog.lowz_B_UnM, 
				catalog.highz,catalog.lowz]
	
	cuts = 'z-cut: %s\t colour-cut (g-r): %s'%(args.zCut,args.cCut)
	outfile_root = join(args.Path,'Wcorr')

	print('CUTTING/SAVING SAMPLES...')
	sample_zs = []
	swot_z = []
	for i, sample in enumerate(samples):
		if args.flipe1 | args.flipe2:
			print('FLIPPING e1 (%s), e2 (%s)' % (args.flipe1, args.flipe2))
		if i<4:
			shapes=1
		else:
			shapes=0

		new_table,sample_z = catalog.cut_columns(sample, args.H, args.flipe1, args.flipe2, args.Kneighbour, args.R0cut, shapes=shapes, mbias=args.mbias)
		sample_zs.append(sample_z)
		if (args.zCut==None) & ('lowZ' in catalog.labels[i]):
			continue
		sample_num = catalog.save_tables(new_table, outfile_root, catalog.labels[i], args.zCut, args.cCut, args.notes)
		np.savetxt(join(catalog.new_root,catalog.labels[i]+'_galZs.txt'),sample_z)
		if i>3:
			swot_z.append(catalog.save_swotfiles(sample,catalog.labels[i]))
	print('=======================\tSAVING TO DIRECTORY: %s'%catalog.new_root)
	if args.densColours:
		samz_keys = ['z2_r', 'z2_b', 'z1_r', 'z1_b']
		sample_zs = sample_zs[4:-2]
		samz_dict = dict(zip(samz_keys, sample_zs))
		swot_z = swot_z[:-2]
	else:
		samz_keys = ['z2', 'z1']
		sample_zs = sample_zs[8:]
		samz_dict = dict(zip(samz_keys, sample_zs))
		swot_z = swot_z[-2:]

	if args.treecorr:
		import configparser
		cp = configparser.ConfigParser()
		cp.read(args.treecorr)
		top_tc_config = cp._sections
		top_tc_config['tc_config']['min_sep'] = args.rpLims[0]
		top_tc_config['tc_config']['max_sep'] = args.rpLims[1]
		top_tc_config['tc_config']['nbins'] = args.rpBins
		top_tc_config['tc_config']['min_rpar'] = -args.losLim
		top_tc_config['tc_config']['max_rpar'] = args.losLim
		top_tc_config['tc_config']['num_threads'] = args.nproc
		top_tc_config['tc3dcf_kwargs']['nbins_rpar'] = args.losBins * 2
		top_tc_config['tc3dcf_kwargs']['largePi'] = args.largePi
	else:
		top_tc_config = {}

	if args.bootstrap or args.jackknife:
		print('COMPUTING SAMPLE COVARIANCES...')
		if not args.largePiOnly:

			# jkData.shape = (10 subsamples, N patches/cubes)
			# random_cutter = N_patch array of functions, each to be applied to (ra, dec, z) of randoms
			jkData, jkWeights, error_scaling, random_cutter = jk3d.resample_data(catalog.data, catalog.samplecuts, patchside=args.patchSize, zcut=args.zCut, do_sdss=args.SDSS, do_3d=args.jk3d, cube_zdepth=args.cubeZdepth, largePi=0, bitmaskCut=args.bitmaskCut, occ_thresh=args.occ_thresh, SHIFT=SHIFT)
			print('jkData: ', jkData.shape)
			print('jkWeights: ', jkWeights.shape, '\n', jkWeights)
			print('=======================\t=======================\terror_scaling: ', error_scaling)

			# read & downsample randoms, for JK trimming
			jkrandoms = ds.read_randoms(args.Random)[:, :3]

			if args.SDSS:
				rand_colour = fits.open(args.Random)[1].data['color']
				jkrandoms = np.column_stack(( jkrandoms, rand_colour ))

			if (not args.SDSS) & SHIFT:
				random_ra = jkrandoms.T[0].copy()
				shifted_dec = jkrandoms.T[1].copy()
				G12 = (random_ra > 170) & (random_ra < 190)
				shifted_dec = np.where(G12, shifted_dec + 1, shifted_dec)
				jkrandoms[:, 1] = shifted_dec

			zlabel = ('Z_TONRY', 'z')[args.SDSS]
			sample_z = catalog.data[zlabel]
			print('INITIAL jackknife random downsampling..')
			jkrandoms = ds.downsample(jkrandoms, sample_z, 1, 20) # 1 bin, target 20x real density
			jkrandoms = jkrandoms[ (jkrandoms.T[2] >= sample_z.min()) & (jkrandoms.T[2] <= sample_z.max()) ]
			ra, dec, z = jkrandoms.T[:3]

			# merge random patch & redshift cuts
			random_cut_bool = np.ones([len(random_cutter), len(jkrandoms)])
			for i, edge_cut in enumerate(random_cutter):
				# edge_cuts: 0 = patch-cut, 1 = z-cut
				if args.jk3d:
					random_cut_bool[i] = edge_cut[0](ra, dec, z) & edge_cut[1](ra, dec, z)
				else:
					random_cut_bool[i] = edge_cut(ra, dec, z)
			Njkregions = len(jkData[0])

			skinny_patch_cuts = []
			for i,sam in enumerate(jkData):
				# filter empty patches
				skinny_patch_cut = np.array( [ ( x.shape!=(0,) ) for x in sam ], dtype=bool )
				skinny_patch_cuts.append(skinny_patch_cut)
				if i<4:
					popd_sam = sam[ skinny_patch_cut ]
				else: # if density sample
					popd_sam = sam[ skinny_patch_cuts[i-4] ] # match patch cuts to shapes for IA
					if i>=8:
						popd_sam = sam[ skinny_patch_cut ] # no direct match if not using colour-cut densities - hope new jackknife works...
					swot_sam = sam[ skinny_patch_cut ] # save all for clustering
					for j,p in enumerate(swot_sam):
						if (args.zCut==None) & ('lowZ' in catalog.labels[i]): continue
						catalog.save_swotpatches(p, catalog.labels[i], j)

				print('saving sample patches..')
				for j,p in enumerate(popd_sam):
					if (args.zCut==None) & ('lowZ' in catalog.labels[i]): continue
					if i<4: shapes=1
					else: shapes=0
					new_p,patch_z = catalog.cut_columns(p, args.H, args.flipe1, args.flipe2, args.Kneighbour, args.R0cut, shapes=shapes, mbias=args.mbias)
					pDir = catalog.save_patches(new_p, catalog.new_root, catalog.labels[i], j, 0) # save_patches returns str(patchDir)
			del jkData, popd_sam
			gc.collect()

			# make jackknife randoms (&reals) for norm (&swot)
			print('making jackknife samples..') 				# MUST feed this fn randoms & random_cutter (.shape=(Ncubes, Nrandoms)) , or will BREAK!!
			gc.collect()
			if args.densColours:
				for radians_bool, paths_key in [(1, 'all'), (0, 'swot-all')]:
					ds_jkfunc(catalog.new_root, random_cutter=random_cut_bool, empty_patches=skinny_patch_cuts, randoms=jkrandoms, radians=radians_bool, save_jks=1, jk_randoms=1, patch_str='patch', paths=paths_key, largePi=0, sdss=args.SDSS, ccut=args.cCut)
					gc.collect()
			else:
				for radians_bool, paths_key in [(1, 'dc0'), (0, 'swot-dc0')]:
					ds_jkfunc(catalog.new_root, random_cutter=random_cut_bool, empty_patches=skinny_patch_cuts, randoms=jkrandoms, radians=radians_bool, save_jks=1, jk_randoms=1, patch_str='patch', paths=paths_key, largePi=0, sdss=args.SDSS, ccut=args.cCut)
					gc.collect()
				#ds_jkfunc(catalog.new_root, random_cutter=random_cut_bool, empty_patches=skinny_patch_cuts, randoms=jkrandoms, radians=1, save_jks=1, jk_randoms=1, patch_str='patch', paths='all', largePi=0, sdss=args.SDSS, ccut=args.cCut)
				#gc.collect()
				#ds_jkfunc(catalog.new_root, random_cutter=random_cut_bool, empty_patches=skinny_patch_cuts, randoms=jkrandoms, radians=0, save_jks=1, jk_randoms=1, patch_str='patch', paths='swot-all', largePi=0, sdss=args.SDSS, ccut=args.cCut)

			jknumbers, jkn_header = [], ''
			for i, lab in enumerate(catalog.labels[:4]):
				pDir = join(catalog.new_root,lab)
				if ((args.zCut==None)&(lab.startswith('low')))|((args.cCut==None)&('Blue' in lab)):
					print('no z/colour-cut; skipping %s..'%lab)
				elif args.jackknife:
					jkWeights_pop = jkWeights[ skinny_patch_cuts[i] ]
					print('===================\t %s reduced jkWeights: '%lab, jkWeights_pop.shape)
					catalog.jackknife_patches(pDir)

					if args.unit_weights:
						print('===================\t FORCING UNIT WEIGHTS FOR JACKKNIFE \t===================')
						jkWeights_pop = np.ones_like(jkWeights_pop, dtype=float)
					np.savetxt(join(catalog.new_root, 'jkWeights_%s.txt'%lab), jkWeights_pop, header='%s\n%i samples'%(lab, len(jkWeights_pop)))

					jknumbers.append(len(jkWeights_pop))
					jkn_header += '%s\t'%lab
					if args.makejk_only:
						print('====================\t====================\t SKIPPING JACKKNIFE CORRELATIONS ====================\t====================\t')
					else:
						catalog.wcorr_jackknife(pDir, args.rpBins, args.rpLims, args.losBins, args.losLim, args.nproc, 0, args.densColours, treecorr=args.treecorr, **top_tc_config)
					catalog.jackknife(pDir, jkWeights_pop, error_scaling, 0)
			np.savetxt(join(catalog.new_root, 'JK_subsample_numbers.txt'), np.array(jknumbers), header=jkn_header, fmt='%i')

		if args.largePi:
			jkData, jkWeights, error_scaling, random_cutter = jk3d.resample_data(catalog.data, catalog.samplecuts, patchside=args.patchSize, zcut=args.zCut, do_sdss=args.SDSS, do_3d=args.jk3d, cube_zdepth=args.cubeZdepth, largePi=1, bitmaskCut=args.bitmaskCut, SHIFT=SHIFT)
			Njkregions = len(jkData[0])

			# read & downsample randoms, for JK trimming
			jkrandoms = ds.read_randoms(args.Random)[:, :3]

			if args.SDSS:
				rand_colour = fits.open(args.Random)[1].data['color']
				jkrandoms = np.column_stack(( jkrandoms, rand_colour ))

			if (not args.SDSS) & SHIFT:
				random_ra = jkrandoms.T[0].copy()
				shifted_dec = jkrandoms.T[1].copy()
				G12 = (random_ra > 170) & (random_ra < 190)
				shifted_dec = np.where(G12, shifted_dec + 1, shifted_dec)
				jkrandoms[:, 1] = shifted_dec

			zlabel = ('Z_TONRY', 'z')[args.SDSS]
			sample_z = catalog.data[zlabel]
			print('INITIAL jackknife random downsampling..')
			jkrandoms = ds.downsample(jkrandoms, sample_z, 1, 20) # 1 bin, target 20x real density
			jkrandoms = jkrandoms[ (jkrandoms.T[2] >= sample_z.min()) & (jkrandoms.T[2] <= sample_z.max()) ]
			ra, dec, z = jkrandoms.T[:3]

			# merge random patch & redshift cuts
			random_cut_bool = np.ones([len(random_cutter), len(jkrandoms)])
			for i, edge_cut in enumerate(random_cutter):
				# edge_cuts: 0 = patch-cut, 1 = z-cut
				if args.jk3d:
					random_cut_bool[i] = edge_cut[0](ra, dec, z) & edge_cut[1](ra, dec, z)
				else:
					random_cut_bool[i] = edge_cut(ra, dec, z)
			Njkregions = len(jkData[0])

			skinny_patch_cuts = []
			for i,sam in enumerate(jkData):

				skinny_patch_cut = np.array( [ ( x.shape!=(0,) ) for x in sam ], dtype=bool )
				skinny_patch_cuts.append(skinny_patch_cut)
				if (i<4) | (i>=8):
					popd_sam = sam[ skinny_patch_cut ]
				else:
					popd_sam = sam[ skinny_patch_cuts[i-4] ]

				for j,p in enumerate(popd_sam):
					if (args.zCut==None) & ('lowZ' in catalog.labels[i]): continue
					if i<4: shapes=1
					else: shapes=0
					new_p,patch_z = catalog.cut_columns(p, args.H, args.flipe1, args.flipe2, args.Kneighbour, args.R0cut, shapes=shapes, mbias=args.mbias)
					pDir = catalog.save_patches(new_p, catalog.new_root, catalog.labels[i], j, 1) # pDir (patch/cube directory) appended with _largePi
				if (3<i<8) | ((i>=8) & args.densColours): # density samples ; gen jk samples
					if (args.zCut==None) & ('lowZ' in catalog.labels[i]): continue
					catalog.jackknife_patches(pDir) # need this function call for WEIGHTS

			# no swot-files for largePi - can't set lower Pi-limit
			if args.densColours:
				ds_jkfunc(catalog.new_root, random_cutter=random_cut_bool, empty_patches=skinny_patch_cuts, randoms=jkrandoms, radians=1, save_jks=0, jk_randoms=1, patch_str='patch', paths='all', largePi=1, sdss=args.SDSS, ccut=args.cCut)
			else:
				ds_jkfunc(catalog.new_root, random_cutter=random_cut_bool, empty_patches=skinny_patch_cuts, randoms=jkrandoms, radians=1, save_jks=1, jk_randoms=1, patch_str='patch', paths='dc0', largePi=1, sdss=args.SDSS, ccut=args.cCut)

			jknumbers, jkn_header = [], ''
			for i, lab in enumerate(catalog.labels[:4]):
				pDir = join(catalog.new_root,lab+'_largePi')
				if ((args.zCut==None) & (lab.startswith('low'))) | ((args.cCut==None) & ('Blue' in lab)):
					print('no z/colour-cut; skipping %s..'%lab)
				elif args.jackknife:
					jkWeights_pop = jkWeights[ skinny_patch_cuts[i] ]
					print('===================\treduced jkWeights: ', jkWeights_pop.shape)
					catalog.jackknife_patches(pDir)
					np.savetxt(join(catalog.new_root, 'jkWeights_%s_largePi.txt'%lab), jkWeights_pop, header='%s largePi\n%i samples'%(lab, len(jkWeights_pop)))
					jknumbers.append(len(jkWeights_pop))
					jkn_header += '%s\t'%lab
					catalog.wcorr_jackknife(pDir, args.rpBins, args.rpLims, args.losBins, args.losLim, args.nproc, 1, args.densColours, treecorr=args.treecorr, **top_tc_config)
					catalog.jackknife(pDir, jkWeights_pop, error_scaling, 1)
			np.savetxt(join(catalog.new_root, 'JK_subsample_numbers_largePi.txt'), np.array(jknumbers), header=jkn_header+'\nlargePi', fmt='%i')

	catalog.make_combos(args.densColours)

	if (args.zCut==None) & (args.cCut==None):
		adjusted_combos = [catalog.wcorr_combos[0]]
	elif (args.zCut==None) & (args.cCut!=None):
		adjusted_combos = catalog.wcorr_combos[:2]
	elif (args.zCut!=None) & (args.cCut==None):
		adjusted_combos = np.array(catalog.wcorr_combos)[np.array([0,2])]
	elif (args.zCut!=None) & (args.cCut!=None):
		adjusted_combos = catalog.wcorr_combos
	print('CHECK THIS z/colour cut-adjusted combos:\n')
	for ac in adjusted_combos:
		print(ac)

	if not args.treecorr:
		catalog.prep_wcorr(catalog.new_root, adjusted_combos, args.rpBins, args.rpLims, args.losBins, args.losLim, args.largePi, args.nproc, args.play, 'real_wcorr')

		if args.wcorr:
			print('QSUBBING REAL_WCORR..')
			if args.play:
				os.system(join(catalog.new_root, 'real_wcorr.sh'))
			else:
				os.system('qsub '+ join(catalog.new_root, 'real_wcorr.sh'))
	else:
		assert args.Random != None, ("TreeCorr module requires randoms")

	if args.Random != None:
		catalog2 = RandomCatalogue(args.Random, args.densColours, args.SDSS, SHIFT=SHIFT)
		randoms_3col = np.column_stack(( catalog2.data[catalog2.headers['ra']], catalog2.data[catalog2.headers['dec']], catalog2.data[catalog2.headers['z']] ))

		for samz_k in samz_keys:
			try:
				sam_z = samz_dict[samz_k]
			except KeyError:
				print('NO KEY: %s - skipping..'%samz_k)
				continue

			if args.SDSS & ('z1' not in samz_k):
				acC = args.cCut # conditional below does not work unless renaming args.cCut...
				if args.densColours & (acC != None):
					# apply colour cut to SDSS randoms - creates much better n(z) wrt real samples
					r_redcut = np.array(catalog2.data['color'] > args.cCut)
					if samz_k.endswith('r'):
						r_3cols = randoms_3col[ r_redcut ]
					elif samz_k.endswith('b'):
						r_3cols = randoms_3col[ ~r_redcut ]
				else:
					r_3cols = randoms_3col.copy()
			else:
				r_3cols = randoms_3col.copy()

			#nbins = (1, 1)[args.SDSS] # currenty favoured: GAMA=1 (or 3?), SDSS=1(if can split red/blue)
			nbins = 1
			ds_randoms = ds.downsample(r_3cols, sam_z, nbin=nbins, target_nz=10)
			catalog2.samples.append( ds_randoms )
			catalog2.samplecounts.append( len(ds_randoms) )

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
		print('CHECK THIS rand combos:\n')
		for rc in rand_combos:
			print(rc)

		if args.treecorr:
			for i in range(len(adjusted_combos)):
				densf, dc, shapesf, sc, outf = adjusted_combos[i]
				rdensf, rdc, rshapesf, rsc, routf = rand_combos[i]
				densf = join(catalog.new_root, densf)
				shapesf = join(catalog.new_root, shapesf)
				rdensf = join(catalog.new_root, rdensf)
				rshapesf = join(catalog.new_root, rshapesf)
				outf = join(catalog.new_root, 'wcorr_'+outf+'.dat')
				# densf, shapesf, drandf, srandf, config, outfile,
				# estim='PW1', np=16, **kwargs(nbins_rpar=30, random_oversampling=10., verbosity=1, {load_RRs, save_RRs -- for PW2 only})
				catalog.run_treecorr(densf, shapesf, rdensf, rshapesf, top_tc_config['tc_config'], outf, np=args.nproc, **top_tc_config['tc3dcf_kwargs'])

			if args.plot:
				os.system(('python /share/splinter/hj/PhD/catalog_sampler.py -Catalog %s '%args.Catalog +
							'-Path %s -bootstrap %s -jackknife %s -SDSS %s '%(catalog.new_root, args.bootstrap, args.jackknife, args.SDSS) +
							'-treecorr %s -plotNow 1 -chiSqu 0' % args.treecorr))
			if args.make_shear:
				os.system('python ShearResp.py %s %s' % (('-gama', '-sdss')[args.SDSS], catalog.new_root))
				if args.largePi:
					os.system('python ShearResp.py %s %s -lpi 1' % (('-gama', '-sdss')[args.SDSS], catalog.new_root))

		elif not args.treecorr:
			catalog2.prep_wcorr(catalog.new_root, rand_combos, args.rpBins, args.rpLims, args.losBins, args.losLim, args.largePi, args.nproc, args.play, 'rand_wcorr')

			with open(join(catalog.new_root, 'rand_wcorr.sh'), 'a') as script:
				if args.plot:
					script.write('\n')
					script.write(
					'\npython /share/splinter/hj/PhD/catalog_sampler.py -Catalog %s -Path %s -plotNow 1 -chiSqu 0 -bootstrap %s -jackknife %s -SDSS %s'%(args.Catalog, catalog.new_root, args.bootstrap, args.jackknife, args.SDSS)
					)
					script.write('\n')

				if args.make_shear:
					script.write('\n')
					script.write('\npython ShearResp.py %s %s' % (('-gama', '-sdss')[args.SDSS], catalog.new_root))
					if args.largePi:
						script.write('\npython ShearResp.py %s %s -lpi 1' % (('-gama', '-sdss')[args.SDSS], catalog.new_root))
					script.write('\n')

			if args.wcorr:
				print('QSUBBING RANDOM_WCORR..')
				if args.play:
					os.system(join(catalog.new_root, 'rand_wcorr.sh'))
				else:
					os.system('qsub '+ join(catalog.new_root, 'rand_wcorr.sh'))

	import datetime
	now = datetime.datetime.now()
	print('FINISHED!')
	with open(join(catalog.new_root,'C-lineArgs_SampleProps.txt'),'a') as script:
		script.write('\n')
		script.write('\n')
		script.write(now.strftime("%Y-%m-%d %H:%M"))
		script.write('\n')
		for argk in np.sort(vars(args).keys()):
			script.write("\n%s = %s"%(argk, vars(args)[argk]))
		script.write('\n')
		[script.write('%s: \t%d, mean(R_mag - 22.): %.4f\n'%(catalog.labels[i],catalog.samplecounts[i],catalog.Rmags[i])) for i in range(len(catalog.labels[:4]))]
		[script.write('%s: \t%d\n'%(catalog.labels[i+4],catalog.samplecounts[i+4])) for i in range(len(catalog.labels[4:]))]
		if args.Random != None:
			[script.write('%s: \t%d\n'%(catalog2.labels[i],catalog2.samplecounts[i])) for i in range(len(catalog2.labels))]
	np.savetxt(join(catalog.new_root, 'R-band_pivot_deltas.txt'), np.array(catalog.Rmags), header='mean differences between sample & pivot R-band abs mag\nignore any lowZ for SDSS\n'+'\t'.join(catalog.labels[:4]), delimiter='\t', fmt='%.4f')
	np.savetxt(join(catalog.new_root, 'SamplePopulations.txt'), np.column_stack((np.array(catalog.samplecounts[:4]), np.array(catalog.samplecounts[4:8]))), header='populations of\nshapes\t\tdensity samples\nignore any lowZ for SDSS\n'+'\t'.join(catalog.labels[:4]), delimiter='\t', fmt='%i')
	pickle.dump(catalog.Sprops, open(join(catalog.new_root, 'SampleProps.p'), 'w'))








