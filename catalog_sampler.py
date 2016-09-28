from __future__ import print_function, division
from astropy.io import fits
import numpy as np
from astropy.io import ascii
from os.path import join, isdir
from os import listdir, mkdir
import os
import argparse
import csv
from astropy import cosmology
from astropy.cosmology import Planck13
import matplotlib.pyplot as plt

class RealCatalogue:

	def __init__(self, path):
		"""""
		read-in catalogue

		"""""
		self.path = path
		hdulist = fits.open(path)
		self.data = hdulist[1].data
		self.columns = hdulist[1].columns.names
		self.labels = ['highZ_Red', 'highZ_Blue', 'lowZ_Red', 'lowZ_Blue', 'highZ', 'lowZ']

	def cut_data(self, pgm_, z_, colour_, *bitmask_): 
		"""""
		cut catalogue according to bitmasks, 
		PGM, & into subsamples
		
		"""""
		assert 'RA_1' in self.columns, "'RA_1' not in columns, see column headers: "+ str(self.columns)
		assert 'DEC_1' in self.columns, "'DEC_1' not in columns, see column headers: "+ str(self.columns)
		assert 'pgm' in self.columns, "'pgm' not in columns, see column headers: "+ str(self.columns)
		assert 'Z_1_1' in self.columns, "'Z_1_1' not in columns, see column headers: "+ str(self.columns)
		assert 'absmag_g_1' in self.columns, "'absmag_g_1' not in columns, see column headers: "+ str(self.columns)
		assert 'absmag_i_1' in self.columns, "'absmag_i_1' not in columns, see column headers: "+ str(self.columns)
		assert 'col3' in self.columns, "'col3' not in columns, see column headers: "+ str(self.columns)

		total_bitmasks = self.data['col3']
		if type(bitmask_) != None:
			bitmask_ = bitmask_[0]
			bitmask_cut = []
			for bit_test in total_bitmasks:
				for mask_ in bitmask_:
					if (mask_ & bit_test == mask_):
						bitmask_cut.append(False)
						break
					if mask_ == bitmask_[-1]:
						bitmask_cut.append(True)
			assert len(bitmask_cut) == len(total_bitmasks), "bitmask testing broken"
			bitmask_cut = np.array(bitmask_cut)
			self.data = self.data[bitmask_cut]
			print('bitmask cut: ', np.unique(bitmask_cut))

		pgm = self.data['pgm']
		pgm_cut = np.array((pgm > pgm_))
		self.data = self.data[pgm_cut]
		print('pgm cut: ', np.unique(pgm_cut))

		# Remove duplicates in RA/DEC:
		coordStrings = ['RA_1', 'DEC_1']
		for i, col in enumerate(coordStrings):
			coords = self.data[col]
			uniqCoords = np.unique(coords, return_inverse=True, return_counts=True)
			inverse = uniqCoords[1]
			count = uniqCoords[2]
			orderedCount = count[inverse]
			duplicateCut = orderedCount == 1
			self.data = self.data[duplicateCut]
			print('Removed %s duplicates in %s' % ((len(duplicateCut)-len(self.data)), col[:-2]))

		self.pre_count = len(self.data)
		z = self.data['z_1_1']
		self.data = self.data[(z >= 0.02)]
		z = self.data['z_1_1']
		self.pre_z = z
		colour = self.data['absmag_g_1'] - self.data['absmag_i_1']


		if colour_ != None:
			colour_cut = np.array((colour > colour_))
			colour_cut_r = np.invert(colour_cut)
		else:
			colour_cut = np.array([True]*len(self.data))
			colour_cut_r = colour_cut
			print('Red catalog == Blue catalog')
		print('c cut: ', colour_, np.unique(colour_cut))
		if z_ != None:
			z_cut = np.array((z > z_))
			z_cut_r = np.invert(z_cut)
		else:
			z_cut = np.array([True]*len(self.data))
			z_cut_r = z_cut
			print('highZ catalog == lowZ catalog')
		print('z cut: ', z_, np.unique(z_cut))

		self.highz_R = self.data[(z_cut & colour_cut)]
		self.highz_B = self.data[(z_cut & colour_cut_r)]
		self.lowz_R = self.data[(z_cut_r & colour_cut)]
		self.lowz_B = self.data[(z_cut_r & colour_cut_r)]
		self.highz = self.data[z_cut]
		self.lowz = self.data[z_cut_r]

		self.samplecounts = [len(self.highz_R), len(self.highz_B),
								len(self.lowz_R), len(self.lowz_B),
								len(self.highz), len(self.lowz)]

		[print('# objects %s: '%self.labels[i], v) for i, v in enumerate(self.samplecounts)]

		self.wcorr_combos = [
		[self.labels[4]+'.asc', self.samplecounts[4], self.labels[0]+'.asc', self.samplecounts[0], 'hiZ_vs_hiZ_R'],
		[self.labels[4]+'.asc', self.samplecounts[4], self.labels[1]+'.asc', self.samplecounts[1], 'hiZ_vs_hiZ_B'],
		[self.labels[5]+'.asc', self.samplecounts[5], self.labels[2]+'.asc', self.samplecounts[2], 'loZ_vs_loZ_R'],
		[self.labels[5]+'.asc', self.samplecounts[5], self.labels[3]+'.asc', self.samplecounts[3], 'loZ_vs_loZ_B']
		]

	def cut_columns(self, subsample, h): 
		"""""
		take subsample data 
		& isolate columns for wcorr

		"""""
		assert 'RA_1' in self.columns, "'RA_1' not in columns, see column headers: "+ str(self.columns)
		assert 'DEC_1' in self.columns, "'DEC_1' not in columns, see column headers: "+ str(self.columns)
		assert 'Z_1_1' in self.columns, "'Z_1_1' not in columns, see column headers: "+ str(self.columns)
		assert 'e1c' in self.columns, "'e1c' not in columns, see column headers: "+ str(self.columns)
		assert 'e2c' in self.columns, "'e2c' not in columns, see column headers: "+ str(self.columns)

		table = subsample
		RA = np.deg2rad(table['RA_1'])
		DEC = np.deg2rad(table['DEC_1'])
		Z = table['Z_1_1']
		# print('Z : ', Z)
		e1 = table['e1c']/table['pgm']
		e2 = table['e2c']/table['pgm']
		e2 *= -1 # for RA increasing leftward, c.f. x-axis increasing rightward
		e_weight = np.array([1]*len(table))

		comov = Planck13.comoving_distance(Z)
		comov *= h

		new_table = np.column_stack((RA, DEC, comov, e1, e2, e_weight))
		return new_table

	def save_tables(self, new_table, outfile_root_, label, z_cut, c_cut):
		"""""
		save subsample tables to ascii

		"""""
		if (z_cut != None) & (c_cut != None):
			outfile_root = outfile_root_ + "_z_" + str(z_cut) + "_c_" + str(c_cut)
		elif z_cut != None:
			outfile_root = outfile_root_ + "_z_" + str(z_cut)
		else:
			outfile_root = outfile_root_

		self.new_root = outfile_root

		if not isdir(outfile_root):
			mkdir(outfile_root)
		ascii.write(new_table, join(outfile_root, label + ".asc"), names=['#RA/rad', '#DEC/rad', '#comov_dist/Mpc/h', '#e1', '#e2', '#e_weight'])
		sample_no = str(label) + " # objects: " + str(len(new_table))
		# print(sample_no)
		return sample_no

	def prep_wcorr(self, files_path, wcorr_combos, rp_bins, rp_lims, los_bins, los_lim, nproc, large_pi, out_sh):

		shell_script = [
		'#!/bin/tcsh',
		'#PBS -q compute',
		'#PBS -N %s'%out_sh,
		'#PBS -l nodes=1',
		'#PBS -l walltime=120:00:00',
		'#PBS -l mem=50gb',
		'',
		'module load dev_tools/nov2014/python-anaconda',
		'',
		'cd $PBS_O_WORKDIR',
		'',
		'setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/share/splinter/hj/PhD/CosmoFisherForecast/bjutils/lib/',
		'',
		'date']

		for combo in wcorr_combos:
			shell_script.append('')
			outfile = combo[4]
			if large_pi == 1:
				# outfile = outfile[:-4]
				outfile += '_largePi'
			shell_script.append('/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr %s %s %s %s %s %s %s %s %s %s %s %s %s 0' %	(files_path, combo[0], combo[1], combo[2], combo[3], rp_bins, rp_lims[0], rp_lims[1], los_bins, los_lim, outfile, nproc, large_pi)
				)
			shell_script.append('')

		shell_script.append('date')

		File = join(files_path, '%s.sh'%out_sh)
		Write = open(File, 'w')
		Text = '\n'.join(shell_script)
		Write.write(str(Text))
		Write.close()

		wcorr_spec = []
		wcorr_spec.append('Comoving transverse separation r_p: %s - %s Mpc/h in %s log-spaced bins'%(rp_lims[0], rp_lims[1], rp_bins))
		wcorr_spec.append('Comoving line-of-sight separation \Pi: %s - %s Mpc/h in %s bins'%(los_lim*(-1), los_lim, los_bins))
		wcorr_spec.append('No. processors: %s'%nproc)
		wcorr_spec.append('Large-Pi systematics testing: %s'%large_pi)

		File = join(files_path, 'wcorr_spec')
		Write = open(File, 'w')
		Text = '\n'.join(wcorr_spec)
		Write.write(str(Text))
		Write.close()

	# def plot_wcorr(self, files_path, combos, ): # NEED TO WORK IN RANDOM-SUBTRACTION SOMEHOW
	# 	data_files = []
	# 	[data_files.append(('wcorr_' + item[4] + '.dat') for item in combos)]
	# 	f, axarr = plt.subplots(2, 2)


class RandomCatalogue(RealCatalogue):

	def __init__(self, path):
		"""""
		read-in catalogue

		"""""
		self.path = path
		hdulist = fits.open(path)
		self.data = hdulist[1].data
		self.columns = hdulist[1].columns.names
		self.labels = ['rand_highZ', 'rand_lowZ']

	def cut_data(self, z_, len_reals, z_reals): 
		"""""
		cut catalogue into redshift subsamples
		
		"""""
		assert 'RA' in self.columns, "'RA' not in columns, see column headers: "+ str(self.columns)
		assert 'DEC' in self.columns, "'DEC' not in columns, see column headers: "+ str(self.columns)
		assert 'Z' in self.columns, "'Z' not in columns, see column headers: "+ str(self.columns)
		assert 'RAND_NUM' in self.columns, "'RAND_NUM' not in columns, see column headers: "+ str(self.columns)

		# cut down (randomly) for correlation
		randNum = self.data['RAND_NUM']
		current = len(randNum)
		target = len_reals*10
		fraction = target/current
		randCut = (randNum > fraction) & (randNum <= 2*fraction)
		self.data = self.data[randCut]
		print('Randoms cut down from %s objects to %s' % (len(randNum), len(self.data)))

		# Remove duplicates in RA/DEC:
		coordStrings = ['RA', 'DEC']
		for i, col in enumerate(coordStrings):
			coords = self.data[col]
			uniqCoords = np.unique(coords, return_inverse=True, return_counts=True)
			inverse = uniqCoords[1]
			count = uniqCoords[2]
			orderedCount = count[inverse]
			duplicateCut = orderedCount == 1
			self.data = self.data[duplicateCut]
			print('Removed %s duplicates in %s' % ((len(duplicateCut)-len(self.data)), col))

		z = self.data['z']
		pre_z_cut = (z >= z_reals.min()) & (z <= z_reals.max())
		self.data = self.data[pre_z_cut]
		z = self.data['z']

		if z_ != None:
			z_cut = np.array((z > z_))
			z_cut_r = np.invert(z_cut)
		else:
			z_cut = np.array([True]*len(self.data))
			z_cut_r = z_cut
			print('highZ catalog == lowZ catalog')

		self.highz = self.data[z_cut]
		self.lowz = self.data[z_cut_r]
		self.samplecounts = [len(self.highz), len(self.lowz)]

		[print('# objects %s: '%self.labels[i], v) for i, v in enumerate(self.samplecounts)]

	def cut_columns(self, subsample, h): 
		"""""
		take subsample data 
		& isolate columns for wcorr

		"""""
		assert 'RA' in self.columns, "'RA' not in columns, see column headers: "+ str(self.columns)
		assert 'DEC' in self.columns, "'DEC' not in columns, see column headers: "+ str(self.columns)
		assert 'Z' in self.columns, "'Z' not in columns, see column headers: "+ str(self.columns)

		table = subsample
		RA = np.deg2rad(table['RA'])
		DEC = np.deg2rad(table['DEC'])
		Z = table['Z']
		e1 = 2*np.random.random(len(table))
		e2 = e1
		e2 *= -1 # for RA increasing leftward, c.f. x-axis increasing rightward
		e_weight = np.array([1]*len(table))

		comov = Planck13.comoving_distance(Z)
		comov *= h

		new_table = np.column_stack((RA, DEC, comov, e1, e2, e_weight))
		return new_table

if __name__ == "__main__":
	parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
	parser.add_argument(
		'Catalog',
		help='full path of REAL catalogue to be sampled into ascii table(s)')
	parser.add_argument(
		'Outfile_root',
		help='full path of destination directory for subsample ascii catalogues, where further directories will be created and appended with "_z_<z_cut>_c_<colour_cut>" if applicable')
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
		help='lowZ vs. highZ redshift threshold, between 0 - 0.5. Omit for no redshift cut',
		default=None)
	parser.add_argument(
		'-cCut',
		type=np.float32,
		help='red vs. blue colour threshold, between 0 - 1.4 (meaningfully, between approx 0.7 - 1.1). Omit for no colour cut',
		default=None)
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
		help='specify no. of (log-spaced) bins in comoving transverse separation r_p (Mpc/h), for measurement of density-shape correlations. Defaults to 10',
		default=10)
	parser.add_argument(
		'-rpLims',
		nargs=2,
		type=np.float32,
		help='specify upper & lower (2 args, space-separated) limit in comoving transverse separation r_p (Mpc/h), for measurement of density-shape correlations. Defaults to 0.3, 60 Mpc/h',
		default=[0.3, 60])
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
		help='specify regular (0) or large-Pi systematics tests (1), defaults to 0',
		default=0)
	parser.add_argument(
		'-wcorr',
		type=int,
		choices=[0,1],
		help='initiate wcorr density-shape correlation measurements (1) or not (0), defaults to 0',
		default=0)
	parser.add_argument(
		'-notes',
		help='notes on any changed wcorr parameters, for appendage to directory name',
		default=None)
	args = parser.parse_args()

	catalog = RealCatalogue(args.Catalog)
	catalog.cut_data(args.pgm_cut, args.zCut, args.cCut, args.bitmaskCut)
	samples = [catalog.highz_R, catalog.highz_B, 									catalog.lowz_R, catalog.lowz_B,										catalog.highz, catalog.lowz]
	# labels = ['highZ_Red', 'highZ_Blue', 'lowZ_Red', 'lowZ_Blue', 'highZ', 'lowZ']
	cuts = 'z-cut: ' + str(args.zCut) + ', colour-cut (g-i): ' + str(args.cCut)
	sample_numbers = [cuts]
	if args.notes != None:
		appendage = '_%s' % str(notes)
		outfile_root = join(args.Outfile_root, 'Wcorr%s'%appendage)
	else:
		outfile_root = join(args.Outfile_root, 'Wcorr')

	for i, sample in enumerate(samples):
		new_table = catalog.cut_columns(sample, args.H)
		sample_num = catalog.save_tables(new_table, outfile_root, catalog.labels[i], args.zCut, args.cCut)
		sample_numbers.append(sample_num)

	File = join(catalog.new_root, 'Sample_popns')
	Write = open(File, "w")
	Text = "\n".join(sample_numbers)
	Write.write(str(Text))
	Write.close()

	catalog.prep_wcorr(catalog.new_root, catalog.wcorr_combos, args.rpBins, args.rpLims, args.losBins, args.losLim, args.nproc, args.largePi, 'real_wcorr')

	if args.wcorr == 1:
		os.system('cd /share/splinter/hj/PhD/CosmoFisherForecast/obstools')
		os.system('gcc ./wcorr.c -fopenmp -lgsl -lgslcblas -lm -I../bjutils/include/ -L../bjutils/lib/ -lbjutils -O3 -Wall -o wcorr')
		os.system('cd ../..')
		os.system('qsub '+ join(catalog.new_root, 'real_wcorr.sh'))

	if args.Random != None:
		catalog2 = RandomCatalogue(args.Random)
		catalog2.cut_data(args.zCut, catalog.pre_count, catalog.pre_z)
		samples = [catalog2.highz, catalog2.lowz]
		# labels = ['highZ_rand', 'lowZ_rand']
		cuts = 'z-cut: ' + str(args.zCut)
		sample_numbers = [cuts]
		# outfile_root = join(args.Outfile_root, 'Random_subsample')

		for i, sample in enumerate(samples):
			new_table = catalog2.cut_columns(sample, args.H)
			sample_num = catalog2.save_tables(new_table, outfile_root, catalog2.labels[i], args.zCut, args.cCut)
			sample_numbers.append(sample_num)

		File = join(catalog.new_root, 'Sample_rand_popns')
		Write = open(File, "w")
		Text = "\n".join(sample_numbers)
		Write.write(str(Text))
		Write.close()

		rand_combos = [
		[catalog2.labels[0]+'.asc', catalog2.samplecounts[0], catalog.labels[0]+'.asc', catalog.samplecounts[0], 'rand_hiZ_vs_hiZ_R'],
		[catalog2.labels[0]+'.asc', catalog2.samplecounts[0], catalog.labels[1]+'.asc', catalog.samplecounts[1], 'rand_hiZ_vs_hiZ_B'],
		[catalog2.labels[1]+'.asc', catalog2.samplecounts[1], catalog.labels[2]+'.asc', catalog.samplecounts[2], 'rand_loZ_vs_loZ_R'],
		[catalog2.labels[1]+'.asc', catalog2.samplecounts[1], catalog.labels[3]+'.asc', catalog.samplecounts[3], 'rand_loZ_vs_loZ_B']
		]

		catalog2.prep_wcorr(catalog.new_root, rand_combos, args.rpBins, args.rpLims, args.losBins, args.losLim, args.nproc, args.largePi, 'rand_wcorr')

		if args.wcorr == 1:
			os.system('qsub '+ join(catalog.new_root, 'rand_wcorr.sh'))
























