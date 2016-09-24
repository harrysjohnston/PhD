from __future__ import print_function, division
from astropy.io import fits
import numpy as np
from astropy.io import ascii
from os.path import join, isdir
from os import listdir, mkdir
import argparse
import csv
from astropy import cosmology
from astropy.cosmology import Planck13

class RealCatalogue:

	def __init__(self, path):
		"""""
		read-in catalogue

		"""""
		self.path = path
		hdulist = fits.open(path)
		self.data = hdulist[1].data
		self.columns = hdulist[1].columns.names

	def cut_data(self, pgm_, z_, colour_, *bitmask_): 
		"""""
		cut catalogue according to bitmasks, 
		PGM, & into subsamples
		
		"""""
		assert 'pgm' in self.columns, "'pgm' not in columns, see column headers: "+ str(self.columns)
		assert 'Z_1_1' in self.columns, "'Z_1_1' not in columns, see column headers: "+ str(self.columns)
		assert 'absmag_g_1' in self.columns, "'absmag_g_1' not in columns, see column headers: "+ str(self.columns)
		assert 'absmag_i_1' in self.columns, "'absmag_i_1' not in columns, see column headers: "+ str(self.columns)
		assert 'col3' in self.columns, "'col3' not in columns, see column headers: "+ str(self.columns)


		# if type(bitmask_) != None:
		# 	bitmask_ = bitmask_[0]
		# 	total_bitmasks = self.data['col3']
		# 	bitmask_cut = []
		# 	for bit_test in total_bitmasks:
		# 		for mask_ in bitmask_:
		# 			if (mask_ & bit_test == mask_):
		# 				bitmask_cut.append(False)
		# 				break
		# 			if mask_ == bitmask_[-1]:
		# 				bitmask_cut.append(True)
		# 	assert len(bitmask_cut) == len(total_bitmasks), "bitmask testing broken"
		# 	print('bitmask cut: ', np.unique(bitmask_cut))

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

		z = self.data['z_1_1']
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

		print('len highz_r: ', len(self.highz_R))
		print('len highz_b: ', len(self.highz_B))
		print('len lowz_r: ', len(self.lowz_R))
		print('len lowz_b: ', len(self.lowz_B))
		

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
		print('Z : ', Z)
		e1 = table['e1c']/table['pgm']
		e2 = table['e2c']/table['pgm']
		e2 *= -1 # for RA increasing leftward, c.f. x-axis increasing rightward
		e_weight = np.array([1]*len(table))

		comov = Planck13.comoving_distance(Z)
		comov *= h

		new_table = np.column_stack((RA, DEC, comov, e1, e2, e_weight))
		return new_table

	def save_tables(self, new_table, outfile_root, label):
		"""""
		save subsample tables to ascii

		"""""
		if not isdir(outfile_root):
			mkdir(outfile_root)
		ascii.write(new_table, join(outfile_root, label + ".asc"), names=['#RA/rad', '#DEC/rad', '#comov_dist/Mpc/h', '#e1', '#e2', '#e_weight'])
		sample_no = str(label) + " # objects: " + str(len(new_table))
		print(sample_no)
		return sample_no

class RandomCatalogue(RealCatalogue):

	def __init__(self, path):
		"""""
		read-in catalogue

		"""""
		self.path = path
		hdulist = fits.open(path)
		self.data = hdulist[1].data
		self.columns = hdulist[1].columns.names

	def cut_data(self, z_): 
		"""""
		cut catalogue into redshift subsamples
		
		"""""
		assert 'Z' in self.columns, "'Z' not in columns, see column headers: "+ str(self.columns)
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
		help='full path of catalogue to be sampled into ascii table(s).')
	parser.add_argument(
		'Type',
		choices=['real', 'random'],
		help="Type of catalogue to be converted to ascii; 'real' or 'random'.")
	parser.add_argument(
		'-pgm_cut',
		type=np.float32,
		help='Shear polarizability cut, defaults to 0.1',
		default=0.1)
	parser.add_argument(
		'-z_cut',
		type=np.float32,
		help='lowZ vs. highZ redshift threshold, between 0 - 0.5. Omit for no redshift cut.',
		default=None)
	parser.add_argument(
		'-c_cut',
		type=np.float32,
		help='Red vs. Blue colour threshold, between 0 - 1.4 (meaningfully, between approx 0.7 - 1.1). Omit for no colour cut.',
		default=None)
	parser.add_argument(
		'-bitmask_cut',
		nargs='*',
		type=int,
		help='Array of bitmask IDs (powers of 2) to exclude from catalogue, eg. [4, 16, 4096,...].',
		default=None)
	parser.add_argument(
		'-_h',
		type=list,
		help='reduced Planck constant, defaults to 0.7',
		default=0.7)
	parser.add_argument(
		'Outfile_root',
		help='full path of destination directory for subsample ascii catalogues.')
	args = parser.parse_args()

	if args.Type == 'real':
		catalog = RealCatalogue(args.Catalog)
		catalog.cut_data(args.pgm_cut, args.z_cut, args.c_cut, args.bitmask_cut)
		samples = [catalog.highz_R, catalog.highz_B, 									catalog.lowz_R, catalog.lowz_B,										catalog.highz, catalog.lowz]
		labels = ['highZ_Red', 'highZ_Blue', 'lowZ_Red', 'lowZ_Blue', 'highZ', 'lowZ']
		sample_numbers = []

		for i, sample in enumerate(samples):
			new_table = catalog.cut_columns(sample, args._h)
			sample_num = catalog.save_tables(new_table, args.Outfile_root, (labels[i]))
			sample_numbers.append(sample_num)

		File = join(args.Outfile_root, 'Sample_numbers')
		Write = open(File, "w")
		Text = "\n".join(sample_numbers)
		Write.write(str(Text))
		Write.close()

	if args.Type == 'random':
		catalog = RandomCatalogue(args.Catalog)
		catalog.cut_data(args.z_cut)
		samples = [catalog.highz, catalog.lowz]
		labels = ['highZ', 'lowZ']
		sample_numbers = []

		for i, sample in enumerate(samples):
			new_table = catalog.cut_columns(sample, args._h)
			sample_num = catalog.save_tables(new_table, args.Outfile_root, (labels[i]))
			sample_numbers.append(sample_num)

		File = join(args.Outfile_root, 'Sample_numbers')
		Write = open(File, "w")
		Text = "\n".join(sample_numbers)
		Write.write(str(Text))
		Write.close()






















