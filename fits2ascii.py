from __future__ import print_function, division
from astropy.io import fits
import numpy as np
import csv
from astropy import cosmology
from astropy.cosmology import Planck13
from astropy.io import ascii
from os.path import join, isdir
from os import listdir, mkdir
import argparse

def main(catalog, colnames, step, outdir):
	# read table
	hdulist = fits.open(catalog)
	data = hdulist[1].data
	columns = hdulist[1].columns.names

	# cut data
	if 'e1c' in columns:
		pgm_cut = data['pgm'] > 0.1
		data = data[pgm_cut]
	else:
		rand_cut = (data['RAND_NUM'] > 0.01) & (data['RAND_NUM'] <= 0.02)
		data = data[rand_cut]
		x = len(data)//10
		if step == 10:
			step *= x
			data = data[step-x:]
		else:
			step *= x
			data = data[step-x:step]

	# compute comoving distances
	z_colname = colnames[2]
	z_col = data[z_colname]
	comov = Planck13.comoving_distance(z_col)
	h = 0.7
	comov *= h

	# isolate columns
	RA = np.array(data[str(colnames[0])])
	DEC = np.array(data[str(colnames[1])])
	RA = np.deg2rad(RA)
	DEC = np.deg2rad(DEC)
	if 'e1c' in columns:
		pgm = data['pgm']
		e1 = np.array(data[str(colnames[3])])/pgm
		e2 = np.array(data[str(colnames[4])])/pgm
	else:
		e1 = 2*np.random.random(len(data)) #- 1
		e2 = 2*np.random.random(len(data)) #- 1
	e_weight = np.array([1]*len(data))

	# stack columns & save as ascii table
	table = np.column_stack((RA,DEC,comov,e1,e2,e_weight))
	for i in np.arange(len(table[0])):
		table = table[table[:,i]!=0.]

	if 'e1c' in columns:
		ascii.write(table, outdir, names=['#RA/rad', '#DEC/rad', '#comov_dist/Mpc/h', '#e1', '#e2', '#e_weight'])
	else:
		if not isdir(outdir):
			mkdir(outdir)
		
		outfile = join(outdir, "step" + str(step) + "_rand.ascii")
		ascii.write(table, outfile, names=['#RA/rad', '#DEC/rad', '#comov_dist/Mpc/h', '#e1', '#e2', '#e_weight'])

	print('# objects = ', len(data))
	return None

if __name__ == "__main__":
	parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
	parser.add_argument(
		'catalog',
		help='path/catalog.fits to be converted to ascii')
	parser.add_argument(
		'colnames',
		nargs=5,
		help='[RA, DEC, redshift, e1, e2] column headers')
	parser.add_argument(
		'--step',
		type=int,
		choices=range(1, 11),
		help='int from 1-10, for parallelisation',
		default=None)
	parser.add_argument(
		'outdir',
		help='complete path of destination directory for random steps, or outfilename for real catalog')
	args = parser.parse_args()
	# catalog = '/share/data1/kids/catalogs/randoms/RandomsWindowedV01.fits'
	# colnames = ['RA', 'DEC', 'Z']
	# outdir = '/share/splinter/ug_hj/PhD/RandomsWindowedV01.ascii'
	main(args.catalog, args.colnames, args.step, args.outdir)
