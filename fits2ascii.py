from __future__ import print_function, division
from astropy.io import fits
import numpy as np
import csv
from astropy import cosmology
from astropy.cosmology import Planck13
from astropy.io import ascii

def main(catalogue, colnames):
	# read table
	hdulist = fits.open(catalogue)
	data = hdulist[1].data
	columns = hdulist[1].columns.names

	# cut data
	if 'e1c' in columns:
		pgm_cut = data['pgm'] > 0.1
		data = data[pgm_cut]

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
		e1 = 2*np.random.random(len(data)) - 1
		e2 = 2*np.random.random(len(data)) - 1
	e_weight = np.array([1]*len(data))

	# stack columns & save as ascii table
	table = np.column_stack((RA,DEC,comov,e1,e2,e_weight))
	out_cat = catalogue[:-4] + 'ascii'
	ascii.write(table, out_cat, names=[
		'RA/rad','DEC/rad','comov_dist/h^{-1}Mpc','e1','e2','e_weight'])

	print('# objects = ', len(data))
	return None

if __name__ == "__main__":
	catalogue = '/share/data1/kids/catalogues/randoms/RandomsWindowedV01.fits'
	colnames = ['RA', 'DEC', 'Z']
	main(catalogue, colnames)
