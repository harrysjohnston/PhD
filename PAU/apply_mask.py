# coding: utf-8
from functions import *
from astropy import wcs
import argparse
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'fnames',
		nargs='*',
		help='give individual fits-file paths for masking')
	parser.add_argument(
		'-mask',
		nargs='*',
		default=['W3.16bit.5arcs.reg2.fits'],
		help='give path(s) to mask file(s)')
	args = parser.parse_args()

	for fname in args.fnames:

		rand = fopen(fname)
		rand_radec = np.column_stack((rand['ra'], rand['dec']))
		cut = np.zeros_like(rand_radec.T[0], dtype=bool)

		for maskpath in args.mask:
			hdu = fits.open(maskpath)[0]
			mask = hdu.data
			w = wcs.WCS(hdu)
			pixcrd = w.wcs_world2pix(rand_radec, 0, ra_dec_order=1)
			pixcrd = np.array(np.round(pixcrd), dtype=int)
			c1 = (pixcrd.T[1] >= 0) & (pixcrd.T[1] < mask.shape[0])
			c0 = (pixcrd.T[0] >= 0) & (pixcrd.T[0] < mask.shape[1])
			idx = np.where(c0 & c1)[0]
			pixcrd_in_mask = pixcrd[c0 & c1]
			vals_at_rand_radec = mask[pixcrd_in_mask[:,1], pixcrd_in_mask[:,0]] 
			cut[idx] |= (vals_at_rand_radec == 0)

		new_rand = rand[cut]
		t = Table(new_rand)
		t.write(fname, overwrite=1)

