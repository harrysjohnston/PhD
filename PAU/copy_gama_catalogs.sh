#!/bin/tcsh

foreach x ($argv)
	cp ../Wcorr${x}fullGAMA_rb_c0.66/highZ_Red.fits playground/GAMA/gama_red_d2${x}.fits
	cp ../Wcorr${x}fullGAMA_rb_c0.66/highZ_Red_UnMasked.fits playground/GAMA/gama_red_d1${x}.fits
	cp ../Wcorr${x}fullGAMA_rb_c0.66/rand_highZ_Red.fits playground/GAMA/gama_red_r1${x}.fits

	cp ../Wcorr${x}fullGAMA_rb_c0.66/highZ_Blue.fits playground/GAMA/gama_blue_d2${x}.fits
	cp ../Wcorr${x}fullGAMA_rb_c0.66/highZ_Blue_UnMasked.fits playground/GAMA/gama_blue_d1${x}.fits
	cp ../Wcorr${x}fullGAMA_rb_c0.66/rand_highZ_Blue.fits playground/GAMA/gama_blue_r1${x}.fits

	python rename_GAMA_randoms_columns.py

	cp playground/GAMA/gama_red_r1${x}.fits playground/GAMA/gama_red_r2${x}.fits
	cp playground/GAMA/gama_blue_r1${x}.fits playground/GAMA/gama_blue_r2${x}.fits
end

