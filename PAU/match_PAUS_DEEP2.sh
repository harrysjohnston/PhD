#!/bin/tcsh

stilts -verbose tmatch2 matcher=sky params=1 join=all1 \
	in1='PAUS_KSB.fits' values1='alpha_j2000 delta_j2000' \
	in2='zcat.deep2.dr4.fits' values2='RA DEC' \
	out='PAUS_KSB.fits1' ofmt=fits
mv PAUS_KSB.fits1 PAUS_KSB.fits

