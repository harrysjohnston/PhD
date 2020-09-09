#!/bin/tcsh
# params=<match-radius/arcsec>
source ~/.login
pau
stilts tcat ifmt=ascii ofmt=fits in="`echo HenkCats2/mos_W3*cat`" out="W3_KSB.fits1"
stilts tpipe in="W3_KSB.fits1" out="W3_KSB.fits" \
	cmd='addcol raok ram' \
	cmd='addcol decok decm' \
	cmd='addcol e1ok e1c' \
	cmd='addcol e2ok e2c' \
	cmd='addcol pgok pgm'
stilts -verbose tmatch2 matcher=sky params=1 \
	in1='PAUS_cut.fits' values1='alpha_j2000 delta_j2000' \
	in2='W3_KSB.fits' values2='raok decok' \
	out='PAUS_KSB.fits_1' ofmt=fits
stilts -verbose tmatch2 matcher=exact \
	in1='PAUS_KSB.fits_1' values1='numeric_id' \
	in2='PAU_W3_zb.fits' values2='ref_id' \
	out='PAUS_KSB.fits_2' ofmt=fits
stilts -verbose tpipe in='PAUS_KSB.fits_2' \
	cmd='addcol "ref_id" "ref_id_1"' \
	cmd='addcol "qz" "qz_1"' \
	cmd='addcol "zb" "zb_1"' \
	cmd='delcols "*_1 *_2"' \
	out='PAUS_KSB.fits' ofmt=fits
python make_shear.py PAUS_KSB.fits

rm PAUS_KSB.fits_*
rm W3_KSB.fits1

