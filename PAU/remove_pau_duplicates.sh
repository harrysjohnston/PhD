#!/bin/tcsh

stilts -verbose tpipe \
	in='PAUS.fits' out='PAUS.fits1' ofmt=fits \
	cmd='uniq ref_id'
stilts -verbose tmatch2 matcher=exact \
	in1='PAUS.fits1' values1='ref_id' \
	in2='PAU_W3_zb.fits' values2='ref_id' \
	out='PAUS.fits2' ofmt=fits
stilts -verbose tpipe \
	in='PAUS.fits2' out='PAUS.fits3' ofmt=fits \
	cmd='addcol "ref_id" "ref_id_1"' \
	cmd='delcols "ref_id_1 ref_id_2"'
stilts -verbose tmatch2 matcher=exact \
    in1='PAUS.fits3' values1='ref_id' \
    in2='W3catalog.fits' values2='ref_id' \
    out='PAUS.fits4'
stilts -verbose tpipe in='PAUS.fits4' \
    cmd='addcol "ref_id" "ref_id_1"' \
    cmd='addcol "lp_mu" "lp_mu_1"' \
    cmd='addcol "lp_mg" "lp_mg_1"' \
    cmd='addcol "lp_mr" "lp_mr_1"' \
    cmd='addcol "lp_mi" "lp_mi_1"' \
    cmd='addcol "bcnz_zb" "bcnz_zb_1"' \
    cmd='addcol "class3" "class_3"' \
    cmd='addcol "class2" "class_2"' \
    cmd='delcols "*_1 *_2 class_3"' \
    out='PAUS.fits5' ofmt=fits
stilts -verbose tmatch2 matcher=exact suffix1='' suffix2='_normi' join=all1 \
	in1='PAUS.fits5' values1='ref_id' \
	in2='3class_inorm.txt' values2='id' ifmt2='ascii' \
	ocmd='delcols id' \
	ocmd='replacecol class3_normi "NULL_class3_normi ? class3 : class3_normi"' \
	out='PAUS_wcolumns.fits'

rm PAUS.fits[12345]

