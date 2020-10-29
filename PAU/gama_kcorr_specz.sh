#!/bin/tcsh
source ~/.login
pau

if ( $# < 2 ) then
	echo '\
1: do kcorrect\
2: collect kcorrections\n'
	exit 1
endif

set kcorrect = $1
set collect = $2

set cat = "$PHD_PATH/SMLambdarApMatchedPhotom.fits"
set zmin = 0.0
set zmax = 0.65
# for k-corrections
set magcols = 'MAG_AUTO_u MAG_AUTO_g MAG_AUTO_r MAG_AUTO_i MAG_AUTO_z'
set magerrs = 'MAGERR_AUTO_u MAGERR_AUTO_g MAGERR_AUTO_r MAGERR_AUTO_i MAGERR_AUTO_z'
#set magcols = 'obsmag_u obsmag_g obsmag_r obsmag_i obsmag_z'
#set magerrs = 'delobsmag_u delobsmag_g delobsmag_r delobsmag_i delobsmag_z'
set bandcol = 5
set zstep = 0.025
# for zmax calculation & cloning
set appmagcols = 'PETROMAG_R'
set maglims = 19.8
set z = 'Z_TONRY'
set id = 'CATAID'

set maggies = `echo $cat | sed 's/.fits/.maggies.dat/g'`
set kcorrs_out = `echo $cat | sed 's/.fits/.kcorrs/g'`
set recmaggies = "`echo $maggies | sed 's/.dat/_z\*/g'`"

echo
echo galaxy catalogue = 		$cat
echo observed magnitudes = 		$magcols
echo obs. magnitude errors = 	$magerrs
echo magnitude limit = 			$maglims
echo detection-band magnitudes =	$appmagcols
echo redshifts = 				$z
echo IDs = 						$id
echo minimum redshift = 		$zmin
echo maximum redshift = 		$zmax
echo step in redshift = 		$zstep
echo maggies column \# = 		$bandcol
echo maggies @ 					$maggies
echo kcorrs @ 					$kcorrs_out
echo

if ( $kcorrect == 1) then
	python get_maggies.py $cat -magcols $magcols -magerrs $magerrs -z $z -id $id
	./parallel_shift_maggies.sh $maggies $zmin $zmax $zstep
endif

if ( $collect == 1) then
	./collect_kcorrections.sh $maggies $bandcol $kcorrs_out $recmaggies
	python convert_kcorrs_to_h5.py $kcorrs_out
endif



