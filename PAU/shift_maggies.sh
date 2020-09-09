#!/bin/tcsh

if ( $# < 2 ) then
	echo '\
1: maggies file (ID-ordered)\
2: move galaxies to this redshift\
3: optional, shift band-pass to this redshift\n'
	exit 1
endif

set maggies = $1
set redshift = $2
if ( $# > 2 ) then
	set bpshift = $3
else
	set bpshift = 0.0
endif

set out = `echo $maggies | sed "s/\.dat/_z${redshift}\.bs${bpshift}&/g"`

echo
echo data = $maggies
echo to redshift = $redshift
echo with band-pass shift = $bpshift
echo outfile = $out
echo

# remove ID column for kcorrect call
awk '{ $1=""; print substr($0,2) }' $maggies > ${maggies}_${redshift}

cat ${maggies}_${redshift} | fit_coeffs | reconstruct_maggies --redshift $redshift --band-shift $bpshift > ! $out
rm ${maggies}_${redshift}


