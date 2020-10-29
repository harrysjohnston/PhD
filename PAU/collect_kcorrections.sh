#!/bin/tcsh
if ( $# < 4 ) then
	echo '\
1: input maggies data, with IDs\
2: column number of chosen band (base-1)\
3: output filename\
4+: k-corrected maggies files\n'
	exit 1
endif


set input = $1
set bandcol = $2
set bandcol1 = ($bandcol - 1)
set output = $3
set maggies = "$argv[4-]"

echo
echo input = $input
echo column = $bandcol
echo output = $output
echo maggies files:
echo $maggies | sed 's/ /\n/g'
echo

set zero_maggies = `echo $maggies | sed 's/ /\n/g' | head -n 1`

awk '{ print $1 }' $input > idcol
awk '{ print $'"$bandcol"' }' $zero_maggies > maggies_at_zero

rm $output
foreach mf ( $maggies )
	awk '{ print $'"$bandcol1"' }' $mf > maggies_at_z
	awk '{ print $1 }' $mf > redshiftcol
	paste idcol redshiftcol maggies_at_z maggies_at_zero | awk '{ print $1 " " $2 " " ( $3 / $4 ) }' >> $output
	echo ${mf} = done
end

# maggy_ratios give the k-correction from z (given by column) to z=0
sed -i '1s/^/ID z maggy_ratio\n/' $output

rm idcol
rm maggies_at_zero
rm maggies_at_z
rm redshiftcol


