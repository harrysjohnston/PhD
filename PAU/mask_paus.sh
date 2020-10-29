#!/bin/tcsh
  
if ($# < 1) then
    echo "Give randoms fits file(s) to have W3 masking applied"
    exit 1
endif

python apply_mask.py $1
python $KIDS_PATH/mask_randoms.py $1 -mask W3_2048mask.fits -hardcut 1 -colnames ra dec
python extra_masking_paus.py $1


