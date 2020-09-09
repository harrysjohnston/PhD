#!/bin/tcsh

set dir = TEST_corrs
python jk_wclust.py $dir/wgg_test_d1.fits $dir/wgg_test_d1.fits $dir/wgg_test_r1.fits $dir/wgg_test_r1.fits 0 'tag' 6 50


