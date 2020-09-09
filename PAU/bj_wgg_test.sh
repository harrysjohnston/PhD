#!/bin/tcsh

#/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr_clustering $PAU_PATH/OUTPUTS_specphot_test wgg_DspecRspec_d1.asc 181292 wgg_DspecRspec_d1.asc 181292 wgg_DspecRspec_r1.asc 1812222 wgg_DspecRspec_r1.asc 1812222 1 11 0.1 60 15 60 bjwclust 0 0 16

/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr $PAU_PATH/OUTPUTS_specphot_test wgp_DspecRspec_d1.asc 163828 wgp_DspecRspec_d2.asc 100672 11 0.1 60 15 60 bjwcorr 16 0 0
/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr $PAU_PATH/OUTPUTS_specphot_test wgp_DspecRspec_r1.asc 1638145 wgp_DspecRspec_d2.asc 100672 11 0.1 60 15 60 bjwcorr_rand 16 0 0


