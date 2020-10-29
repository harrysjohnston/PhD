#!/bin/tcsh
#PBS -q compute
#PBS -N PAU_ia
#PBS -l nodes=1:ppn=16
#PBS -l mem=120gb
#PBS -l walltime=72:00:00
#PBS -M zcaphjo@ucl.ac.uk
#PBS -m a

source ~/.login

setenv NP 16

setenv TC_PW1 "$PHD_PATH/TC_top_config_IA_PW1.ini"
set cs = "$PHD_PATH/catalog_sampler.py"
set rc = "$PHD_PATH/run_clustering.py"
set temp = "@$PAU_PATH/CSconfig_PAUtemplate.ini"

python $cs $temp -Catalog $PAU_PATH/OUTPUTS_flipg1/wgp_red_pau_d2.fits      -Random $PAU_PATH/OUTPUTS_flipg1/wgp_red_pau_r2.fits -notes PAUred -jackknife 0 -jk3d 0 -nproc ${NP} #-treecorr $TC_PW1
python $cs $temp -Catalog $PAU_PATH/OUTPUTS_flipg1/wgp_blue_pau_d2.fits      -Random $PAU_PATH/OUTPUTS_flipg1/wgp_blue_pau_r2.fits -notes PAUblue -jackknife 0 -jk3d 0 -nproc ${NP} #-treecorr $TC_PW1

#python $cs $temp -Catalog $PAU_PATH/PAU_BCNZ_W3KSB.fits      -Random $PAU_PATH/PAU_BCNZ_W3KSB_UniformRandoms.fits -notes PAUtotal -jackknife 1 -jk3d 0 -nproc ${NP} -treecorr $TC_PW1
#python $cs $temp -Catalog $PAU_PATH/PAU_BCNZ_W3KSB_red.fits  -Random $PAU_PATH/PAU_BCNZ_W3KSB_UniformRandoms_red.fits -notes PAUred -jackknife 1   -jk3d 0 -nproc ${NP} -treecorr $TC_PW1
#python $cs $temp -Catalog $PAU_PATH/PAU_BCNZ_W3KSB_blue.fits -Random $PAU_PATH/PAU_BCNZ_W3KSB_UniformRandoms_blue.fits -notes PAUblue -jackknife 1  -jk3d 0 -nproc ${NP} -treecorr $TC_PW1

#python $rc $PAU_PATH/Wcorr_PAUtotal_allz/ highZ_Red_UnMasked.asc -main 1 -swot 1 -swotcov 1 -rpBins 12 -rpLims 0.01 60
#python $rc $PAU_PATH/Wcorr_PAUred_allz/ highZ_Red_UnMasked.asc -main 1 -swot 1 -swotcov 1 -rpBins 12 -rpLims 0.01 60
#python $rc $PAU_PATH/Wcorr_PAUblue_allz/ highZ_Red_UnMasked.asc -main 1 -swot 1 -swotcov 1 -rpBins 12 -rpLims 0.01 60


