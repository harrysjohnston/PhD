#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=wclusPAU
#SBATCH -c 16
#SBATCH -t 7-00:00:00
#SBATCH --mem=120GB

date

source ~/.login

set tag = $1
set type = $2
set Gpath = '/share/splinter/hj/PhD/PAU/playground/'
set wclust = '/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr_clustering'

if ($type == 'all') then
	$wclust $Gpath \
		wgp_unwindowed_d1_${tag}.asc `wc -l $Gpath/wgp_unwindowed_d1_${tag}.asc | awk '{print $1-1}'` \
		wgp_unwindowed_d2_${tag}.asc `wc -l $Gpath/wgp_unwindowed_d2_${tag}.asc | awk '{print $1-1}'` \
		wgp_unwindowed_r1_${tag}.asc `wc -l $Gpath/wgp_unwindowed_r1_${tag}.asc | awk '{print $1-1}'` \
		wgp_unwindowed_r1_${tag}.asc `wc -l $Gpath/wgp_unwindowed_r1_${tag}.asc | awk '{print $1-1}'` \
		1 6 0.1 50.0 15 60.0 ps_${tag} 0 0 16

else if ($type == 'red') then
	$wclust $Gpath \
		wgp_red_unwindowed_d1_${tag}.asc `wc -l $Gpath/wgp_red_unwindowed_d1_${tag}.asc | awk '{print $1-1}'` \
		wgp_red_unwindowed_d2_${tag}.asc `wc -l $Gpath/wgp_red_unwindowed_d2_${tag}.asc | awk '{print $1-1}'` \
		wgp_red_unwindowed_r1_${tag}.asc `wc -l $Gpath/wgp_red_unwindowed_r1_${tag}.asc | awk '{print $1-1}'` \
		wgp_red_unwindowed_r1_${tag}.asc `wc -l $Gpath/wgp_red_unwindowed_r1_${tag}.asc | awk '{print $1-1}'` \
		1 6 0.1 50.0 15 60.0 ps_${tag}_rr 0 0 16

else if ($type == 'blue') then
	$wclust $Gpath \
		wgp_blue_unwindowed_d1_${tag}.asc `wc -l $Gpath/wgp_blue_unwindowed_d1_${tag}.asc | awk '{print $1-1}'` \
		wgp_blue_unwindowed_d2_${tag}.asc `wc -l $Gpath/wgp_blue_unwindowed_d2_${tag}.asc | awk '{print $1-1}'` \
		wgp_blue_unwindowed_r1_${tag}.asc `wc -l $Gpath/wgp_blue_unwindowed_r1_${tag}.asc | awk '{print $1-1}'` \
		wgp_blue_unwindowed_r1_${tag}.asc `wc -l $Gpath/wgp_blue_unwindowed_r1_${tag}.asc | awk '{print $1-1}'` \
		1 6 0.1 50.0 15 60.0 ps_${tag}_bb 0 0 16
endif

date

