#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=wclusGAMA
#SBATCH -c 12
#SBATCH -t 7-00:00:00
#SBATCH --mem=90GB

date

source ~/.login

set tag = $1
set type = $2
set Gpath = '/share/splinter/hj/PhD/PAU/playground/GAMA'
set wclust = '/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr_clustering'

if ($type == 'all') then
	$wclust $Gpath \
		wgp_unwindowed_d1_${tag}.asc `wc -l $Gpath/wgp_unwindowed_d1_${tag}.asc | awk '{print $1-1}'` \
		wgp_unwindowed_d2_${tag}.asc `wc -l $Gpath/wgp_unwindowed_d2_${tag}.asc | awk '{print $1-1}'` \
		wgp_unwindowed_r1_${tag}.asc `wc -l $Gpath/wgp_unwindowed_r1_${tag}.asc | awk '{print $1-1}'` \
		wgp_unwindowed_r1_${tag}.asc `wc -l $Gpath/wgp_unwindowed_r1_${tag}.asc | awk '{print $1-1}'` \
		1 11 0.1 60.0 15 60.0 ps_${tag} 0 0 12

else if ($type == 'red') then
	$wclust $Gpath \
		wgp_red_unwindowed_d1_${tag}.asc `wc -l $Gpath/wgp_red_unwindowed_d1_${tag}.asc | awk '{print $1-1}'` \
		wgp_red_unwindowed_d2_${tag}.asc `wc -l $Gpath/wgp_red_unwindowed_d2_${tag}.asc | awk '{print $1-1}'` \
		wgp_red_unwindowed_r1_${tag}.asc `wc -l $Gpath/wgp_red_unwindowed_r1_${tag}.asc | awk '{print $1-1}'` \
		wgp_red_unwindowed_r1_${tag}.asc `wc -l $Gpath/wgp_red_unwindowed_r1_${tag}.asc | awk '{print $1-1}'` \
		1 11 0.1 60.0 15 60.0 ps_${tag}_rr 0 0 12

else if ($type == 'blue') then
	$wclust $Gpath \
		wgp_blue_unwindowed_d1_${tag}.asc `wc -l $Gpath/wgp_blue_unwindowed_d1_${tag}.asc | awk '{print $1-1}'` \
		wgp_blue_unwindowed_d2_${tag}.asc `wc -l $Gpath/wgp_blue_unwindowed_d2_${tag}.asc | awk '{print $1-1}'` \
		wgp_blue_unwindowed_r1_${tag}.asc `wc -l $Gpath/wgp_blue_unwindowed_r1_${tag}.asc | awk '{print $1-1}'` \
		wgp_blue_unwindowed_r1_${tag}.asc `wc -l $Gpath/wgp_blue_unwindowed_r1_${tag}.asc | awk '{print $1-1}'` \
		1 11 0.1 60.0 15 60.0 ps_${tag}_bb 0 0 12
endif

date

