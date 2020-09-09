#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=wcorrPAU
#SBATCH -c 16
#SBATCH -t 7-00:00:00
#SBATCH --mem=120GB

date

source ~/.login

set tag = $1

# data
/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr /share/splinter/hj/PhD/PAU/playground wgp_unwindowed_d1_${tag}.asc 182583 wgp_unwindowed_d2_${tag}.asc 182583 11 0.1 60.0 15 60.0 ps_${tag} 16 0 0
/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr /share/splinter/hj/PhD/PAU/playground wgp_blue_unwindowed_d1_${tag}.asc 182583 wgp_blue_unwindowed_d2_${tag}.asc 135785 11 0.1 60.0 15 60.0 ps_${tag}_bb 16 0 0
/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr /share/splinter/hj/PhD/PAU/playground wgp_red_unwindowed_d1_${tag}.asc 182583 wgp_red_unwindowed_d2_${tag}.asc 46798 11 0.1 60.0 15 60.0 ps_${tag}_rr 16 0 0

# rand
/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr /share/splinter/hj/PhD/PAU/playground wgp_unwindowed_r1_${tag}.asc 1826493 wgp_unwindowed_d2_${tag}.asc 182583 11 0.1 60.0 15 60.0 rand_ps_${tag} 16 0 0
/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr /share/splinter/hj/PhD/PAU/playground wgp_blue_unwindowed_r1_${tag}.asc 1825926 wgp_blue_unwindowed_d2_${tag}.asc 135785 11 0.1 60.0 15 60.0 rand_ps_${tag}_bb 16 0 0
/share/splinter/hj/PhD/CosmoFisherForecast/obstools/wcorr /share/splinter/hj/PhD/PAU/playground wgp_red_unwindowed_r1_${tag}.asc 1825651 wgp_red_unwindowed_d2_${tag}.asc 46798 11 0.1 60.0 15 60.0 rand_ps_${tag}_rr 16 0 0


date
