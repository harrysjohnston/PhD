#!/bin/tcsh

# SET JOB-ARRAY TO 1-AT-A-TIME BEFORE RUNNING THIS
#set args = "--dependency=afterany:41080"
set args = "--reservation=hj_3"

# specphot_test_corrs.sh:
# clustering (1) or IA (2)
# output directory
# windowed randoms (1) or unwindowed (0)
# dynamic Pi-binning (1) or uniform (0)
# McNaught-style k-corrected randoms (0) or z-dependent 'uniq' (1)
# use standard method (0) or drawn z-spec method (1)"

#set suff = "_zph_5e6w"
# zph randoms
#sbatch $args specphot_test_corrs.sh 1 OUTPUTS_specphot_test_dynPiBins${suff} 0 1 0 1
#sbatch $args specphot_test_corrs.sh 2 OUTPUTS_specphot_test_dynPiBins${suff} 0 1 0 1
#sbatch $args specphot_test_corrs.sh 1 OUTPUTS_specphot_test_dynPiBins_windowed${suff} 5e6 1 0 1
#sbatch $args specphot_test_corrs.sh 2 OUTPUTS_specphot_test_dynPiBins_windowed${suff} 5e6 1 0 1
#sbatch $args specphot_test_corrs.sh 1 OUTPUTS_specphot_test_unifPiBins${suff} 0 0 0 1
#sbatch $args specphot_test_corrs.sh 2 OUTPUTS_specphot_test_unifPiBins${suff} 0 0 0 1
#sbatch $args specphot_test_corrs.sh 1 OUTPUTS_specphot_test_unifPiBins_windowed${suff} 4e6 0 0 1
#sbatch $args specphot_test_corrs.sh 2 OUTPUTS_specphot_test_unifPiBins_windowed${suff} 4e6 0 0 1

# DSPECRSPEC CONFIG HAS EXTRA CORRELATION - NEED TO EDIT SBATCH SCRIPT

set suff = ""
#srun -c1 --mem=50GB -p PREEMPT python skyknife.py DspecRspec_config.ini >& log_skyknife_DspecRspec
## standard randoms
#sbatch $args specphot_test_corrs.sh 1 OUTPUTS_specphot_test_dynPiBins${suff} 0 1 0 0
#sbatch $args specphot_test_corrs.sh 2 OUTPUTS_specphot_test_dynPiBins${suff} 0 1 0 0
#sbatch $args specphot_test_corrs.sh 1 OUTPUTS_specphot_test_dynPiBins_windowed${suff} 4e6 1 0 0
#sbatch $args specphot_test_corrs.sh 2 OUTPUTS_specphot_test_dynPiBins_windowed${suff} 4e6 1 0 0
sbatch $args specphot_test_corrs.sh 1 OUTPUTS_specphot_test_unifPiBins${suff} 0 0 0 0
sbatch $args specphot_test_corrs.sh 2 OUTPUTS_specphot_test_unifPiBins${suff} 0 0 0 0
#sbatch $args specphot_test_corrs.sh 1 OUTPUTS_specphot_test_unifPiBins_windowed${suff} 4e6 0 0 0
#sbatch $args specphot_test_corrs.sh 2 OUTPUTS_specphot_test_unifPiBins_windowed${suff} 4e6 0 0 0

