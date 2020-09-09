#!/bin/tcsh
#SBATCH -p COMPUTE
#SBATCH --job-name=paucats
#SBATCH -c 1
#SBATCH --mem=47G
#SBATCH -t 3-00:00:00

# remove duplicates and match with photo-z + Cigale catalogs
./remove_pau_duplicates.sh >& log_refresh_paucats

# define red_sequences
python make_pau_catalog.py >>& log_refresh_paucats

# match with KSB shapes catalog -- MUST EDIT z-limits HERE IF CHANGED FROM 0.1 -- 0.8
./match_PAUS_W3_KSB.sh >>& log_refresh_paucats

# match with DEEP2 catalog
./match_PAUS_DEEP2.sh >>& log_refresh_paucats


