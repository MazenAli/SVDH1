#!/bin/bash

### SGE options ################################################################
# use bash as shell
#$ -S /bin/bash
# use current (submit) directory as working directory
#$ -cwd

#$ -N PmQm1

#$ -q xeon_e5-2630v3_2x8.q
#$ -l h_rt=5:00:00
#$ -l h_rss=10G
#$ -l h_vmem=10G
#$ -j n

# Send notification emails to:
#$ -M mazen.ali@uni-ulm.de
# Send on beginning, end, abortion, suspension:
#$ -m abes
#
################################################################################

# init environment
. /etc/profile.d/modules.sh
module load matlab/r2016b

# File
INPUTFILE="svdPmQm.m"

# directories
BASE=$(pwd)

# debug output
echo "running on          $HOSTNAME"
echo "base is set to      $BASE"

# prepare scratch directory
cd $BASE

# run job
echo "running job:"
matlab < $INPUTFILE
echo "Job exited with status $?"

# clean up
echo "Done, moving scratch/result data to local..."

exit 0
