#!/bin/bash
#SBATCH --ntasks 4
#SBATCH --time 48:0:0
#SBATCH --qos bbdefault
#SBATCH --mail-type ALL

set -e

module purge; module load bluebear
module load apps/matlab/r2017a

cd /rds/projects/2017/gallagmt-rapid-sperm-capture/ALHM/relaxingTests/relaxingTest_intDynamics
matlab -nodisplay -r relaxingTest_intDynamics