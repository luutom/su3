#!/usr/bin/env bash

#
# This is an example submit script to demonstrate how to run
# the myChroma executable
#

#SBATCH --job-name="su3_meas"
#SBATCH --account=jias41
#SBATCH --nodes=8
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --output=./out.%j
#SBATCH --error=./err.%j
#SBATCH --mail-user=t.luu@fz-juelich.de
#SBATCH --mail-type=FAIL
#SBATCH --partition=booster
##########SBATCH --cpus-per-task=68

set -e

maxConfig=4640
delta=100
let nextcfg=$1+$delta

# submit job-dependent script if possible
if [ "$nextcfg" -le "$maxConfig" ]
then
    echo "sbatch --dependency=${SLURM_JOB_ID} su3_submit_ms-0.013.sh $nextcfg"
    sbatch --dependency=${SLURM_JOB_ID} su3_submit_ms-0.013.sh $nextcfg
fi

# get all the relevant packages and so on. . .
ml Architecture/jurecabooster
source $HOME/scripts/activate_booster_chroma.sh
ml SciPy-Stack

cd $PROJECT_cjias41/luu/su3  # executable is here. . .

#parameters
beta=3.6
ms=-0.013
mc=0.25
cfg=${1}

#make directory where the contractions will be stored
direc=$PROJECT_cjias41/luu/su3/data2/b${beta}_ms${ms}_mc${mc}_cfg_${cfg}
mkdir -p ${direc}

#generate xml scripts
python3 genXML.py ${ms} ${mc} ${beta} ${cfg} ${direc} > ${direc}/input.xml

#export OMP_NUM_THREADS=68
export KMP_AFFINITY=scatter

nnodes=8
ncores=2

if [[ $ncores -eq 2 ]]; then
    nmpi=$(( 64 * nnodes ))
    export OMP_NUM_THREADS=4
    #MG_ARGS="-geom 2 4 4 16 -c 2 -by 4 -bz 4 -sy 1 -sz 2 -minct1 -pxy 0 -pxyz 0"
    MG_ARGS="-geom 2 4 4 16 -c 2 -by 4 -bz 4 -sy 1 -sz 2 -minct1 -pxy 1 -pxyz 0"
elif [[ $ncores -eq 4 ]]; then
    nmpi=$(( 32 * nnodes ))
    export OMP_NUM_THREADS=8
    MG_ARGS="-geom 1 4 4 16 -c 4 -by 4 -bz 4 -sy 1 -sz 2 -minct1 -pxy 0 -pxyz 0"
fi

srun -N$nnodes -n$nmpi --cpu-bind=cores ./myChroma_KNL $MG_ARGS -i ${direc}/input.xml -o ${direc}/output.xml -l ${direc}/log.xml > ${direc}/std.out

