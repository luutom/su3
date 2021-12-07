#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --account=jjsc37
#SBATCH --partition=dc-gpu
#SBATCH --gres=gpu:4
BASE=`pwd -P`


if [ $# -ne 2 ]; 
  then 
  echo "Need iconf nsrc_step"
  exit 1
fi

. env-JURDC-LLVM11-CUDA11-PS.sh

THIS_SCRIPT=job.3bar_GPU_nsrc_mom.sh
LIMIT=`cat conf.list | wc -l`
BASE=`pwd -P`
n=$1
CONF=`head -n ${n} conf.list | tail -n 1`

nsrc=`tail -n 1 ${BASE}/progress/${CONF}.txt`
echo ${nsrc}
if [ -z "${nsrc}" ]
then
  nsrc=0
fi

nsrc=$((nsrc+1))
echo ${nsrc}
cd ${BASE}
nsrc_step=$2
XML=`./xml_3bar_GPU_nsrc_mom.sh ${n} ${nsrc} $nsrc_step`
CUDA_VISIBLE_DEVICES=0,1,2,3 srun --cpu-bind=map_ldoms:5,7,1,3 --output=3bar_GPU_${n}_${nsrc}.out -n 16 ./flowChroma.gpu -i ${XML} -geom 1 1 2 8

echo $((nsrc+nsrc_step-1)) >> ${BASE}/progress/${CONF}.txt

#n=$((n+30))
#if [ "$n" -le "$LIMIT" ]
#then
#  sbatch ${THIS_SCRIPT} ${n} ${nsrc_step}
#fi
