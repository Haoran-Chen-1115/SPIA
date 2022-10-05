#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH --partition=GPU36
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --gres=gpu:4
#SBATCH --qos=low

# 导入运行环境
module load matlab/R2019a 

# 生成machinefile
NCL=$1
srun hostname -s | sort -n > slurm_${NCL}.hosts

# MPI跨节点运行
matlab -nodesktop -nosplash -nodisplay -r "NCL=$NCL;Tbar_ND;quit" > log_$NCL
