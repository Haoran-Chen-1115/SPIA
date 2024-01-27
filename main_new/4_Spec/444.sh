#!/bin/bash
#SBATCH -o log_%1
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --gres=gpu:1
#SBATCH --qos=low

# 导入运行环境
module load matlab/R2022a 

# 生成machinefile
NCL=$1

# MPI跨节点运行
matlab -nodesktop -nosplash -nodisplay -r "NCL=$NCL;Spec_main;quit" >> log_$NCL
