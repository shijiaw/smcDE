#!/bin/bash 
#SBATCH -t 0-23:00
#SBATCH --account=def-liang-ab
#SBATCH --mem-per-cpu=8G      # memory; default unit is megabytes


cd /project/6003576/shijia57/ODE/sigma_DDE/01

#echo "Using R:"
module load r/3.5.0



#echo "Starting R at `date`."
R --vanilla < main.R  
#echo "Finished R  at `date`."