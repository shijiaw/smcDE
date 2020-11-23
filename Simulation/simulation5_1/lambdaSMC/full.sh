#!/bin/bash 
#SBATCH -t 2-20:00
#SBATCH --account=def-liang-ab
#SBATCH --mem-per-cpu=8G      # memory; default unit is megabytes


cd /project/6003576/shijia57/ODE/newLambdaKnots2/lambdaSMC/SMCFullB1

#echo "Using R:"
module load r/3.5.0



#echo "Starting R at `date`."
R --vanilla < main.R  
#echo "Finished R  at `date`."