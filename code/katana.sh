#!/bin/bash

#PBS -l ncpus=16
#PBS -l mem=124gb
#PBS -l walltime=12:00:00
#PBS -j oe 
#PBS -J 1-300
#PBS -m a
#PBS -M peng.zhong@unsw.edu.au 

cd ${PBS_O_WORKDIR}

module load gsl/2.7.1 
module load gmp/6.2.1 
module load r/4.3.1

Rscript code/simulation_studies.R "id=${PBS_ARRAY_INDEX};thres=0.99;m=1000;computer=\"hpc\"" 
Rscript code/simulation_studies.R "id=${PBS_ARRAY_INDEX};thres=0.95;m=1000;computer=\"hpc\"" 
Rscript code/simulation_studies.R "id=${PBS_ARRAY_INDEX};thres=0.9;m=1000;computer=\"hpc\"" 
Rscript code/simulation_studies.R "id=${PBS_ARRAY_INDEX};thres=0.85;m=1000;computer=\"hpc\""
Rscript code/simulation_studies.R "id=${PBS_ARRAY_INDEX};thres=0.8;m=1000;computer=\"hpc\"" 