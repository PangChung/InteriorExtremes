#!/bin/bash

#PBS -l ncpus=16
#PBS -l mem=124gb
#PBS -l walltime=12:00:00
#PBS -j oe 
#PBS -J 1-100
#PBS -m ae
#PBS -M peng.zhong@unsw.edu.au 

cd ${PBS_O_WORKDIR}

module load gsl/2.7.1 
module load gmp/6.2.1 
module load r/4.3.1

Rscript code/simulation_studies.R "id=${PBS_ARRAY_INDEX};m=1000;computer=\"hpc\";model=\"logskew\";d=15;basis.idx=1" 
Rscript code/simulation_studies.R "id=${PBS_ARRAY_INDEX};m=1000;computer=\"hpc\";model=\"logskew\";d=15;basis.idx=2" 
#Rscript code/simulation_studies.R "id=${PBS_ARRAY_INDEX};m=5000;computer=\"hpc\";model=\"truncT\";d=8" 