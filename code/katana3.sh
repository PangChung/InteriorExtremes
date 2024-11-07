#!/bin/bash

#PBS -l ncpus=16
#PBS -l mem=124gb
#PBS -l walltime=12:00:00
#PBS -j oe 
#PBS -J 1-300
#PBS -m ae
#PBS -M peng.zhong@unsw.edu.au 

cd ${PBS_O_WORKDIR}

module load gsl/2.7.1 
module load gmp/6.2.1 
module load r/4.3.1

# Rscript $script
# Rscript code/application_florida.R "computer=\"hpc\";idx.jack=${PBS_ARRAY_INDEX};method=\"Nelder-Mead\";id=${id}" 
# Rscript $script "computer=\"hpc\";idx.jack=${PBS_ARRAY_INDEX};method=\"L-BFGS-B\";id=${id}" 
# Rscript $script "${inputs};id=${PBS_ARRAY_INDEX};computer=\"hpc\"" 
# Rscript code/simulation_studies_comp.R "id=${PBS_ARRAY_INDEX};computer=\"hpc\";d=15;m=500"
# Rscript code/simulation_studies_comp2.R "id=${PBS_ARRAY_INDEX};computer=\"hpc\";d=15;m=500"
Rscript code/simulation_studies_comp2_2.R "id=${PBS_ARRAY_INDEX};computer=\"hpc\";d=15;m=500"