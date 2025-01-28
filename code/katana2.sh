#!/bin/bash

#PBS -l ncpus=16
#PBS -l mem=124gb
#PBS -l walltime=12:00:00
#PBS -j oe 
#PBS -J 1-4
##PBS -m ae
##PBS -M peng.zhong@unsw.edu.au 

cd ${PBS_O_WORKDIR}

module load gsl/2.7.1 
module load gmp/6.2.1 
module load r/4.3.1

# Rscript $script "computer=\"hpc\";id=${PBS_ARR2AY_INDEX}" 
# Rscript code/application_florida.R "computer=\"hpc\";idx.jack=${PBS_ARRAY_INDEX};method=\"Nelder-Mead\";id=${id}" 
# Rscript code/application_florida.R "computer=\"hpc\";idx.jack=0;method=\"Nelder-Mead\";id=${PBS_ARRAY_INDEX}" 
Rscript $script "computer=\"hpc\";idx.jack=0;method=\"Nelder-Mead\";id=${PBS_ARRAY_INDEX}" 
# Rscript code/application.R "id=${PBS_ARRAY_INDEX};computer=\"hpc\""
# Rscript code/application.R "computer=\"hpc\""
# Rscript code/simulation_studies_comp.R "id=${PBS_ARRAY_INDEX};computer=\"hpc\";d=15;m=100"
# for i in {1..300}; do
#     Rscript code/simulation_studies_comp.R "id=${i};computer=\"hpc\";d=15;m=100"
#     echo $i
# done
# Rscript code/simulation_studies_comp2_2.R "id=${PBS_ARRAY_INDEX};computer=\"hpc\";d=15;m=2000"