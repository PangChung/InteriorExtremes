#!/bin/bash

#PBS -l ncpus=16
#PBS -l mem=124gb
#PBS -l walltime=12:00:00
#PBS -j oe 
##PBS -J 1-300
##PBS -m ae 
##PBS -M peng.zhong@unsw.edu.au 

module load gsl/2.7.1 
module load gmp/6.2.1 
module load r/4.3.1

Rscript code/compress_file.R

