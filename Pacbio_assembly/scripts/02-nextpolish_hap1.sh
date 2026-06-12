#!/bin/bash
#SBATCH -J polish1
#SBATCH -o log/%x.o%j
#SBATCH -e log/%x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 128


nextPolish nextpolish_hap1.cfg

