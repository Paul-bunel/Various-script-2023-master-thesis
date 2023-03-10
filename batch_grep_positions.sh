#!/bin/bash -l
#SBATCH -A snic2022-22-1048

srun grep_positions.sh ../AfricanNeo_and_Public_DB_minN10_Jan2022.map
