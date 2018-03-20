#!/bin/bash

PARA=$(basename "$1")
SDF=$(basename "$2")
CHG=$(basename "$3")

PARA_name=${PARA:4:-4}
date=$(date +%Y_%m_%d_%H_%M_%S_$PARA_name)
ssh dargen3@tarkil.grid.cesnet.cz "cd calculator_charges ; mkdir $date" 




printf "Copying of data to MetaCentrum...\n"
scp "$1" "$2" "$3" dargen3@tarkil.grid.cesnet.cz:/storage/praha1/home/dargen3/calculator_charges/$date
printf "Copying was sucessfull.\n\n\n"
printf "Connection to MetaCentrum...\n"
ssh dargen3@tarkil.grid.cesnet.cz "cd calculator_charges/para_submit; ./para_submit.sh $PARA $SDF $CHG $date $5 \"${4}\" -l select=1:ncpus=$5:mem=$5gb:scratch_local=$5gb -l walltime=100:00:00"
printf "Job is in planning system.\n\n\n"


