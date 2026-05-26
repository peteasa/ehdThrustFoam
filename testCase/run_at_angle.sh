#!/bin/bash

# helper script to run the simulation updating the flow angle in 0/U
# ./run_at_angle.sh 0.1

cd "$(dirname "$0")"
echo "Running airfoil simulation at $1 degrees at location $(pwd)"

# dont forget to run (once only)
#source $HOME/OpenFOAM/OpenFOAM-13/etc/bashrc

./Allcleansimulation

foamDictionary 0/U -entry angle -set $1
./Allrun
if [[ -z "${SILENT}" ]]; then
    echo paraFoam -builtin
fi
