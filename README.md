This OpenFOAM custom solver and validation case were developed for OpenFOAM-13.  The validation case is inspired by Haofeng Xu, Yiou He, Kieran L. Strobel, Christopher K. Gilmore, Sean P. Kelley, Cooper C. Hennick, Thomas Sebastian, Mark R. Woolston, David J. Perreault & Steven R. H. Barrett "Flight of an aeroplane with solid-state propulsion" https://www.nature.com/articles/s41586-018-0707-9

## Instructions

One time only activities
```
source ${HOME}/OpenFOAM/OpenFOAM-13/etc/bashrc
cd ehdThrustFoam
./Allwmake
```
Then prepare model for the test case
```
cd testCase
./generate_gmsh_unstructured.sh
```
Now change any parameters in the 0.orig, system/*, constant/* folder and run the simulation
```
./run_at_angle.sh 10.0
```
View the results with paraFoam (see [OpenFOAM user guide](https://doc.cfd.direct/openfoam/user-guide-v13/paraview#dx39-204001)) and with gnuplot (see minmax.gp).

## Installation

Requires:

OpenFOAM-13 https://github.com/OpenFOAM/OpenFOAM-13

Julia scripting language to create the model https://julialang.org/

Gmsh to generate the model mesh https://gitlab.onelab.info/gmsh/gmsh.git

