This an OpenFOAM-13 custom solver for electrohydrodynamic (EHD) ionic wind thrust simulation using coupled Navier-Stokes and charge transport.  The validation case is inspired by Haofeng Xu et al. "Flight of an aeroplane with solid-state propulsion" https://www.nature.com/articles/s41586-018-0707-9.  Read the documentation for this custom solver [here](https://paracpg.gitlab.io/wikis/cfdpages.html).

## Features

- simulates corona discharge thrust generation
- target application is for electrohydrodynamic propulsion devices
- custom incompressible solver also modelling temperature and physical pressure
- Electrostatic body forces
- multiple charged species
- coupled drift-diffusion charge transport
- temperature and pressure dependent mobility and diffusion coefficients
- body-force coupling to Navier - Stokes transport
- turbulent flow
- parallel execution

## Instructions

One time per session activities
```
source ${HOME}/OpenFOAM/OpenFOAM-13/etc/bashrc
cd ehdThrustFoam
./Allwmake
```
Then prepare model mesh for the test case
```
cd testCase
./generate_gmsh_unstructured.sh
```
That's the preparation done, now change any parameters in the 0.orig, system/*, constant/* folder and run the simulation for a particular background wind direction
```
./run_at_angle.sh 10.0
```
View the results with paraFoam (see [OpenFOAM user guide](https://doc.cfd.direct/openfoam/user-guide-v13/paraview#dx39-204001)) and with gnuplot (see minmax.gp).

Once complete you can stop and restart the simulation, or re-run with different parameters set, all without changing the model.

## Actively seaking collaboration and or feedback!

You are welcome to provide feedback by responding to an existing issue or creating a new issue on the [issue board](https://gitlab.com/paracpg/ehdthrustfoam/-/boards).

I invite you to help

- extend the multi species equations adding additional plasma chemistry
- create new models for verification
- develop an equivalent compressible solver
- review and improve my documentation for the plasma physics

## Installation

Requires:

OpenFOAM-13 https://github.com/OpenFOAM/OpenFOAM-13

If you provide your own mesh then the following are optional, otherwise these are the tools I used to create the mesh

Julia scripting language to create the model https://julialang.org/

Gmsh to generate the model mesh https://gitlab.onelab.info/gmsh/gmsh.git

The development documentation for this solver is on [gitlab](https://paracpg.gitlab.io/wikis/cfdpages.html)