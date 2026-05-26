#!/bin/bash

##########################################################
# create the geometry (one time only)
# ./generate_gmsh_unstructured.sh
##########################################################

# dont forget to run (once only)
# source $HOME/OpenFOAM/OpenFOAM-13/etc/bashrc

MESH=mesh
MESHNAME=ehdthrust
MESHGEN=$MESHNAME$MESH

julia $MESHGEN.jl
gmsh -3 $MESHNAME.geo -format msh2
gmshToFoam $MESHNAME.msh

function setBoundaryType {
    foamDictionary constant/polyMesh/boundary -entry entry0/$1/type -set "$2"
    if [[ "$3" != "" ]]; then
        foamDictionary constant/polyMesh/boundary -entry entry0/$1/inGroups -set "$3"
    fi
}

if [[ "$MESHNAME" = "ehdthrust" ]]; then
    setBoundaryType frontAndBackPlanes empty
    setBoundaryType INLET wall
    setBoundaryType OUTLET wall
    setBoundaryType WALL wall
    setBoundaryType PELEMENT wall
    setBoundaryType NELEMENT wall
fi

