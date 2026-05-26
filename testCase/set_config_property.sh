#!/bin/bash

##########################################################
# helper script used by julia script ehdthrustmesh.jl
# used to save simulation settings
##########################################################

echo $(jq ".$1=\"$2\"" config.json) > config.json
