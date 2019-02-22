#!/usr/bin/env bash
# bash script to recalculate all TD01 files including the colocation with model data

# work in data TD01: emep domain, colocate emep data
./Aeolus2Netcdf_TD01_colocate.sh

# just convert data to netcdf and plot maps and curtain
./Aeolus2Netcdf_TD01.sh