#!/bin/bash
#$ -l arch=intel*
#$ -l mem=12G
#$ -l rmem=12G
##$ -P mhd

module load apps/idl/8.5

#source /data/cs1mkg/tools/vapor/vapor-2.2.2/bin/vapor-setup.sh
#export IDL_STARTUP=runconvert_wash
export IDL_STARTUP=runconvert_wash_vtk

idl
