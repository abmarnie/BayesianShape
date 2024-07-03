#!/bin/bash

# Moves a slurm log to <file>.log, where <file>.h5 is the HDF5 file associated 
# with the job.

if [ $# -eq 0 ]; then
  echo "A Slurm log must be provided"
else
  fl=$1
  
  # Validate the input file.
  if [[ ! -e $fl ]]; then 
    echo "error: input file does not exist: $fl" 
    exit 1
  fi

  # Find the hdf5 file.
  h5fl=$( grep '^Wrote.*h5' $fl | sed 's|Wrote:* *\(.*h5\).*$|\1|' )

  # Repace the .h5 extension with .log.
  outfl=$( echo $h5fl | sed 's/h5$/log/' )
  
  # Don't overwrite output file.
  if [[ -e $outfl ]]; then 
    echo "error: target file exists: $outfl" 
    exit 1
  fi
  
  # `mv` the user tell them that we did so.
  cmd="mv $fl $outfl"
  echo $cmd
  $cmd
fi
