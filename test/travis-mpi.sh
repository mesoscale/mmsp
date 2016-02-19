#!/bin/sh
# Install the specified MPI environment in Travis-CI
# Copied from https://github.com/elemental/Elemental
# Questions/comments to trevor.keller@gmail.com (Trevor Keller)
set -e
case $1 in
  mpich2) set -x;
    sudo apt-get install -q mpich2 libmpich2-3 libmpich2-dev;;
  openmpi) set -x;
    sudo apt-get install openmpi-bin openmpi-dev;;
  *)
    echo "Unknown MPI implementation:" $1; exit 1;;
esac
