#!/usr/bin/env bash

# step 1: download and install miniconda

# https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
tag="MacOSX"
WORKSPACE=/tmp
MINICONDA_VER="latest"
THIS_DIR=$(pwd)

curl -L -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-$MINICONDA_VER-$tag-x86_64.sh
bash miniconda.sh -b -p $WORKSPACE/miniconda

# step 2: setup channels
$WORKSPACE/miniconda/bin/conda config --system --add channels defaults
$WORKSPACE/miniconda/bin/conda config --system --add channels conda-forge
$WORKSPACE/miniconda/bin/conda config --system --add channels bioconda
$WORKSPACE/miniconda/bin/conda config --system --add channels nyuad-cgsb 


# step 4: configure local channel
$WORKSPACE/miniconda/bin/conda index $WORKSPACE/miniconda/conda-bld/linux-64 $WORKSPACE/miniconda/conda-bld/osx-64 $WORKSPACE/miniconda/conda-bld/noarch
$WORKSPACE/miniconda/bin/conda config --system --add channels file://$WORKSPACE/miniconda/conda-bld

# step 5: cleanup
$WORKSPACE/miniconda/bin/conda clean -y --all
$WORKSPACE/miniconda/bin/conda -y conda-build 
rm miniconda.sh

cp -rf conda_recipe /tmp
/tmp/conda_recipe/conda_build.sh
