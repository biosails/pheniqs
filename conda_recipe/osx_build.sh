#!/usr/bin/env bash

echo "Building conda package on osx"

set -x -e

tag="MacOSX"
WORKSPACE="/tmp"
MINICONDA_VER="latest"
THIS_DIR=$(pwd)

# step 1: download and install miniconda
# https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
curl -L -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-$MINICONDA_VER-$tag-x86_64.sh
bash miniconda.sh -b -p $WORKSPACE/miniconda

# step 2: setup channels
$WORKSPACE/miniconda/bin/conda config --system --add channels defaults
$WORKSPACE/miniconda/bin/conda config --system --add channels conda-forge
$WORKSPACE/miniconda/bin/conda config --system --add channels bioconda
$WORKSPACE/miniconda/bin/conda config --system --add channels nyuad-cgsb

# step 3: install conda build
$WORKSPACE/miniconda/bin/conda install -y conda-build anaconda-client conda-verify


# step 4: configure local channel
$WORKSPACE/miniconda/bin/conda index $WORKSPACE/miniconda/conda-bld
# $WORKSPACE/miniconda/bin/conda index $WORKSPACE/miniconda/conda-bld/linux-64 $WORKSPACE/miniconda/conda-bld/osx-64 $WORKSPACE/miniconda/conda-bld/noarch
$WORKSPACE/miniconda/bin/conda config --system --add channels file://$WORKSPACE/miniconda/conda-bld

# step 5: cleanup
$WORKSPACE/miniconda/bin/conda clean -y --all
rm miniconda.sh

export PATH=$WORKSPACE/miniconda/bin:$PATH

cp -rf conda_recipe /tmp
cd /tmp/conda_recipe
/tmp/conda_recipe/conda_build.sh

exit 0
