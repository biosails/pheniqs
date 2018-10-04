#!/usr/bin/env bash

# set -x -e

export WORKSPACE="/tmp/conda_recipe"

if [ "$TRAVIS_OS_NAME" == "linux" ]; then
	# WORKSPACE='/bioconda'
	# File mounts with docker don't seem to work the same as they do everywhere else
	# So we are just copying the recipe to tmp and working from there
	mkdir -p /tmp/conda_recipe
	cp -rf /bioconda/* /tmp/conda_recipe
	cd /tmp/conda_recipe
	export WORKSPACE="/tmp/conda_recipe"

	mkdir -p /opt/conda/conda-bld/{noarch,linux-64,osx-64}
	/opt/conda/bin/conda index $WORKSPACE/miniconda/conda-bld

else
	export WORKSPACE="/tmp/conda_recipe"
fi


UPLOAD_ARGS=""
if [ "$TRAVIS_BRANCH" == "master" ] && [ "$TRAVIS_BUILD_STAGE_NAME" == "Deploy" ]; then
	echo "Uploading package after successful build"
	UPLOAD_ARGS="--token $ANACONDA_API_TOKEN"
	conda config --set anaconda_upload yes
fi

DATE=$(date +"%Y%m%d")

cd $WORKSPACE

## Create a latest release
sed -i.bak 's/THIS_VERSION/latest/' latest/meta.yaml
rm latest/meta.yaml.bak
conda config --add channels conda-forge
conda config --add channels bioconda
conda build  $UPLOAD_ARGS $WORKSPACE/latest
conda build purge
cp meta.yaml latest/

## Create a version from this datetime
sed -i.bak "s/THIS_VERSION/$DATE/" latest/meta.yaml
rm latest/meta.yaml.bak
conda config --add channels conda-forge
conda config --add channels bioconda
conda build $UPLOAD_ARGS $WORKSPACE/latest

## Purge all the builds
conda build purge
