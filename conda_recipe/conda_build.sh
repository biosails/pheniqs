#!/usr/bin/env bash

set -x -e

env | grep TRAVIS

if [ "$TRAVIS_BRANCH" == "master" ]; then
	DATE=$(date +"%Y%m%d%H%M")

	echo $ANACONDA_API_TOKEN

	cd /bioconda

	## Create a latest release
	sed -i 's/THIS_VERSION/latest/' latest/meta.yaml
	conda config --add channels bioconda
	conda config --set anaconda_upload yes
	conda build --token $ANACONDA_API_TOKEN /bioconda/latest 
	cp meta.yaml latest/

	## Create a version from this datetime
	sed -i "s/THIS_VERSION/$DATE/" latest/meta.yaml
	conda config --add channels bioconda
	conda config --set anaconda_upload yes
	conda build --token $ANACONDA_API_TOKEN /bioconda/latest 
fi

