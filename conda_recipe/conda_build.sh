#!/usr/bin/env bash

set -x -e

if [ "$TRAVIS_OS_NAME" == "linux" ]; then
	WORKSPACE='/bioconda'
else
	WORKSPACE='/tmp/conda_recipe'
fi

if [ "$TRAVIS_BRANCH" == "master" ]; then
	ls -lahR

	DATE=$(date +"%Y%m%d%H%M")

	cd $WORKSPACE

	## Create a latest release
	sed -i.bak 's/THIS_VERSION/latest/' latest/meta.yaml
	rm latest/meta.yaml.bak
	conda config --add channels bioconda
	conda config --set anaconda_upload yes
	conda build --token $ANACONDA_API_TOKEN $WORKSPACE/latest
        conda build purge
	cp meta.yaml latest/

	## Create a version from this datetime
	sed -i.bak "s/THIS_VERSION/$DATE/" latest/meta.yaml
	rm latest/meta.yaml.bak
	conda config --add channels bioconda
	conda config --set anaconda_upload yes
        conda build purge
	conda build --token $ANACONDA_API_TOKEN $WORKSPACE/latest
fi
