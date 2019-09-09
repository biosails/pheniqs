#!/usr/bin/env bash

echo "Building conda package on linux"

set -x -e

if [ "$TRAVIS_BRANCH" == "master" ]; then
	echo "Building and uploading pheniqs for linux conda package"
	docker run \
		--user 'root' \
		-e ANACONDA_API_TOKEN=$ANACONDA_API_TOKEN \
		-e TRAVIS_OS_NAME=$TRAVIS_OS_NAME \
		-e TRAVIS_BRANCH=$TRAVIS_BRANCH \
		-e TRAVIS_BUILD_STAGE_NAME=$TRAVIS_BUILD_STAGE_NAME \
		-it -v $(pwd)/conda_recipe:/bioconda condaforge/linux-anvil \
		/bioconda/conda_build.sh
else
	echo "Building pheniqs for linux conda package"
	docker run \
		--user 'root' \
		-e TRAVIS_OS_NAME=$TRAVIS_OS_NAME \
		-e TRAVIS_BRANCH=$TRAVIS_BRANCH \
		-e TRAVIS_BUILD_STAGE_NAME=$TRAVIS_BUILD_STAGE_NAME \
		-it -v $(pwd)/conda_recipe:/bioconda condaforge/linux-anvil \
		/bioconda/conda_build.sh
fi
