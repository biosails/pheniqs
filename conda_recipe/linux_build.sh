#!/usr/bin/env bash

echo "Building conda package on linux"

set -x -e

cp -rf conda_recipe /tmp/
cd /tmp

if [ "$TRAVIS_BRANCH" == "master" ]; then
	docker run \
		--user 'root' \
		-e ANACONDA_API_TOKEN=$ANACONDA_API_TOKEN \
		-e TRAVIS_OS_NAME=$TRAVIS_OS_NAME \
		-e TRAVIS_BRANCH=$TRAVIS_BRANCH \
		-it -v /tmp/conda_recipe:/bioconda condaforge/linux-anvil \
		/bioconda/conda_build.sh
else
	echo "This is not the master branch, no build"
	exit 0
fi
