#!/usr/bin/env bash

docker run \
	-e ANACONDA_API_TOKEN=$ANACONDA_API_TOKEN \
	-e TRAVIS_OS_NAME=$TRAVIS_OS_NAME \
	-e TRAVIS_BRANCH=$TRAVIS_BRANCH \
	-it -v `pwd`:/bioconda condaforge/linux-anvil /bioconda/conda_build.sh
