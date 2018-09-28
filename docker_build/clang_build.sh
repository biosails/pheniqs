#!/usr/bin/env bash

echo "Building clang 3.9 from ubuntu"

set -x -e

rm -rf install.sh || echo "No install file found"
export PREFIX=/tmp/pheniqs
export PWD=$(pwd)

cat <<EOF >>install.sh
cd /tmp
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b

apt-get update -y
apt-get install -y gnupg gnupg2 gnupg1 wget
wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key |  apt-key add -
apt-add-repository "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-6.0 main"
apt-get install -y clang-3.9 libclang-3.9-dev

##Begin pheniqs install
export LD_LIBRARY_PATH="${PREFIX}/lib"
make all PREFIX=${PREFIX}
make install PREFIX=${PREFIX}

./test/BDGGG/run.sh

EOF

chmod 777 install.sh

docker run \
	-it -v $(pwd):$(pwd) ubuntu \
	$(pwd)/install.sh
