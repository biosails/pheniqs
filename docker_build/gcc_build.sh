#!/usr/bin/env bash

echo "Building with gcc"

set -x -e

rm -rf install.sh || echo "No install file found"
#export PREFIX=/tmp/pheniqs
export PWD=$(pwd)

cat <<EOF >>install.sh
#!/usr/bin/env bash

set -x -e

echo "Installing system packages"
apt-get update -y
apt-get install -y rsync

echo "Installing miniconda"
cd /tmp
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p /tmp/miniconda3

echo "Beginning pheniqs install"
cd $PWD

LD_LIBRARY_PATH="\${HOME}/.pheniqs/travis/install/lib"
export PATH=/tmp/miniconda3/bin/:\$PATH

./tool/ppkg.py -v info build test/build.json
make all PREFIX="\${HOME}/.pheniqs/travis/install"
./pheniqs --version
./pheniqs --help
./pheniqs demux --help
./test/BDGGG/run.sh

EOF

chmod 777 install.sh


docker run \
 	--rm \
	-it -v $(pwd):$(pwd) gcc:4.9 \
	$(pwd)/install.sh
