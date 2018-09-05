#!/usr/bin/env bash

echo "Building with gcc"

set -x -e

rm -rf install.sh || echo "No install file found"
export PREFIX=/tmp/pheniqs
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
export LD_LIBRARY_PATH="${PREFIX}/lib"
export PATH=/tmp/miniconda3/bin/:\$PATH

./tool/ppkg.py -v debug build test/build.json
make PREFIX=${PREFIX}
make install PREFIX=${PREFIX}

pheniqs demux --config test/BDGGG/BDGGG_interleave.json --validate --distance
pheniqs demux --config test/BDGGG/BDGGG_interleave.json --compile >> test/BDGGG/BDGGG.log 2>&1
pheniqs demux --config test/BDGGG/BDGGG_interleave.json >> test/BDGGG/BDGGG.log 2>&1
pheniqs demux --config test/BDGGG/BDGGG_annotated.json --validate --distance
pheniqs demux --config test/BDGGG/BDGGG_annotated.json --compile >> test/BDGGG/BDGGG.log 2>&1
pheniqs demux --config test/BDGGG/BDGGG_annotated.json >> test/BDGGG/BDGGG.log 2>&1
EOF

chmod 777 install.sh


docker run \
 	--rm \
	-it -v $(pwd):$(pwd) gcc:4.9 \
	$(pwd)/install.sh
