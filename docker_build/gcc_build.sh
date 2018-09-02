#!/usr/bin/env bash

echo "Building with gcc"

set -x -e

echo <<EOF >> install.sh
#!/usr/bin/env bash

set -x -e
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b

##Begin pheniqs install
export PREFIX=/tmp/pheniqs
export LD_LIBRARY_PATH="${PREFIX}/lib"
make all PREFIX=${PREFIX}
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
	-it -v $(pwd)/$(pwd) gcc:4.9 \
	$(pwd)/install.sh
