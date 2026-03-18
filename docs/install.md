---
layout: default
title: "Build and Install"
permalink: /install
id: install
---

* placeholder
{:toc}

Pheniqs is distributed as precompiled binaries, which may be installed with package managers. It may also be built from source code. Both methods are described below.

---

# Prebuilt Binaries

## Installing with homebrew on MacOS
**coming soon**

## Installing with conda

Stable realeases of Pheniqs are [available on bioconda](https://anaconda.org/bioconda/pheniqs). 

### *One time setup - Install Miniconda*
The easiest way to do this is to head on over to [**Anaconda**](https://conda.io/miniconda.html) and select the correct distribution for Python3.

### *One time setup - Configure your Conda channels*
Many groups contribute software to **Conda**. Each of these groups corresponds to a different channel. Bioconda is a well known channel for distributing bioinformatics software. It depends on **conda-forge**, another group for distributing more general software, including R and Python packages.

>```shell
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

### *Install Pheniqs*
Simply install pheniqs using the conda package manager.

>```shell
# Installs pheniqs in your 'global' conda install
conda install pheniqs
# Installs pheniqs in an isolated environment - recommended
conda create -n pheniqs pheniqs
```

If you want to live on the bleeding edge, you can also install pheniqs from our anaconda channel, nyuad-cgsb.

>```shell
conda install -c nyuad-cgsb pheniqs/latest
```

---

# Build from Source

## *Dependencies*
Pheniqs depends on [HTSlib](http://www.htslib.org), [RapidJSON](http://rapidjson.org) and [zlib](https://zlib.net). HTSLib further depends on [bzip2](http://www.bzip.org), [LZMA](https://tukaani.org/xz) and optionally [libdeflate](https://github.com/ebiggers/libdeflate) for improved gzip compressed FASTQ manipulation. Pheniqs requires [HTSLib version 1.8](https://github.com/samtools/htslib/releases/tag/1.8) or later and [RapidJSON version 1.1.0](https://github.com/Tencent/rapidjson/releases/tag/v1.1.0) or later. The versions packaged in most linux distributions are very outdated and cannot be used to build Pheniqs.

## *Building with CMake*

CMake is the recommended build system. It requires CMake 3.10 or later and can automatically download and build vendored copies of all dependencies if system libraries are not available.

### *Quick build against system libraries*

Install dependencies first:

On Ubuntu/Debian:

>```shell
apt-get install -y \
build-essential \
cmake \
rapidjson-dev \
libhts-dev \
liblzma-dev \
libdeflate-dev \
libbz2-dev \
libssl-dev
```

On macOS with Homebrew:

>```shell
brew install cmake zlib bzip2 rapidjson xz htslib libdeflate
```
(libcrypto is not needed on macOS — CommonCrypto is part of the system.)

Then build:

>```shell
mkdir build && cd build
cmake ..
cmake --build . -j$(nproc)
cmake --install . --prefix /usr/local
```

### *Fully vendored static build*

To produce a portable, statically linked binary with all dependencies built from source — no system libraries required:

>```shell
mkdir build && cd build
cmake -DBUILD_STATIC=ON \
      -DVENDOR_ZLIB=ON \
      -DVENDOR_BZIP2=ON \
      -DVENDOR_LZMA=ON \
      -DVENDOR_LIBDEFLATE=ON \
      -DVENDOR_HTSLIB=ON \
      -DVENDOR_RAPIDJSON=ON ..
cmake --build . -j$(nproc)
```

Dependencies are downloaded and built automatically into `build/external/`. This does not require elevated permissions and is ideal for cluster or cloud environments. Vendoring htslib, zlib, bzip2, or xz additionally requires `autoconf`, `automake`, and `libtool`.

For the full list of CMake options, see [CMAKE_BUILD_GUIDE.md](CMAKE_BUILD_GUIDE.md) in the repository root.

## *Building with the Makefile*

A legacy Makefile is also provided. It builds pheniqs against existing system dependencies and does not support vendoring or automatic dependency download. Execute `make && make install`. You can run `make help` for instructions.

### *Dependencies on Ubuntu*
All Pheniqs build dependencies are available on [Ubuntu 20.04 Focal Fossa](http://releases.ubuntu.com/20.04) and can be installed with:

>```shell
apt-get install -y \
build-essential \
rapidjson-dev \
libhts-dev \
liblzma-dev \
libdeflate-dev \
libbz2-dev \
libssl-dev
```

### *Dependencies on MacOS*
All Pheniqs build dependencies are available on [homebrew](https://brew.sh) and can be installed with:

>```shell
brew install \
zlib \
bzip2 \
rapidjson \
xz \
htslib \
libdeflate
```

If you want to build Pheniqs against a specific root you may provide a `PREFIX` parameter, but notice that you need to specify it on each make invocation, for instance `make PREFIX=/usr/local && make install PREFIX=/usr/local`. You can tell `make` which compiler to use by setting the `CXX` parameter.
