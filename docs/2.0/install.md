layout: page
title: "Install"
permalink: /install/

<section id="navigation">
    <ul>
        <li><a href="/pheniqs/2.0/background.html">Background</a></li>
        <li><a href="/pheniqs/2.0/cli.html">Cli</a></li>
        <li><a href="/pheniqs/2.0/configuration.html">Configuration</a></li>
        <li><a href="/pheniqs/2.0/glossary.html">Glossary</a></li>
        <li><a href="/pheniqs/2.0/">Home</a></li>
        <li><a href="/pheniqs/2.0/install.html">Install</a></li>
        <li><a href="/pheniqs/2.0/manual.html">Manual</a></li>
        <li><a href="/pheniqs/2.0/overview.html">Overview</a></li>
        <li><a href="/pheniqs/2.0/quickstart.html">Quickstart</a></li>
        <li><a href="/pheniqs/2.0/tokenization.html">Tokenization</a></li>
        <li><a href="/pheniqs/2.0/vignettes.html">Vignettes</a></li>
        <li><a href="/pheniqs/2.0/workflow.html">Workflow</a></li>
        <li><a class="github"   href="http://github.com/biosails/pheniqs">View on GitHub</a></li>
    </ul>
    <div class="clear" />
</section>


# Installation
{:.page-title}

* placeholder
{:toc}

Pheniqs is distributed as precompiled binaries, which may be installed with package managers. It may also be built from source code. Both methods are described below.

---

# Prebuilt Binaries

## Installing with homebrew on MacOS
**coming soon**

## Installing with conda

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

## *Building with `pheniqs-build-api.py`*
Pheniqs comes bundled with a Python3 helper tool called `pheniqs-build-api.py`. To build an entire virtual root of all the dependencies and compile a [statically linked](https://en.wikipedia.org/wiki/Static_library), portable, binary snapshot of the latest code against them simply execute `./tool/pheniqs-build-api.py build build/trunk_static.json` in the code root folder. The `build` folder contains several other configurations for official releases. Building with `pheniqs-build-api.py` does not require elevated permissions and is ideal for building an executable on cluster environments.

>```shell
% ./tool/pheniqs-build-api.py build build/trunk_static.json
INFO:Package:unpacking zlib 1.2.11
INFO:Package:configuring make environment zlib 1.2.11
INFO:Package:building with make zlib 1.2.11
INFO:Package:installing with make zlib 1.2.11
INFO:Package:unpacking bz2 1.0.6
INFO:Package:building with make bz2 1.0.6
INFO:Package:installing with make bz2 1.0.6
INFO:Package:unpacking xz 5.2.4
INFO:Package:configuring make environment xz 5.2.4
INFO:Package:building with make xz 5.2.4
INFO:Package:installing with make xz 5.2.4
INFO:Package:unpacking libdeflate 1.0
INFO:Package:building with make libdeflate 1.0
INFO:Package:unpacking htslib 1.9
INFO:Package:configuring make environment htslib 1.9
INFO:Package:building with make htslib 1.9
INFO:Package:installing with make htslib 1.9
INFO:Package:unpacking rapidjson 1.1.0
INFO:Package:downloaded archive saved pheniqs 2.0-trunk None
INFO:Package:unpacking pheniqs 2.0-trunk
INFO:Package:building with make pheniqs 2.0-trunk
INFO:Package:installing with make pheniqs 2.0-trunk
```

When `pheniqs-build-api.py` is done you may inspect your binary, statically linked builds made with `pheniqs-build-api.py` will also report the versions of all built in libraries.

>```shell
% ./bin/trunk_static/install/bin/pheniqs --version
pheniqs version 2.0.4
zlib 1.2.11
bzlib 1.0.6
xzlib 5.2.4
libdeflate 1.0
rapidjson 1.1.0
htslib 1.9
```

You can check that your binary indeed does not link against any of the dependencies dynamically with `otool` on MacOs:

>```shell
% otool -L ./bin/trunk_static/install/bin/pheniqs
./bin/trunk_static/install/bin/pheniqs:
	/usr/lib/libc++.1.dylib (compatibility version 1.0.0, current version 400.9.4)
	/usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1252.200.5)
```

Or `ldd` on Ubuntu:

>```shell
% ldd pheniqs
	linux-vdso.so.1 =>  (0x00007ffff3300000)
	libstdc++.so.6 => /usr/lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007f6910e2d000)
	libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007f6910b24000)
	libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007f691090e000)
	libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007f69106f1000)
	libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007f6910327000)
	/lib64/ld-linux-x86-64.so.2 (0x00007f69111af000)
```

## *Building with the Makefile*

Pheniqs does not use automake and so does not have a configure stage. The provided Makefile will build pheniqs against existing dependencies, if they are already present. Simply execute `make && make install`. You can execute `make help` for some general instructions.

### *Dependecies on Ubuntu*
All Pheniqs build dependencies are available on [Ubuntu 19.04 Disco Dingo](http://releases.ubuntu.com/19.04) and can be installed with:

>```shell
apt-get install -y \
build-essential \
rapidjson-dev \
libhts-dev \
liblzma-dev \
libdeflate-dev \
libbz2-dev
```

If you want to build Pheniqs against a specific root you may provide a `PREFIX` parameter, but notice that you need to specify it on each make invocation, for instance `make PREFIX=/usr/local && make install PREFIX=/usr/local`. Pheniqs is regularly tested on several versions of both [Clang](https://clang.llvm.org) and [GCC](https://gcc.gnu.org), you can tell `make` which compiler to use by setting the `CXX` parameter. See [travis](https://travis-ci.org/biosails/pheniqs) for a comprehensive list and test results.
