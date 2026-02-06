# Building Pheniqs with CMake

This project has been converted from a Makefile-based build system to CMake. CMake is a more portable and maintainable build system that works across different operating systems and compilers.

## Quick Start

### Prerequisites

Before building Pheniqs, ensure you have the following libraries installed:

- **zlib**: https://zlib.net
- **bzip2**: http://www.bzip.org
- **xz/lzma**: https://tukaani.org/xz
- **htslib**: http://www.htslib.org
- **libdeflate** (optional): https://github.com/ebiggers/libdeflate

On Ubuntu/Debian:
```bash
sudo apt-get install zlib1g-dev libbz2-dev liblzma-dev libhts-dev libdeflate-dev
```

On macOS with Homebrew:
```bash
brew install zlib bzip2 xz htslib libdeflate
```

Also ensure you have CMake 3.10+ installed:
```bash
cmake --version
```

### Basic Build

```bash
# Create a build directory
mkdir build
cd build

# Configure the build (Release build is recommended for production)
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build the project
cmake --build . -j4
```

The compiled executable will be in `build/pheniqs`.

### Quick Build with Default Settings

If system libraries are available, the simplest command is:

```bash
mkdir build && cd build
cmake .. && cmake --build . -j$(nproc)
```

This will automatically use system-installed libraries for zlib, bzip2, xz, htslib, and libdeflate if available.

### Installation

```bash
# From the build directory
cmake --install . --prefix /usr/local
```

Or manually:
```bash
make install
```

### Uninstall

An `uninstall` target is provided to remove files that were installed by `cmake --install`.

The build configures a small uninstall script at configure time which removes the files listed in the `install_manifest.txt` produced during install. The uninstaller looks for the manifest in the build directory first, and falls back to `${CMAKE_INSTALL_PREFIX}/install_manifest.txt`.

To run the uninstaller from your build directory:

```bash
cmake --build . --target uninstall
```

If you installed to a system prefix (for example `/usr/local`) as root, run the uninstall with `sudo`:

```bash
sudo cmake --build . --target uninstall
```

Files added to the source tree to support this:

- `cmake/cmake_uninstall.cmake.in` — uninstall script template configured into the build tree
- `CMakeLists.txt` — now configures the uninstall script and exposes the `uninstall` custom target


## Configuration Options

You can pass options to cmake during configuration:

```bash
# Specify installation prefix
cmake -DCMAKE_INSTALL_PREFIX=/opt/pheniqs ..

# Build with static libraries
cmake -DBUILD_STATIC=ON ..

# Release build with optimizations (recommended for production)
cmake -DCMAKE_BUILD_TYPE=Release ..

# Specify a different compiler
cmake -DCMAKE_CXX_COMPILER=/usr/bin/g++-11 ..
```

### Vendoring Options

Pheniqs supports optional vendored (locally-built) versions of dependencies. All are disabled by default, so the build will use system-installed libraries when available:

```bash
# Force building zlib locally
cmake -DVENDOR_ZLIB=ON ..

# Force building bzip2 locally
cmake -DVENDOR_BZIP2=ON ..

# Force building xz/liblzma locally
cmake -DVENDOR_LZMA=ON ..

# Force building htslib locally
cmake -DVENDOR_HTSLIB=ON ..

# Force building libdeflate locally
cmake -DVENDOR_LIBDEFLATE=ON ..

# Force building RapidJSON locally (header-only)
cmake -DVENDOR_RAPIDJSON=ON ..
```

**Note**: Building vendored dependencies requires `autoconf`, `automake`, and `libtool`:
```bash
sudo apt-get install autoconf automake libtool  # Ubuntu/Debian
brew install autoconf automake libtool          # macOS
```

## Common Tasks

### View all available targets
```bash
cmake --build . --target help
```

Or from the build directory:
```bash
make help
```

### Running tests
```bash
cmake --build . --target test
```

### Cleaning build artifacts
```bash
rm -rf build
```

Or to clean specific test results:
```bash
cmake --build . --target clean.test.pheniqs.BDGGG
```

### Install without using user directories
```bash
cmake --install . --prefix ./install
```

## Environment Detection

CMake will automatically:

1. Detect your system (Linux, macOS, etc.)
2. Find installed libraries in standard locations
3. Determine the appropriate compiler flags
4. Generate `version.h` from git describe (if available)
5. Generate `configuration.h` from `configuration.json`
6. Generate zsh completion script `_pheniqs`

### Library Detection Priority

CMake searches for libraries in the following order:
1. System default paths (`/usr/lib`, `/usr/local/lib`, etc.)
2. Paths specified by `CMAKE_PREFIX_PATH`
3. If not found, and vendoring is enabled, build locally

## Build Recommendations

### For Production/Release

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j$(nproc)
cmake --install . --prefix /usr/local
```

This uses system libraries when available and applies optimization flags.

### For Development

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build . -j$(nproc)
```

This includes debugging symbols but may be slower.

### For Distribution/All Features

If you want to include all dependencies in the build:

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DVENDOR_ZLIB=ON \
      -DVENDOR_BZIP2=ON \
      -DVENDOR_LZMA=ON \
      -DVENDOR_HTSLIB=ON \
      -DBUILD_STATIC=ON ..
cmake --build . -j$(nproc)
```

This produces a self-contained, portable binary.

## Differences from the Makefile

| Feature | Makefile | CMake |
|---------|----------|-------|
| Build directory | Current directory | Separate `build/` directory |
| Configuration | `make config` | `cmake --system-information` |
| Help | `make help` | `cmake --build . --target pheniqs-help` |
| Version detection | Via git describe | Automatic via CMake |
| Installation | `make install` | `cmake --install .` |
| Parallel builds | `make -j4` | `cmake --build . -j4` |

## Troubleshooting

### CMake not found
Install CMake 3.10 or later:
```bash
sudo apt-get install cmake  # Ubuntu/Debian
brew install cmake          # macOS
```

### Build downloads dependencies despite system libraries being available

This typically means the CMake cache is stale. Clean and reconfigure:
```bash
rm -rf build
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j4
```

If libraries are still not found, verify they're installed:
```bash
# Check for zlib
pkg-config --cflags --libs zlib
# Check for bzip2
dpkg -l | grep bzip2-dev  # or rpm -qa | grep bzip2-devel
# Check for xz
dpkg -l | grep liblzma-dev  # or rpm -qa | grep xz-devel
```

### Vendored build fails with "configure: not found"

If building vendored dependencies, ensure autotools are installed:
```bash
sudo apt-get install autoconf automake libtool  # Ubuntu/Debian
brew install autoconf automake libtool          # macOS
```

### Library not found
If CMake can't find a required library, specify its location:
```bash
cmake -DCMAKE_PREFIX_PATH=/opt/htslib ..
```

### Compiler not found
If you want to use a specific compiler:
```bash
cmake -DCMAKE_CXX_COMPILER=/path/to/compiler ..
```

### Python script error when generating configuration.h
Ensure Python 3 is installed and accessible:
```bash
which python3
```

### RapidJSON compilation errors

If you encounter RapidJSON-related compilation errors like `memcpy` warnings or API mismatches:
- The recommended solution is to use the system or default RapidJSON (header-only)
- Do not use `-DVENDOR_RAPIDJSON=ON` unless necessary
- If vendoring is required, ensure compiler warnings are compatible with v1.1.0

## For Developers

### Modifying the build system

The main CMake configuration is in `CMakeLists.txt`. Key sections:

- **Project setup**: Lines 1-70 (versions, compiler setup)
- **Dependencies**: Lines 72-140 (finding libraries)
- **Source files**: Lines 165-195 (list of sources)
- **Custom generation**: Lines 197-255 (version.h, configuration.h, _pheniqs)
- **Install rules**: Lines 285-290

To add new source files, simply add them to the `PHENIQS_SOURCES` list in CMakeLists.txt.

### Adding new make targets to CMake equivalents

To add new custom targets, use:
```cmake
add_custom_target(my-target
    COMMAND some-command arg1 arg2
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS pheniqs
)
```

## Additional Resources

- CMake Documentation: https://cmake.org/documentation/
- CMake Tutorial: https://cmake.org/cmake-tutorial/
- Finding Packages with CMake: https://cmake.org/cmake/help/latest/command/find_package.html
