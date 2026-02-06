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

# Configure the build
cmake ..

# Build the project
cmake --build . -j4
```

The compiled executable will be in `build/pheniqs`.

### Installation

```bash
# From the build directory
cmake --install . --prefix /usr/local
```

Or manually:
```bash
make install
```

## Configuration Options

You can pass options to cmake during configuration:

```bash
# Specify installation prefix
cmake -DCMAKE_INSTALL_PREFIX=/opt/pheniqs ..

# Build with static libraries
cmake -DBUILD_STATIC=ON ..

# Specify a different compiler
cmake -DCMAKE_CXX_COMPILER=/usr/bin/g++-11 ..

# Release build with optimizations
cmake -DCMAKE_BUILD_TYPE=Release ..
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
