Set of random number generators

# At CMC

## Build dependencies

- CMake 3.20+

## Environment

Load the right environment, depending on the architecture you need.  This
will load the specified compiler and its parameters, and set the
`EC_CMAKE_MODULE_PATH` variable for the `cmake_rpn` modules.

- Example for ppp6/sc6 and icelake specific architecture:

```
. r.load.dot mrd/rpn/code-tools/latest/env/rhel-8-icelake-64@inteloneapi-2025.1.0
```

- Example for generic architecture on ppp6/sc6

```
. r.load.dot mrd/rpn/code-tools/latest/env/rhel-8-amd64-64@inteloneapi-2025.1.0
```

- Example for GNU on any architecture:

```
. r.load.dot mrd/rpn/code-tools/latest/env/gnu
```

## Build and install

```
mkdir build
cd build
cmake .. -DBUILD_TESTING=true -DCMAKE_INSTALL_PREFIX=${your_choice}
make -j
# to launch tests
make check
```

# Outside CMC (external users)

## Build dependencies

- CMake 3.20+

`cmake_rpn` is included as a git submodule.  Please clone with the
`--recurse --remote-submodules` options, or run `git submodule update --init
--remote` in the git repo after having cloned.

## Build and install

```
mkdir -p build
cd build
cmake .. -DBUILD_TESTING=true -DCMAKE_INSTALL_PREFIX=${your_choice}
make -j
# to launch tests
make check
# to install
make install
```

# Documentation

See documentation on functions in doc/randomgeneric.html

# Troubleshooting

If you encounter problems with the tests, it may because stack size is
limited on your computer. You could have to change its limit, for example:
```
ulimit -s 4000000
```
