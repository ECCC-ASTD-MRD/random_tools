# random_tools

Set of random number generators

# Instructions in a nutshell

# At CMC only

```
# clone random_tools repository:
git clone --branch dev --recurse-submodules git@gitlab.science.gc.ca:RPN-SI/random_tools.git
cd random_tools
# use the appropriate setup file found in ci-env directory, for example:
. ./ci-env/latest/ubuntu-18.04-amd-64/gnu-9.3.0.sh
mkdir -p build
cd build
cmake .. -DBUILD_TESTING=true -DCMAKE_INSTALL_PREFIX=${your_choice}
make -j
# to launch tests
make check
# to install
make install
# to prepare a ssm package, use make package
make package
```

# Outside CMC (external users)

```
# clone random_tools repository:
# ci-env submodules is not available outside CMC and has to be excluded
git clone --branch dev --recurse-submodules=":(exclude)ci-env" https://github.com/ECCC-ASTD-MRD/random_tools
cd random_tools

# if you have already made a clone of random_tools without the
# exclusion of ci-env listed above, then use the following command:
git -c submodule."ci-env".update=none submodule update --init --recursive

# load compiler and cmake version 3.16 minimum
# build and compile random_tools
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
ulimit -s unlimited
```
