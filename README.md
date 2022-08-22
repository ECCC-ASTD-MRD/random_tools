# random_tools

Set of random number generators

# Instructions in a nutshell

# At CMC only

```
git clone --branch dev --recurse-submodules git@gitlab.science.gc.ca:RPN-SI/random_tools.git
cd random_tools
# use the appropriate setup file in ci-env
mkdir -p build
cd build
cmake -DBUILD_TESTING=true -DCMAKE_INSTALL_PREFIX=${your_choice} ..
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
git clone --branch dev https://github.com/ECCC-ASTD-MRD/random_tools
cd random_tools
git -c submodule."ci-env".update=none submodule update --init --recursive
mkdir -p build
cd build
cmake -DBUILD_TESTING=true -DCMAKE_INSTALL_PREFIX=${your_choice} ..
make -j
# to launch tests
make check
# to install
make install
```

# Troubleshooting

If you encounter problems with the tests, it may because stack size is
limited on your computer. You could have to change its limit:
```
ulimit -s unlimited
```
