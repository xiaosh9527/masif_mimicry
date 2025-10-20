#!/bin/bash

# Install these dependencies in the root directory of the repository
# before running the postprocessing and truncation pipeline

#--------------------------------
# Installation of EvoEF2
#--------------------------------
# Clone from github repo
echo "Cloning EvoEF2 from github repo"
git clone https://github.com/tommyhuangthu/EvoEF2.git

# build EvoEF2.cpp to executable
echo "Building EvoEF2.cpp to executable"
cd EvoEF2
chmod +x ./build.sh
./build.sh

# Test if EvoEF2 executable is built successfully
./EvoEF2 --help

echo "Done"

#--------------------------------
# Installation of Stride
#--------------------------------
# wget the stride executable
echo "Downloading Stride executable"
cd ..
mkdir -p stride
cd stride
wget https://webclu.bio.wzw.tum.de/stride/stride.tar.gz

# untar the stride executable
echo "Untaring Stride executable"
tar -xzf stride.tar.gz

# compile the stride executable
make

# test if stride executable is downloaded successfully
./stride --help

echo "Done"