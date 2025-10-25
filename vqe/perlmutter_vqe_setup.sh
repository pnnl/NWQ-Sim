#!/bin/bash

FOLDER_NAME=${1:-"NWQ-Sim"}
echo "Install NWQ-Sim in folder: $FOLDER_NAME"

git clone https://github.com/pnnl/NWQ-Sim.git $FOLDER_NAME
cd ~/$FOLDER_NAME

# Initialize all submodules recursively
git submodule update --init --recursive

# Verify critical submodules are present
if [ ! -f "vqe/pybind11/CMakeLists.txt" ]; then
    echo "ERROR: pybind11 submodule not initialized. Re-attempting..."
    git submodule update --init vqe/pybind11
    if [ ! -f "vqe/pybind11/CMakeLists.txt" ]; then
        echo "FATAL: Failed to initialize pybind11 submodule"
        exit 1
    fi
fi

if [ ! -f "vqe/nlopt/CMakeLists.txt" ]; then
    echo "ERROR: nlopt submodule not initialized. Re-attempting..."
    git submodule update --init vqe/nlopt
    if [ ! -f "vqe/nlopt/CMakeLists.txt" ]; then
        echo "FATAL: Failed to initialize nlopt submodule"
        exit 1
    fi
fi

source ~/"$FOLDER_NAME"/environment/setup_perlmutter.sh
module load python

cd ~/$FOLDER_NAME

# Create build directory
mkdir -p ~/$FOLDER_NAME/build
cd ~/$FOLDER_NAME/build

# Run CMake with pybind11 submodule support
cmake .. \
    -DCMAKE_C_COMPILER=cc \
    -DCMAKE_CXX_COMPILER=CC \
    -DNWQSIM_ENABLE_VQE=ON \
    -DVQE_ENABLE_MPI=ON \
    -DCMAKE_CUDA_HOST_COMPILER=CC \
    -DCMAKE_BUILD_TYPE=Release

make -j16

if [ $? -eq 0 ]; then
    echo "NWQSim VQE CMake configuration completed successfully!"
else
    echo "ERROR: CMake configuration failed"
    exit 1
fi
