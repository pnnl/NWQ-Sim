#!/bin/bash

FOLDER_NAME=${1:-"NWQ-Sim"}
echo "Install in folder: $FOLDER_NAME"

# git clone --recursive -b vqe_memory https://github.com/pnnl/NWQ-Sim.git $FOLDER_NAME

source ~/"$FOLDER_NAME"/environment/setup_perlmutter.sh
module load python

cd ~/"$FOLDER_NAME"

rm -rf build
mkdir build; cd build

cmake .. -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_CUDA_HOST_COMPILER=CC -DCMAKE_BUILD_TYPE=Release && make -j10

echo "NWQSim VQE Mem completed!"