#!/bin/bash

FOLDER_NAME=${1:-"NWQ-Sim"}
echo "Install NWQ-Sim in folder: $FOLDER_NAME"

git clone https://github.com/pnnl/NWQ-Sim.git $FOLDER_NAME
cd ~/$FOLDER_NAME
git submodule update --init --recursive vqe/nlopt
git submodule update --init --recursive vqe/pybind11

source ~/"$FOLDER_NAME"/environment/setup_perlmutter.sh
module load python

cd ~/"$FOLDER_NAME"

mkdir ~/$FOLDER_NAME/build
cd ~/$FOLDER_NAME/build

cmake .. -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DVQE_ENABLE_MPI=ON -DCMAKE_CUDA_HOST_COMPILER=CC -DCMAKE_BUILD_TYPE=Release && make -j16

echo "NWQSim VQE Mem completed!"
