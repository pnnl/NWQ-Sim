# Shared by dmsim and svsim under single
# This is for NERSC Perlmutter 

###################### NVIDIA GPU Configuration #######################
#NVCC =/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/cuda/12.2/bin/nvcc 
NVCC = nvcc
NVCC_FLAGS = -O3 -arch=sm_80 -m64 -std=c++17 -rdc=true --compiler-options -fPIC
NVCC_LIBS = 
#######################################################################

####################### Host CPU Configuration ########################
CC = CC
CC_FLAGS = -O3 -m64 -std=c++17 -fPIC
CC_LIBS = 
#######################################################################

##################### Host Q#/QIR Configuration #######################
QIRCC = /global/common/software/nersc/pe/llvm/17/bin/clang++
QIRCC_FLAGS = -std=c++17 -m64 -O3 -I. -fPIC
QIR_BRIDGE_PUBLIC = /global/homes/a/angli/qsharp-runtime/src/Qir/Runtime/public/
QIR_BRIDGE_TEST = 
QIR_BRIDGE_BUILD = /global/homes/a/angli/qsharp-runtime/src/Qir/Runtime/build/
QIR_BRIDGE_FLAGS = -I. -I$(QIR_BRIDGE_PUBLIC) -L$(QIR_BRIDGE_BUILD)/lib/QIR -L$(QIR_BRIDGE_BUILD)/lib/QSharpCore -L$(QIR_BRIDGE_BUILD)/lib/QSharpFoundation -L$(QIR_BRIDGE_BUILD)/lib/Tracer -lMicrosoft.Quantum.Qir.Runtime -lMicrosoft.Quantum.Qir.QSharp.Core -lMicrosoft.Quantum.Qir.QSharp.Foundation -lMicrosoft.Quantum.Qir.Tracer
#######################################################################

