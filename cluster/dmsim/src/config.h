// ---------------------------------------------------------------------------
// NWQSim: Northwest Quantum Simulation Environment 
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/SV-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
// File: config.h
// Configuration file. Some of the options here are configured by compiler flags.
// ---------------------------------------------------------------------------
#ifndef CONFIG_H_
#define CONFIG_H_
/***********************************************
 * Backend configuration:
 *   * USE_CPU: use CPU backend
 *   * USE_NVGPU: use NVIDIA GPU backend
 *   * USE_AMDGPU: use AMD GPU backend
 *   * USE_AVX512: enable AVX512 vectorization
 * The default is USE_CPU.
 ***********************************************/
//#define USE_CPU
/*#define USE_NVGPU*/
//#define USE_AMDGPU
//#define USE_AVX512

/***********************************************
 * Parallel configuration:
 *    * USE_OMP: use single-node multi-devices
 *    * USE_MPI: use multi-nodes
 * The default is to use a single device.
 ***********************************************/
//#define USE_OMP
//#define USE_MPI

/***********************************************
 * Runtime configuration:
 ***********************************************/
/* Default number of CPUs/GPUs in OMP mode */
#define DEFAULT_OMP_PE 4
/* Default GPU device to be used in single device mode */
#define DEFAULT_SIN_GPU 0
/* Default number of repetitions */
#define DEFAULT_REPETITIONS 1
/* Print out per kernel execution trace */
#define PRINT_KERNEL_TRACE
/* Print out gate trace */
/*#define PRINT_GATE_TRACE */
/* Print out qubit allocation, release and reset trace */
/*#define PRINT_QUBIT_TRACE */
/* Print out circuit trace: reset, optimization, transform, etc.*/
/*#define PRINT_CIRCUIT_TRACE */
/* Print out simulation trace: reset, run, initialize*/
/*#define PRINT_SIM_TRACE */
/* Perform purity check per gate. Useful for debugging. 
 * It breaks when the purity is far from 1. */
/*#define PURITY_CHECK*/
/* Perform GPU error check */
#define GPU_ERROR_CHECK
/* Number of qubit slots. Depends on memory capacity */
#define N_QUBIT_SLOT 32
/* Number of threads per thread block for NVGPU */
#define THREADS_CTA_NVGPU 256 
/* Number of threads per thread block for AMDGPU */
#define THREADS_CTA_AMDGPU 256 
/* Error bar for purity check and other error check */
#define ERROR_BAR (1e-3)
/* Disable noise injection for functional validation */
#define DISABLE_NOISE
/* Disable gate fusion */
/*#define DISABLE_GATE_FUSION*/
/* Enable the double-precision Tensor Core support in A100 GPUs */
/*#define ENABLE_TENSOR_CORES*/
/* NISQ device to be simulated */
/*#define DEVICE_CONFIG_NAME "ibmq_guadalupeConfig"*/
#define DEVICE_CONFIG_NAME "dummy_ibmq11"
/* NISQ deivce config file path */
#define DEVICE_CONFIG_PATH "/global/homes/a/angli/sc22/NWQ-Sim/device/"

/***********************************************
 * Constant configuration:
 ***********************************************/
/* Constant value of PI */
#define PI 3.14159265358979323846
/* Constant value of sqrt(2) */
#define S2I 0.70710678118654752440 
/* Constant value of 0.5 */
#define HALF 0.5

/***********************************************
 * Data type configuration:
 ***********************************************/
namespace NWQSim {
/* Basic data type for indices */
using IdxType = long long int;
/* Basic data type for value */
using ValType = double;
}; //namespace DMSim
#endif
