// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/DM-Sim
// PNNL-IPID: 31919-E, ECCN: EAR99, IR: PNNL-SA-143160
// GitHub repo: http://www.github.com/pnnl/DM-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
// File: util.h
// Define utility functions.
// ---------------------------------------------------------------------------
#ifndef UTIL_H_
#define UTIL_H_

#include <stdio.h>
#include <sys/time.h>
#include <assert.h>
#include <iostream>
#include "config.h"
namespace NWQSim {
#ifdef USE_NVGPU
//==================================== NVGPU =======================================
/***********************************************
 * Error Checking:
 ***********************************************/
// Error checking for CUDA API call
#define cudaSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
inline void __cudaSafeCall(cudaError err, const char *file, const int line)
{
#ifdef GPU_ERROR_CHECK
    if (cudaSuccess != err)
    {
        fprintf(stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString(err));
        exit(-1);
    }
#endif
    return;
}
// Error checking for CUDA kernel call
#define cudaCheckError()  __cudaCheckError(__FILE__, __LINE__)
inline void __cudaCheckError(const char *file, const int line)
{
#ifdef GPU_ERROR_CHECK
    cudaError err = cudaGetLastError();
    if (cudaSuccess != err)
    {
        fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                file, line, cudaGetErrorString(err));
        exit(-1);
    }
    // Expensive checking
    err = cudaDeviceSynchronize();
    if( cudaSuccess != err )
    {
        fprintf(stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
                file, line, cudaGetErrorString(err));
        exit(-1);
    }
#endif
    return;
}
// Error checking for NVSHMEM call
#define NVSHMEM_CHECK(stmt)                                           \
    do {                                                              \
        int result = (stmt);                                          \
        if (NVSHMEMX_SUCCESS != result) {                             \
            fprintf(stderr, "[%s:%d] nvshmem failed with error %d \n",\
                    __FILE__, __LINE__, result);                      \
            exit(-1);                                                 \
        }                                                             \
    } while (0)
/***********************************************
 * Memory allocation and free
 ***********************************************/
//CPU host allocation with pinned memory
#define SAFE_ALOC_HOST(X,Y) cudaSafeCall(cudaMallocHost((void**)&(X),(Y)));
//GPU device allocation
#define SAFE_ALOC_GPU(X,Y) cudaSafeCall(cudaMalloc((void**)&(X),(Y)));
//CPU host free for pinned memory
#define SAFE_FREE_HOST(X) if ((X) != NULL){     \
               cudaSafeCall(cudaFreeHost((X))); \
               (X) = NULL;}
//GPU device free
#define SAFE_FREE_GPU(X) if ((X) != NULL){  \
               cudaSafeCall(cudaFree((X))); \
               (X) = NULL;}
/***********************************************
 * GPU Timer based on CUDA event
 ***********************************************/
typedef struct GPU_Timer
{
    GPU_Timer()
    {
        cudaSafeCall( cudaEventCreate(&this->start) );
        cudaSafeCall( cudaEventCreate(&this->stop) );
    }
    ~GPU_Timer()
    {
        cudaSafeCall( cudaEventDestroy(this->start));
        cudaSafeCall( cudaEventDestroy(this->stop));
    }
    void start_timer() { cudaSafeCall( cudaEventRecord(this->start) ); }
    void stop_timer() { cudaSafeCall( cudaEventRecord(this->stop) ); }
    double measure()
    {
        cudaSafeCall( cudaEventSynchronize(this->stop) );
        float millisconds = 0;
        cudaSafeCall(cudaEventElapsedTime(&millisconds, this->start, this->stop) ); 
        return (double)millisconds;
    }
    cudaEvent_t start;
    cudaEvent_t stop;
} gpu_timer;

#elif defined USE_AMDGPU
//==================================== AMDGPU =======================================
/***********************************************
 * Error Checking:
 ***********************************************/
// Error checking for HIP API call
#define hipSafeCall(err) __hipSafeCall(err, __FILE__, __LINE__)
inline void __hipSafeCall(hipError_t err, const char *file, const int line)
{
#ifdef GPU_ERROR_CHECK
    if ( hipSuccess != err )
    {
        fprintf(stderr, "hipSafeCall() failed at %s:%i : %s\n",
                 file, line, hipGetErrorString(err));
        exit(-1);
    }
#endif
    return;
}
// Error checking for HIP kernel call
#define hipCheckError() __hipCheckError(__FILE__, __LINE__)
inline void __hipCheckError(const char *file, const int line)
{
#ifdef GPU_ERROR_CHECK
    hipError_t err = hipGetLastError();
    if (hipSuccess != err)
    {
        fprintf(stderr, "hipCheckError() failed at %s:%i : %s\n",
                file, line, hipGetErrorString(err));
        exit( -1 );
    }
    // Expensive checking
    err = hipDeviceSynchronize();
    if( hipSuccess != err )
    {
        fprintf( stderr, "hipCheckError() with sync failed at %s:%i : %s\n",
                 file, line, hipGetErrorString( err ) );
        exit( -1 );
    }
#endif
    return;
}
/***********************************************
 * Memory allocation and free
 ***********************************************/
//CPU host allocation
#define SAFE_ALOC_HOST(X,Y) hipSafeCall(hipHostMalloc((void**)&(X),(Y)));
//GPU device allocation
#define SAFE_ALOC_GPU(X,Y) hipSafeCall(hipMalloc((void**)&(X),(Y)));
//CPU host free
#define SAFE_FREE_HOST(X) if ((X) != NULL){ \
               hipSafeCall( hipHostFree((X))); \
               (X) = NULL;}
//GPU device free
#define SAFE_FREE_GPU(X) if ((X) != NULL){ \
               hipSafeCall( hipFree((X))); \
               (X) = NULL;}
/***********************************************
 * GPU Timer based on HIP event
 ***********************************************/
//GPU Timer object definition
typedef struct GPU_Timer
{
    GPU_Timer()
    {
        hipSafeCall( hipEventCreate(&this->start) );
        hipSafeCall( hipEventCreate(&this->stop) );
    }
    ~GPU_Timer()
    {
        hipSafeCall( hipEventDestroy(this->start));
        hipSafeCall( hipEventDestroy(this->stop));
    }
    void start_timer() { hipSafeCall( hipEventRecord(this->start) ); }
    void stop_timer() { hipSafeCall( hipEventRecord(this->stop) ); }
    double measure()
    {
        hipSafeCall( hipEventSynchronize(this->stop) );
        float millisconds = 0;
        hipSafeCall(hipEventElapsedTime(&millisconds, this->start, this->stop) ); 
        return (double)millisconds;
    }
    hipEvent_t start;
    hipEvent_t stop;
} gpu_timer;

#else //CPU
//==================================== CPU =======================================
/***********************************************
 * Memory allocation and free
 ***********************************************/
//CPU allocation
#define SAFE_ALOC_HOST(X,Y) (*(void**)(&(X))) = (void*)malloc(Y); 
//CPU free
#define SAFE_FREE_HOST(X) if ((X) != NULL) { \
               free((X)); \
               (X) = NULL;}
#endif
//==================================== Common =======================================
/***********************************************
 * Error Checking:
 ***********************************************/
// Checking null pointer
#define CHECK_NULL_POINTER(X) __checkNullPointer( __FILE__, __LINE__, (void**)&(X))
inline void __checkNullPointer( const char *file, const int line, void** ptr)
{
    if ((*ptr) == NULL)
    {
        fprintf( stderr, "Error: NULL pointer at %s:%i.\n", file, line);
        exit(-1);
    }
}
/***********************************************
 * CPU Timer based on Linux sys/time.h
 ***********************************************/
//CPU timer
double get_cpu_timer()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    //get current timestamp in milliseconds
    return (double)tp.tv_sec * 1e3 + (double)tp.tv_usec * 1e-3;
}
//CPU timer object definition
typedef struct CPU_Timer
{
    CPU_Timer() { start = stop = 0.0; }
    void start_timer() { start = get_cpu_timer(); }
    void stop_timer() { stop = get_cpu_timer(); }
    double measure() { double millisconds = stop - start; return millisconds; }
    double start;
    double stop;
} cpu_timer;
/***********************************************
 * Printing
 ***********************************************/
//print a binary number
void print_binary(IdxType v, int width)
{
    for (int i=width-1; i>=0; i--) putchar('0' + ((v>>i)&1));
}
//print measurement results for n repetitions
void print_measurement(IdxType* res_state, IdxType n_qubits, int repetition)
{
    assert(res_state != NULL);
    printf("\n===============  Measurement (tests=%d) ================\n", repetition);
    for (int i=0; i<repetition; i++)
    {
        printf("Test-%d: ",i);
        print_binary(res_state[i], n_qubits);
        printf("\n");
    }
}
/***********************************************
 * Runtime:
 ***********************************************/
//Swap two pointers
inline void swap_pointers(ValType** pa, ValType** pb)
{
    ValType* tmp = (*pa); (*pa) = (*pb); (*pb) = tmp;
}
//Verify whether a number is power of 2
inline bool is_power_of_2(int x)
{
    return (x>0 && !(x&(x-1)));
}
//Random value between 0 and 1
inline ValType randomval()
{
    return (ValType)std::rand()/(ValType)RAND_MAX;
}
}; //namespace DMSim
#endif
