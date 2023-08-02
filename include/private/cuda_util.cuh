#pragma once
#include <stdio.h>
#include <sys/time.h>
#include <assert.h>
#include <iostream>

/***********************************************
 * CUDA specific attributes
 ***********************************************/
#define THREADS_CTA_CUDA 256

/***********************************************
 * Macros for CUDA
 ***********************************************/
#define LOCAL_G_CUDA(arr, i) arr[(i)]
#define LOCAL_P_CUDA(arr, i, val) arr[(i)] = val;
#define BARR_CUDA grid.sync();

// CLUSTER BASED
#define PGAS_P(arr, i, val) nvshmem_double_p(&(arr)[(i) & ((m_gpu)-1)], (val), ((i) >> (lg2_m_gpu)))
#define PGAS_G(arr, i) nvshmem_double_g(&(arr)[(i) & ((m_gpu)-1)], ((i) >> (lg2_m_gpu)))
#define BARR_NVSHMEM                         \
    if (threadIdx.x == 0 && blockIdx.x == 0) \
        nvshmem_barrier_all();               \
    grid.sync();

#define LOCAL_G_CUDA_MPI(arr, i) arr[(i) & (m_gpu - 1)]
#define LOCAL_P_CUDA_MPI(arr, i, val) arr[(i) & (m_gpu - 1)] = val;

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
#define cudaCheckError() __cudaCheckError(__FILE__, __LINE__)
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
    if (cudaSuccess != err)
    {
        fprintf(stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
                file, line, cudaGetErrorString(err));
        exit(-1);
    }
#endif
    return;
}
// Error checking for NVSHMEM call
#define NVSHMEM_CHECK(stmt)                                            \
    do                                                                 \
    {                                                                  \
        int result = (stmt);                                           \
        if (NVSHMEMX_SUCCESS != result)                                \
        {                                                              \
            fprintf(stderr, "[%s:%d] nvshmem failed with error %d \n", \
                    __FILE__, __LINE__, result);                       \
            exit(-1);                                                  \
        }                                                              \
    } while (0)
/***********************************************
 * Memory allocation and free
 ***********************************************/
// CPU host allocation with pinned memory
#define SAFE_ALOC_HOST_CUDA(X, Y) cudaSafeCall(cudaMallocHost((void **)&(X), (Y)));
// GPU device allocation
#define SAFE_ALOC_GPU(X, Y) cudaSafeCall(cudaMalloc((void **)&(X), (Y)));
// CPU host free for pinned memory
#define SAFE_FREE_HOST_CUDA(X)           \
    if ((X) != NULL)                     \
    {                                    \
        cudaSafeCall(cudaFreeHost((X))); \
        (X) = NULL;                      \
    }
// GPU device free
#define SAFE_FREE_GPU(X)             \
    if ((X) != NULL)                 \
    {                                \
        cudaSafeCall(cudaFree((X))); \
        (X) = NULL;                  \
    }
/***********************************************
 * GPU Timer based on CUDA event
 ***********************************************/
typedef struct GPU_Timer
{
    GPU_Timer()
    {
        cudaSafeCall(cudaEventCreate(&this->start));
        cudaSafeCall(cudaEventCreate(&this->stop));
    }
    ~GPU_Timer()
    {
        cudaSafeCall(cudaEventDestroy(this->start));
        cudaSafeCall(cudaEventDestroy(this->stop));
    }
    void start_timer() { cudaSafeCall(cudaEventRecord(this->start)); }
    void stop_timer() { cudaSafeCall(cudaEventRecord(this->stop)); }
    double measure()
    {
        cudaSafeCall(cudaEventSynchronize(this->stop));
        float millisconds = 0;
        cudaSafeCall(cudaEventElapsedTime(&millisconds, this->start, this->stop));
        return (double)millisconds;
    }
    cudaEvent_t start;
    cudaEvent_t stop;
} gpu_timer;