#pragma once

#include <hip/hip_runtime.h>

/***********************************************
 * AMD GPU specific attributes
 ***********************************************/
#define THREADS_CTA_HIP 1024

#define LOCAL_G_HIP(arr, i) arr[(i)]
#define LOCAL_P_HIP(arr, i, val) arr[(i)] = val;
#define BARR_HIP grid.sync();

// CLUSTER BASED

#define PGAS_P(arr, i, val) roc_shmem_putmem(&(arr)[(i) & ((m_gpu) - 1)], (val), ((i) >> (lg2_m_gpu)))
#define PGAS_G(arr, i) roc_nvshmem_getmem(&(arr)[(i) & ((m_gpu) - 1)], ((i) >> (lg2_m_gpu)))
#define BARR_ROC_SHMEM                         \
    if (threadIdx.x == 0 && blockIdx.x == 0) \
        roc_shmem_ctx_wg_barrier_all(*p_ctx);               \
    grid.sync();



/***********************************************
 * Error Checking:
 ***********************************************/
// Error checking for HIP API call
#define hipSafeCall(err) __hipSafeCall(err, __FILE__, __LINE__)
inline void __hipSafeCall(hipError_t err, const char *file, const int line)
{
#ifdef GPU_ERROR_CHECK
    if (hipSuccess != err)
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
        exit(-1);
    }
    // Expensive checking
    err = hipDeviceSynchronize();
    if (hipSuccess != err)
    {
        fprintf(stderr, "hipCheckError() with sync failed at %s:%i : %s\n",
                file, line, hipGetErrorString(err));
        exit(-1);
    }
#endif
    return;
}
/***********************************************
 * Memory allocation and free
 ***********************************************/
// CPU host allocation
#define SAFE_ALOC_HOST_HIP(X, Y) hipSafeCall(hipHostMalloc((void **)&(X), (Y)));
// GPU device allocation
#define SAFE_ALOC_GPU_HIP(X, Y) hipSafeCall(hipMalloc((void **)&(X), (Y)));
// CPU host free
#define SAFE_FREE_HOST_HIP(X)          \
    if ((X) != NULL)                   \
    {                                  \
        hipSafeCall(hipHostFree((X))); \
        (X) = NULL;                    \
    }
// GPU device free
#define SAFE_FREE_GPU_HIP(X)       \
    if ((X) != NULL)               \
    {                              \
        hipSafeCall(hipFree((X))); \
        (X) = NULL;                \
    }
/***********************************************
 * GPU Timer based on HIP event
 ***********************************************/
// GPU Timer object definition
typedef struct GPU_Timer
{
    GPU_Timer()
    {
        hipSafeCall(hipEventCreate(&this->start));
        hipSafeCall(hipEventCreate(&this->stop));
    }
    ~GPU_Timer()
    {
        hipSafeCall(hipEventDestroy(this->start));
        hipSafeCall(hipEventDestroy(this->stop));
    }
    void start_timer() { hipSafeCall(hipEventRecord(this->start)); }
    void stop_timer() { hipSafeCall(hipEventRecord(this->stop)); }
    double measure()
    {
        hipSafeCall(hipEventSynchronize(this->stop));
        float millisconds = 0;
        hipSafeCall(hipEventElapsedTime(&millisconds, this->start, this->stop));
        return (double)millisconds;
    }
    hipEvent_t start;
    hipEvent_t stop;
} gpu_timer;

/***********************************************
 * VQE Related Functions
 ***********************************************/

__device__ bool hasEvenParity(unsigned long long x, const size_t *in_bits, const size_t in_bits_size)
{
    size_t count = 0;
    for (size_t i = 0; i < in_bits_size; ++i)
    {
        if (x & (1ULL << in_bits[i]))
        {
            count++;
        }
    }
    return (count % 2) == 0;
}

__device__ double parity(unsigned long long num)
{
    num ^= num >> 32;
    num ^= num >> 16;
    num ^= num >> 8;
    num ^= num >> 4;
    num &= 0xf;
    return (0x6996 >> num) & 1;
}

__global__ void gpu_exp_z_bits(const size_t *in_bits, size_t in_bits_size, const double *sv_real, const double *sv_imag, double *result, const unsigned long long dim)
{
    unsigned long long idx = threadIdx.x + blockIdx.x * blockDim.x + blockIdx.y * blockDim.x * gridDim.x;

    if (idx < dim)
    {
        double res = (hasEvenParity(idx, in_bits, in_bits_size) ? 1.0 : -1.0) *
                     (sv_real[idx] * sv_real[idx] + sv_imag[idx] * sv_imag[idx]);

        atomicAdd(result, res);
    }
}

__global__ void gpu_exp_z(const double *sv_real, const double *sv_imag, double *result, const unsigned long long dim)
{
    unsigned long long idx = threadIdx.x + blockIdx.x * blockDim.x + blockIdx.y * blockDim.x * gridDim.x;

    if (idx < dim)
    {
        bool parityVal = parity(idx);
        double res = (parityVal ? -1.0 : 1.0) * (sv_real[idx] * sv_real[idx] + sv_imag[idx] * sv_imag[idx]);

        atomicAdd(result, res);
    }
}

    
