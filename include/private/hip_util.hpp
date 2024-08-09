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
#define BARR_ROC_SHMEM                         \
    if (threadIdx.x == 0 && blockIdx.x == 0) \
        roc_shmem_ctx_wg_barrier_all(*p_ctx);               \
    grid.sync();

#define LOCAL_G_HIP_MPI(arr, i) arr[(i) & (m_gpu - 1)]
#define LOCAL_P_HIP_MPI(arr, i, val) arr[(i) & (m_gpu - 1)] = val;


#define PGAS_OG(arr, i) roc_shmem_ctx_double_g(*p_ctx, &(arr)[(i) & ((m_gpu) - 1)], ((i) >> (lg2_m_gpu)))
#define PGAS_OP(arr, i, val) roc_shmem_ctx_double_p(*p_ctx, &(arr)[(i) & ((m_gpu) - 1)], (val), ((i) >> (lg2_m_gpu)))

#define PGAS_G(arr, i)  ((((i)>>(lg2_m_gpu))==i_proc) ? (LOCAL_G_HIP_MPI(arr, i)) : (PGAS_OG(arr, i)));
#define PGAS_P(arr, i, val) if (((i)>>(lg2_m_gpu))==i_proc) { \
    LOCAL_P_HIP_MPI(arr, i, val);} else { PGAS_OP(arr, i, val); roc_shmem_ctx_quiet(*p_ctx);}


#define PGAS_GG(dst,src,i) if (((i)>>(lg2_m_gpu))==i_proc) { \
    *(dst) = src[(i)&(m_gpu-1)];\
} else { \
    IdxType i_gpu = ((i)>>lg2_m_gpu); \
    IdxType i_idx = ((i)&(m_gpu-1)); \
    roc_shmem_ctx_getmem(*p_ctx, &m_imag[i_idx], &src[i_idx], sizeof(ValType), i_gpu);\
    *(dst) = m_imag[i_idx]; } 
    

#define PGAS_PP(dst,i,val) if (((i)>>(lg2_m_gpu))==i_proc) { \
    dst[(i)&(m_gpu-1)]=val;\
} else { \
    IdxType i_gpu = ((i)>>lg2_m_gpu); \
    IdxType i_idx = ((i)&(m_gpu-1)); \
    m_imag[i_idx] = val; \
    roc_shmem_ctx_putmem(*p_ctx, &dst[i_idx], &m_imag[i_idx], sizeof(ValType), i_gpu);\
    roc_shmem_ctx_quiet(*p_ctx);}



                        

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

__global__ void gpu_exp_z(const double *sv_real, const double *sv_imag, double *result, const unsigned long long dim, const int offset)
{
    extern __shared__ double sdata[]; // Shared memory to store partial sums
    int tid = threadIdx.x;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    double temp_sum = 0.0;

    // Each thread computes its part
    if (index < dim)
    {
        bool parityVal = parity(index + offset);
        double magnitude_sq = sv_real[index] * sv_real[index] + sv_imag[index] * sv_imag[index];
        temp_sum += (parityVal ? -1.0 : 1.0) * magnitude_sq;
    }

    // Store the computed value in shared memory and synchronize
    sdata[tid] = temp_sum;
    __syncthreads();

    // Reduction in shared memory
    for (int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }

    // Only the thread 0 of each block writes back to global memory
    if (tid == 0)
    {
        atomicAdd(result, sdata[0]);
    }
}

__global__ void gpu_exp_z_bits(const size_t *in_bits, size_t in_bits_size, const double *sv_real, const double *sv_imag, double *result, const unsigned long long dim, const int offset)
{
    extern __shared__ double sdata[]; // Shared memory to store partial sums
    int tid = threadIdx.x;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    double temp_sum = 0.0;

    // Each thread computes its part
    if (index < dim)
    {
        bool parityVal = hasEvenParity(index + offset, in_bits, in_bits_size);

        double magnitude_sq = sv_real[index] * sv_real[index] + sv_imag[index] * sv_imag[index];
        temp_sum += (parityVal ? -1.0 : 1.0) * magnitude_sq;
    }

    // Store the computed value in shared memory and synchronize
    sdata[tid] = temp_sum;
    __syncthreads();

    // Reduction in shared memory
    for (int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }

    // Only the thread 0 of each block writes back to global memory
    if (tid == 0)
    {
        atomicAdd(result, sdata[0]);
    }
}
