#pragma once

#include <hip/hip_runtime.h>

/***********************************************
 * AMD GPU specific attributes
 ***********************************************/
#define THREADS_CTA_HIP 256


#define LOCAL_G_HIP(arr, i) arr[(i)]
#define LOCAL_P_HIP(arr, i, val) arr[(i)] = val;
#define BARR_HIP grid.sync();

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