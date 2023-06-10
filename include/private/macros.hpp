#pragma once

#include <stdlib.h>

namespace NWQSim
{
/***********************************************
 * Memory allocation and free
 ***********************************************/
// CPU allocation
#define SAFE_ALOC_HOST(X, Y) (*(void **)(&(X))) = (void *)malloc(Y);
// CPU free
#define SAFE_FREE_HOST(X) \
    if ((X) != NULL)      \
    {                     \
        free((X));        \
        (X) = NULL;       \
    }

//==================================== Common =======================================
/***********************************************
 * Error Checking:
 ***********************************************/
// Checking null pointer
#define CHECK_NULL_POINTER(X) __checkNullPointer(__FILE__, __LINE__, (void **)&(X))
    inline void __checkNullPointer(const char *file, const int line, void **ptr)
    {
        if ((*ptr) == NULL)
        {
            fprintf(stderr, "Error: NULL pointer at %s:%i.\n", file, line);
            exit(-1);
        }
    }

/***********************************************
 * Key Macros CPU
 ***********************************************/
#define PUT(arr, i, val) (arr[(i)] = (val)) 
#define GET(arr, i) (arr[(i)]) 
#define BARR  \
    while (0) \
    {         \
    };

#define BARR_MPI MPI_Barrier(MPI_COMM_WORLD);

    /***********************************************
     * Bitwise Macros
     ***********************************************/

// For C2 and C4 gates
#define DIV2E(x, y) ((x) >> (y))
#define MOD2E(x, y) ((x) & (((IdxType)1 << (y)) - (IdxType)1))
#define EXP2E(x) ((IdxType)1 << (x))
#define SV4IDX(x) (((x >> 1) & 1) * EXP2E(qubit0) + ((x & 1) * EXP2E(qubit1)))
#define SV16IDX(x) (((x >> 3) & 1) * EXP2E(qubit0) + ((x >> 2) & 1) * EXP2E(qubit1) + ((x >> 1) & 1) * EXP2E(qubit2) + ((x & 1) * EXP2E(qubit3)))

// FOR MPI
#define LOCAL_G(arr, i) arr[(i) & (m_cpu - 1)]
#define LOCAL_P(arr, i, val) arr[(i) & (m_cpu - 1)] = val;

} // namespace NWQSim