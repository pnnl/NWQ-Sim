#pragma once

#include <assert.h>
#include <iostream>
#include <sys/time.h>
#include <vector>

/***********************************************
 * Constant configuration:
 ***********************************************/
/* Constant value of PI */
#define PI 3.14159265358979323846
/* Constant value of sqrt(2) */
#define S2I 0.70710678118654752440
/* Constant value of 0.5 */
#define HALF 0.5

#define DEFAULT_REPETITIONS 1024

/* Error bar for purity check and other error check */
#define ERROR_BAR (1e-3)

namespace NWQSim
{
    /* Basic data type for indices */
    using IdxType = long long int;
    /* Basic data type for value */
    using ValType = double;

    bool hasEvenParity(unsigned long long x, const std::vector<size_t> &in_qubitIndices)
    {
        size_t count = 0;
        for (const auto &bitIdx : in_qubitIndices)
        {
            if (x & (1ULL << bitIdx))
            {
                count++;
            }
        }
        return (count % 2) == 0;
    };

    /***********************************************
     * CPU Timer based on Linux sys/time.h
     ***********************************************/
    // CPU timer
    double get_cpu_timer()
    {
        struct timeval tp;
        gettimeofday(&tp, NULL);
        // get current timestamp in milliseconds
        return (double)tp.tv_sec * 1e3 + (double)tp.tv_usec * 1e-3;
    }
    // CPU timer object definition
    typedef struct CPU_Timer
    {
        CPU_Timer() { start = stop = 0.0; }
        void start_timer() { start = get_cpu_timer(); }
        void stop_timer() { stop = get_cpu_timer(); }
        double measure()
        {
            double millisconds = stop - start;
            return millisconds;
        }
        double start;
        double stop;
    } cpu_timer;
    /***********************************************
     * Printing
     ***********************************************/
    // print a binary number
    void print_binary(IdxType v, int width)
    {
        for (int i = width - 1; i >= 0; i--)
            putchar('0' + ((v >> i) & 1));
    }
    // print measurement results for n repetitions
    void print_measurement(IdxType *res_state, IdxType n_qubits, int repetition)
    {
        assert(res_state != NULL);
        printf("\n===============  Measurement (tests=%d) ================\n", repetition);
        for (int i = 0; i < repetition; i++)
        {
            printf("Test-%d: ", i);
            print_binary(res_state[i], n_qubits);
            printf("\n");
        }
    }
    /***********************************************
     * Runtime:
     ***********************************************/
    // Swap two pointers
    inline void swap_pointers(ValType **pa, ValType **pb)
    {
        ValType *tmp = (*pa);
        (*pa) = (*pb);
        (*pb) = tmp;
    }
    // Verify whether a number is power of 2
    inline bool is_power_of_2(int x)
    {
        return (x > 0 && !(x & (x - 1)));
    }
    // Random value between 0 and 1
    inline ValType randomval()
    {
        return (ValType)std::rand() / (ValType)RAND_MAX;
    }
} // namespace NWQSim
