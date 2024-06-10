#pragma once

#include <assert.h>
#include <iostream>
#include <sys/time.h>
#include <vector>
#include <chrono>
#include <iomanip>

/***********************************************
 * Constant configuration:
 ***********************************************/
/* Constant value of PI */
#define PI 3.14159265358979323846
/* Constant value of sqrt(2) */
#define S2I 0.70710678118654752440
/* Constant value of 0.5 */
#define HALF 0.5

/* Error bar for purity check and other error check */
#define ERROR_BAR (1e-3)

#define PRINT_PROGRESS_BAR

#define PURITY_CHECK

namespace NWQSim
{
    /* Basic data type for indices */
    using IdxType = long long int;
    /* Basic data type for value */
    using ValType = double;

    inline std::string formatDuration(std::chrono::seconds input_seconds)
    {
        using namespace std::chrono;
        hours hrs = duration_cast<hours>(input_seconds % 24h);
        minutes mins = duration_cast<minutes>(input_seconds % 1h);
        seconds secs = duration_cast<seconds>(input_seconds % 1min);
        return (hrs.count() < 10 ? "0" : "") + std::to_string(hrs.count()) + ":" +
               (mins.count() < 10 ? "0" : "") + std::to_string(mins.count()) + ":" +
               (secs.count() < 10 ? "0" : "") + std::to_string(secs.count());
    }

    inline void printProgressBar(int current, int total, std::chrono::time_point<std::chrono::steady_clock> start_time)
    {
        const int barWidth = 50;
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
        auto estimated = (total > current && current > 0) ? elapsed * total / current : elapsed;
        auto remaining = estimated - elapsed;

        std::cout << "\033[1;34m[";
        int pos = barWidth * current / total;
        for (int i = 0; i < barWidth; ++i)
        {
            if (i < pos)
                std::cout << "=";
            else if (i == pos)
                std::cout << ">";
            else
                std::cout << " ";
        }

        std::cout << "] "
                  << "\033[1;32m" << int(current * 100.0 / total) << " % "
                  << "\033[1;33mElapsed: " << formatDuration(std::chrono::seconds(elapsed))
                  << " Estimated: " << formatDuration(std::chrono::seconds(estimated))
                  << " Remaining: " << formatDuration(std::chrono::seconds(remaining)) << "\033[0m  \r";
        std::cout.flush();
    }
    inline bool hasEvenParity(unsigned long long x, const std::vector<size_t> &in_qubitIndices)
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
    }
    inline bool hasEvenParity(unsigned long long x, size_t maxsize)
    {
        size_t count = 0;
        for (size_t bitIdx = 0; bitIdx < maxsize; bitIdx++)
        {
            if (x & (1ULL << bitIdx))
            {
                count++;
            }
        }
        return (count % 2) == 0;
    }

    /***********************************************
     * CPU Timer based on Linux sys/time.h
     ***********************************************/
    // CPU timer
    inline double get_cpu_timer()
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
    inline void print_binary(IdxType v, int width)
    {
        for (int i = width - 1; i >= 0; i--)
            putchar('0' + ((v >> i) & 1));
    }
    // print measurement results for n repetitions
    inline void print_measurement(IdxType *res_state, IdxType n_qubits, int repetition)
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
