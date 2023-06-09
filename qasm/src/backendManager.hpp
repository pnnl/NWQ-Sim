#pragma once

#include "public/state.hpp"
#include "public/util.hpp"

#include "svsim/sv_cpu.hpp"

#ifdef OMP_ENABLED
#include "svsim/sv_omp.hpp"
#endif

#ifdef MPI_ENABLED
#include "svsim/sv_mpi.hpp"
#endif

#ifdef CUDA_ENABLED
#include "svsim/sv_nvgpu.cuh"
#endif

#ifdef CUDA_MPI_ENABLED
#include "svsim/sv_nvgpu_mpi.cuh"
#endif

#include <iostream>
#include <memory>

#include <string>
#include <algorithm>
#include <cctype>

class BackendManager
{
public:
    // Variadic template to act like printf
    template <typename... Args>
    static void safe_print(const char *format, Args... args)
    {
#ifdef MPI_ENABLED
        int flag;
        MPI_Initialized(&flag);
        if (flag)
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0)
            {
                printf(format, args...);
            }
        }
        else
        {
            printf(format, args...);
        }
#else
        printf(format, args...);
#endif
    }

    static void print_available_backends()
    {
        safe_print("Available backends:\n");
        safe_print("- CPU\n");
#ifdef OMP_ENABLED
        safe_print("- OpenMP\n");
#endif

#ifdef MPI_ENABLED
        safe_print("- MPI\n");
#endif

#ifdef CUDA_ENABLED
        safe_print("- NVGPU\n");
#endif

#ifdef CUDA_MPI_ENABLED
        safe_print("- NVGPU_MPI\n");
#endif
    }

    static std::shared_ptr<NWQSim::QuantumState> create_state(std::string backend, NWQSim::IdxType numQubits)
    {

        if (backend == "CPU")
        {
            return std::make_shared<NWQSim::SV_CPU>(numQubits);
        }

#ifdef OMP_ENABLED
        else if (backend == "OPENMP")
        {
            return std::make_shared<NWQSim::SV_OMP>(numQubits);
        }
#endif

#ifdef MPI_ENABLED
        else if (backend == "MPI")
        {
            return std::make_shared<NWQSim::SV_MPI>(numQubits);
        }
#endif

#ifdef CUDA_ENABLED
        else if (backend == "NVGPU")
        {
            return std::make_shared<NWQSim::SV_NVGPU>(numQubits);
        }
#endif

#ifdef CUDA_MPI_ENABLED
        else if (backend == "NVGPU_MPI")
        {
            return std::make_shared<NWQSim::SV_NVGPU_MPI>(numQubits);
        }
#endif
        else
        {
            safe_print("Invalid backend name. Please use one of the available backends. (Case insensitive)\n");
            print_available_backends();
            exit(1);
        }
    }
};
