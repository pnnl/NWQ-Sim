#pragma once

#include "state.hpp"
#include "nwq_util.hpp"

#include "svsim/sv_cpu.hpp"
#include "dmsim/dm_cpu.hpp"

#ifdef OMP_ENABLED
#include "svsim/sv_omp.hpp"
#endif

#ifdef MPI_ENABLED
#include "svsim/sv_mpi.hpp"
#endif

#ifdef CUDA_ENABLED
#include "svsim/sv_nvgpu.cuh"
#include "dmsim/dm_nvgpu.cuh"
#endif

#ifdef CUDA_MPI_ENABLED
#include "svsim/sv_nvgpu_mpi.cuh"
#include "dmsim/dm_nvgpu_mpi.cuh"
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

    static std::shared_ptr<NWQSim::QuantumState> create_state(std::string backend, NWQSim::IdxType numQubits, std::string simulator_method = "sv")
    {
        // Convert to uppercase
        std::transform(backend.begin(), backend.end(), backend.begin(),
                       [](unsigned char c)
                       { return std::toupper(c); });

        if (backend == "CPU")
        {
            if (simulator_method == "sv")
                return std::make_shared<NWQSim::SV_CPU>(numQubits);
            else
                return std::make_shared<NWQSim::DM_CPU>(numQubits);
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
            if (simulator_method == "sv")
                return std::make_shared<NWQSim::SV_NVGPU>(numQubits);
            else
                return std::make_shared<NWQSim::DM_NVGPU>(numQubits);
        }
#endif

#ifdef CUDA_MPI_ENABLED
        else if (backend == "NVGPU_MPI")
        {
            if (simulator_method == "sv")
                return std::make_shared<NWQSim::SV_NVGPU_MPI>(numQubits);
            else
                return std::make_shared<NWQSim::DM_NVGPU_MPI>(numQubits);
        }
#endif
        else if (backend == "LIST")
        {
            print_available_backends();
            exit(0);
        }
        else
        {
            safe_print("Invalid backend name: %s. Please use one of the available backends. (Case insensitive)\n", backend.c_str());
            print_available_backends();
            exit(1);
        }
    }
};
