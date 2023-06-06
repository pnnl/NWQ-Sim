#pragma once

#include "state.hpp"
#include "util.hpp"
#include "svsim/sv_cpu.hpp"

#ifdef OMP_ENABLED
#include "svsim/sv_omp.hpp"
#endif

#ifdef MPI_ENABLED
#include "svsim/sv_mpi.hpp"

#include <mpi.h>
#endif

#ifdef CUDA_ENABLED
#include "svsim/sv_nvgpu.cuh"
#endif

#ifdef CUDA_MPI_ENABLED
#include "svsim/sv_nvgpu_mpi.cuh"
#endif

#include <iostream>
#include <memory>
// #include <fmt/core.h>

enum class BackendType
{
    CPU,
    OpenMP,
    MPI,
    NVGPU,
    NVGPU_MPI
};

class BackendManager
{
public:
    // Variadic template to act like printf
    template <typename... Args>
    static void safePrint(const char *format, Args... args)
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

    static void printAvailableBackends()
    {
        safePrint("Available backends:\n");
        safePrint("- CPU\n");
#ifdef OMP_ENABLED
        safePrint("- OpenMP\n");
#endif

#ifdef MPI_ENABLED
        safePrint("- MPI\n");
#endif

#ifdef CUDA_ENABLED
        safePrint("- NVGPU\n");
#endif

#ifdef NVSHMEM_ENABLED
        safePrint("- NVGPU_MPI\n");
#endif
    }

    static std::unique_ptr<NWQSim::QuantumState> createBackend(BackendType backendType, NWQSim::IdxType numQubits)
    {
        switch (backendType)
        {
        case BackendType::CPU:
            return std::make_unique<NWQSim::SV_CPU>(numQubits);

#ifdef OMP_ENABLED
        case BackendType::OpenMP:
            return std::make_unique<NWQSim::SV_OMP>(numQubits);
#endif

#ifdef MPI_ENABLED
        case BackendType::MPI:
            return std::make_unique<NWQSim::SV_MPI>(numQubits);
#endif

#ifdef CUDA_ENABLED
        case BackendType::NVGPU:
            return std::make_unique<NWQSim::SV_NVGPU>(numQubits);
#endif

#ifdef NVSHMEM_ENABLED
        case BackendType::NVGPU_MPI:
            return std::make_unique<NWQSim::SV_NVGPU_MPI>(numQubits);
#endif

        default:
            safePrint("Invalid backend type. Please use one of the available backends.\n");
            return nullptr;
        }
    }
};
