#pragma once

#include "state.hpp"
#include "nwq_util.hpp"

#include "svsim/sv_cpu.hpp"
#include "dmsim/dm_cpu.hpp"
#include "stabsim/stab_cpu.hpp"

#ifdef OMP_ENABLED
#include "svsim/sv_omp.hpp"
#endif

#ifdef MPI_ENABLED
#include "svsim/sv_mpi.hpp"
#endif

#ifdef CUDA_ENABLED
#include "svsim/sv_cuda.cuh"
#include "dmsim/dm_cuda.cuh"
#include "stabsim/stab_cuda.cuh"
#endif

#ifdef CUDA_MPI_ENABLED
#include "svsim/sv_cuda_mpi.cuh"
#include "dmsim/dm_cuda_mpi.cuh"
#endif

#ifdef HIP_ENABLED
#include "svsim/sv_hip.hpp"
#include "dmsim/dm_hip.hpp"
#endif

#include <iostream>
#include <memory>

#include <string>
#include <algorithm>
#include <cctype>

// TNSIM backends
#ifdef ITENSOR_ENABLED
#include "tnsim/tn_itensor.hpp"
#endif

// Both itensor and TAMM have Error (may be a better way to do this but this works)
#ifdef Error
#undef Error
#endif

#ifdef TAMM_ENABLED
#include "tnsim/tn_tamm.hpp"
#endif

class BackendManager
{
public:
    static void print_available_backends()
    {
        NWQSim::safe_print("Available backends:\n");
        NWQSim::safe_print("- CPU\n");
#ifdef OMP_ENABLED
        NWQSim::safe_print("- OpenMP\n");
#endif

#ifdef MPI_ENABLED
        NWQSim::safe_print("- MPI\n");
#endif

#ifdef CUDA_ENABLED
        NWQSim::safe_print("- NVGPU\n");
#endif

#ifdef CUDA_MPI_ENABLED
        NWQSim::safe_print("- NVGPU_MPI\n");
#endif

#ifdef HIP_ENABLED
        NWQSim::safe_print("- AMDGPU\n");
#endif

#ifdef TAMM_ENABLED
        NWQSim::safe_print("- TN_TAMM_CPU");
        NWQSim::safe_print("- TN_TAMM_GPU");
#endif
    }

    static std::shared_ptr<NWQSim::QuantumState> create_state(std::string backend, NWQSim::IdxType numQubits, std::string simulator_method = "SV", int max_dim = 100, double sv_cutoff = 0.0)
    {
        // Convert to uppercase
        std::transform(backend.begin(), backend.end(), backend.begin(),
                       [](unsigned char c)
                       { return std::toupper(c); });

        // Convert to uppercase
        std::transform(simulator_method.begin(), simulator_method.end(), simulator_method.begin(),
                       [](unsigned char c)
                       { return std::toupper(c); });
        if (backend == "CPU")
        {
            if (simulator_method == "SV")
                return std::make_shared<NWQSim::SV_CPU>(numQubits);
            else if (simulator_method == "DM")
                return std::make_shared<NWQSim::DM_CPU>(numQubits);
            else if (simulator_method == "STAB")
                return std::make_shared<NWQSim::STAB_CPU>(numQubits);
#ifdef ITENSOR_ENABLED
            else if (simulator_method == "TN")
                return std::make_shared<NWQSim::TN_ITENSOR>(numQubits, max_dim, sv_cutoff);
#endif
            else
            {
                NWQSim::safe_print("Invalid simulator method name: %s. Please use one of the available methods: SV, DM, STAB. (Case insensitive)\n", simulator_method.c_str());
                exit(1);
            }
        }

#ifdef TAMM_ENABLED
        else if (backend.rfind("TN_TAMM", 0) == 0)
        {
            return std::make_shared<NWQSim::TN_TAMM>(numQubits, max_dim, sv_cutoff, backend);
        }
#endif

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
            if (simulator_method == "SV")
                return std::make_shared<NWQSim::SV_CUDA>(numQubits);
            else if (simulator_method == "DM")
                return std::make_shared<NWQSim::DM_CUDA>(numQubits);
            else if (simulator_method == "STAB")
                return std::make_shared<NWQSim::STAB_CUDA>(numQubits);
            else
            {
                NWQSim::safe_print("Invalid simulator method name: %s. Please use one of the available methods: SV, DM. (Case insensitive)\n", simulator_method.c_str());
                exit(1);
            }
        }
#endif

#ifdef HIP_ENABLED
        else if (backend == "AMDGPU")
        {
            if (simulator_method == "SV")
                return std::make_shared<NWQSim::SV_HIP>(numQubits);
            else
                return std::make_shared<NWQSim::DM_HIP>(numQubits);
        }
#endif

#ifdef CUDA_MPI_ENABLED
        else if (backend == "NVGPU_MPI")
        {
            if (simulator_method == "SV")
                return std::make_shared<NWQSim::SV_CUDA_MPI>(numQubits);
            else
                return std::make_shared<NWQSim::DM_CUDA_MPI>(numQubits);
        }
#endif
        else if (backend == "LIST")
        {
            print_available_backends();
            exit(0);
        }
        else
        {
            NWQSim::safe_print("Invalid backend name: %s. Please use one of the available backends. (Case insensitive)\n", backend.c_str());
            print_available_backends();
            exit(1);
        }
    }
};
