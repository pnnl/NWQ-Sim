#pragma once

#include "vqe_state.hpp"
#include "backendManager.hpp"
#include "svsim_vqe/sv_cpu_vqe.hpp"


#ifdef MPI_ENABLED
#include "svsim_vqe/sv_mpi_vqe.hpp"
#endif

#ifdef CUDA_ENABLED
#include "svsim_vqe/sv_cuda_vqe.cuh"
#endif

#ifdef CUDA_MPI_ENABLED
#include "svsim_vqe/sv_cuda_mpi_vqe.cuh"
#endif


class VQEBackendManager: public BackendManager
{
public:
    
    static std::shared_ptr<NWQSim::VQE::VQEState> create_vqe_solver(std::string backend, 
                                                                    const std::string& config,
                                                                    std::shared_ptr<NWQSim::VQE::Ansatz> a, 
                                                                    std::shared_ptr<NWQSim::VQE::Hamiltonian> h, 
                                                                    nlopt::algorithm optimizer_algorithm,
                                                                    NWQSim::VQE::Callback _callback,
                                                                    NWQSim::IdxType seed = 0,
                                                                    NWQSim::VQE::OptimizerSettings opt_settings = NWQSim::VQE::OptimizerSettings()) {
        // Convert to uppercase
        std::transform(backend.begin(), backend.end(), backend.begin(),
                       [](unsigned char c)
                       { return std::toupper(c); });
        if (backend == "CPU")
        {
            return std::make_shared<NWQSim::VQE::SV_CPU_VQE>(a, h, optimizer_algorithm, _callback, config, seed, opt_settings);
        }

#ifdef MPI_ENABLED
        else if (backend == "MPI")
        {
            return std::make_shared<NWQSim::VQE::SV_MPI_VQE>(a, h, optimizer_algorithm, _callback, config, seed, opt_settings);
        }
#endif

#ifdef CUDA_ENABLED
        else if (backend == "NVGPU")
        {
            return std::make_shared<NWQSim::VQE::SV_CUDA_VQE>(a, h, optimizer_algorithm, _callback, config, seed, opt_settings);
        }
#endif

#ifdef CUDA_MPI_ENABLED
        else if (backend == "NVGPU_MPI")
        {
            return std::make_shared<NWQSim::VQE::SV_CUDA_MPI_VQE>(a, h, optimizer_algorithm, _callback, config, seed, opt_settings);
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
