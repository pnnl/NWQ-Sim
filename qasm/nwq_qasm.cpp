#include "backendManager.hpp"

#include <memory>
#include <string>

#include "state.hpp"
#include "util.hpp"

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " [backend]\n";
        return 1;
    }

    std::string backendName(argv[1]);
    BackendType backendType = BackendType::CPU;

    if (backendName == "CPU")
    {
        backendType = BackendType::CPU;
    }
    else if (backendName == "OpenMP")
    {
        backendType = BackendType::OpenMP;
    }
    else if (backendName == "MPI")
    {
        backendType = BackendType::MPI;
    }
    else if (backendName == "NVGPU")
    {
        backendType = BackendType::NVGPU;
    }
    else if (backendName == "NVGPU_MPI")
    {
        backendType = BackendType::NVGPU_MPI;
    }
    else
    {
        std::cerr << "Invalid backend name\n";
        return 1;
    }

// If MPI or NVSHMEM backend, initialize MPI
#ifdef MPI_ENABLED
    if (backendType == BackendType::MPI || backendType == BackendType::NVGPU_MPI)
    {
        MPI_Init(&argc, &argv);
    }
#endif

    // Create the backend
    std::unique_ptr<NWQSim::QuantumState> backend = BackendManager::createBackend(backendType, 10);
    if (!backend)
    {
        std::cerr << "Failed to create backend\n";
        return 1;
    }

    // Use the backend...
    BackendManager::safePrint("Created backend: %s\n", backendName.c_str());

// Finalize MPI if necessary
#ifdef MPI_ENABLED
    if (backendType == BackendType::MPI || backendType == BackendType::NVGPU_MPI)
    {
        MPI_Finalize();
    }
#endif

    return 0;
}
