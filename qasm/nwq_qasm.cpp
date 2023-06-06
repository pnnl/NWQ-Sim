#include <memory>
#include <string>
#include <algorithm>
#include <cctype>

/********** UPDATE THE INCLUDE PATH FOR LOCAL HEADER FILE HERE ************/
#include "src/qasm_parser.hpp"
#include "src/parser_util.hpp"
#include "src/backendManager.hpp"
/**************************************************************************/

/********** LOAD THE SIMULATION BACKEND ************/
#include "state.hpp"
#include "util.hpp"
/**************************************************************************/

using namespace NWQSim;
ValType pass_threshold = 0.99;

int main(int argc, char **argv)
{
    IdxType total_shots = 16384;
    // bool run_with_basis = false;
    std::string backend_name = "CPU";

    if (cmdOptionExists(argv, argv + argc, "-shots"))
    {
        const char *shots_str = getCmdOption(argv, argv + argc, "--shots");
        total_shots = stoi(shots_str);
    }

    // if (cmdOptionExists(argv, argv + argc, "-basis"))
    // {
    //     run_with_basis = true;
    // }

    if (cmdOptionExists(argv, argv + argc, "-backend"))
    {
        backend_name = std::string(getCmdOption(argv, argv + argc, "-backend"));
    }

    std::transform(backend_name.begin(), backend_name.end(), backend_name.begin(),
                   [](unsigned char c)
                   { return std::toupper(c); });

// If MPI or NVSHMEM backend, initialize MPI
#ifdef MPI_ENABLED
    if (backend_name == "MPI" || backend_name == "NVGPU_MPI")
    {
        MPI_Init(&argc, &argv);
    }
#endif

    if (cmdOptionExists(argv, argv + argc, "-q"))
    {
        const char *filename = getCmdOption(argv, argv + argc, "-q");

        qasm_parser parser(filename);

        // Create the backend
        std::shared_ptr<NWQSim::QuantumState> state = BackendManager::create_state(backend_name, parser.num_qubits());
        if (!state)
        {
            std::cerr << "Failed to create backend\n";
            return 1;
        }
        map<string, IdxType> *counts = parser.execute(state, total_shots);

        std::vector<size_t> in_bits(parser.num_qubits());
        for (size_t i = 0; i < parser.num_qubits(); ++i)
        {
            in_bits[i] = i;
        }

        BackendManager::safe_print("exp-val-z: %f\n", state->get_exp_z(in_bits));

        delete counts;
    }

// Finalize MPI if necessary
#ifdef MPI_ENABLED
    if (backend_name == "MPI" || backend_name == "NVGPU_MPI")
    {
        MPI_Finalize();
    }
#endif

    return 0;
}
