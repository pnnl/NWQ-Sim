#include <memory>
#include <string>
#include <algorithm>
#include <cctype>

/********** UPDATE THE INCLUDE PATH FOR LOCAL HEADER FILE HERE ************/
#include "src/qasm_parser.hpp"
#include "src/parser_util.hpp"
/**************************************************************************/

/********** LOAD THE SIMULATION BACKEND ************/
#include "backendManager.hpp"
#include "state.hpp"
#include "nwq_util.hpp"
/**************************************************************************/

using namespace NWQSim;
ValType pass_threshold = 0.98;
ValType run_brnchmark(std::string backend_name, IdxType index, IdxType total_shots, std::string simulation_method, bool is_basis);

int main(int argc, char **argv)
{
    IdxType total_shots = 16384;
    bool run_with_basis = false;
    std::string backend_name = "CPU";
    std::string simulation_method = "sv";

    if (cmdOptionExists(argv, argv + argc, "-shots"))
    {
        const char *shots_str = getCmdOption(argv, argv + argc, "-shots");
        total_shots = stoi(shots_str);
    }

    if (cmdOptionExists(argv, argv + argc, "-basis"))
    {
        run_with_basis = true;
    }

    if (cmdOptionExists(argv, argv + argc, "-backend"))
    {
        backend_name = std::string(getCmdOption(argv, argv + argc, "-backend"));
    }

    if (cmdOptionExists(argv, argv + argc, "-sim"))
    {
        simulation_method = std::string(getCmdOption(argv, argv + argc, "-sim"));
    }

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
        std::shared_ptr<NWQSim::QuantumState> state = BackendManager::create_state(backend_name, parser.num_qubits(), simulation_method);
        if (!state)
        {
            std::cerr << "Failed to create backend\n";
            return 1;
        }
        state->print_config(simulation_method);
        map<string, IdxType> *counts = parser.execute(state, total_shots);

        std::vector<size_t> in_bits;

        for (int i = 0; i < parser.num_qubits(); ++i)
        {
            in_bits.push_back(i);
        }

        ValType exp_val_z = state->get_exp_z(in_bits);
        ValType exp_val_z_all = state->get_exp_z();
        for (size_t i = 0; i < parser.num_qubits(); ++i)
        {
            in_bits[i] = i;
        }

        if (state->i_proc == 0)
        {
            print_counts(counts, total_shots);
            // fflush(stdout);
            printf("exp-val-z: %f\n", exp_val_z);
            printf("exp-val-z all: %f\n", exp_val_z_all);
        }

        delete counts;
    }

    if (cmdOptionExists(argv, argv + argc, "-t"))
    {
        int benchmark_index = stoi(getCmdOption(argv, argv + argc, "-t"));

        ValType fidelity = run_brnchmark(backend_name, benchmark_index, total_shots, simulation_method, run_with_basis);

        BackendManager::safe_print("Fidelity between NWQSim and Qiskit Execution: %.4f\n", fidelity);
    }

    if (cmdOptionExists(argv, argv + argc, "-a"))
    {
        bool passed = true;
        for (int benchmark_index = 12; benchmark_index < 36; benchmark_index++)
        {
            ValType fidelity = run_brnchmark(backend_name, benchmark_index, total_shots, simulation_method, run_with_basis);
            if (fidelity < pass_threshold)
            {
                BackendManager::safe_print("Benchmark %d fidelity: %.4f Failed!\n", benchmark_index, fidelity);
                passed = false;
            }
        }
        if (passed)
            BackendManager::safe_print("All benchmarks passed!\n");
        else
            BackendManager::safe_print("TESTING FAILED!\n");
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

ValType run_brnchmark(std::string backend_name, IdxType index, IdxType total_shots, std::string simulation_method, bool is_basis)
{
    stringstream ss_file, ss_result;

    if (is_basis)
    {
        ss_file << "../data/benchmarks_basis/circuits/" << index << ".qasm";
        ss_result << "../data/benchmarks_basis/results/" << index << "_result.txt";
    }
    else
    {
        ss_file << "../data/benchmarks/circuits/" << index << ".qasm";
        ss_result << "../data/benchmarks/results/" << index << "_result.txt";
    }

    ifstream resultFile(ss_result.str().c_str());
    if (!resultFile)
    {
        printf("%s\n", "Could not open result file\n");
        return -1;
    }

    qasm_parser parser(ss_file.str().c_str());

    // Create the backend
    std::shared_ptr<NWQSim::QuantumState> state = BackendManager::create_state(backend_name, parser.num_qubits(), simulation_method);
    if (!state)
    {
        std::cerr << "Failed to create backend\n";
        return 1;
    }
    state->print_config(simulation_method);
    map<string, IdxType> *svsim_counts = parser.execute(state, total_shots);

    map<string, ValType> ref_probs;
    map<string, ValType> svsim_probs;
    string line;
    ValType ref_norm = 0, svsim_norm = 0;

    while (getline(resultFile, line))
    {
        vector<string> ref_result = split(line, ' ');
        string outcome = ref_result[0];
        ValType ref_prob = stod(ref_result[1]);
        ref_probs.insert({outcome, ref_prob});
        ref_norm += ref_prob * ref_prob;
    }

    for (auto it = svsim_counts->begin(); it != svsim_counts->end(); ++it)
    {
        const auto &key = it->first;
        const auto &val = it->second;

        ValType out_prob = val / static_cast<ValType>(total_shots);
        svsim_norm += out_prob * out_prob;
        svsim_probs.insert(std::make_pair(key, out_prob));
    }

    ref_norm = sqrt(ref_norm);
    svsim_norm = sqrt(svsim_norm);

    ValType fidelity = 0;

    for (auto it = svsim_probs.begin(); it != svsim_probs.end(); ++it)
    {
        const auto &key = it->first;
        const auto &val = it->second;

        if (ref_probs.count(key))
            fidelity += (val / svsim_norm) * (ref_probs.at(key) / ref_norm);
    }

    return fidelity;
}
