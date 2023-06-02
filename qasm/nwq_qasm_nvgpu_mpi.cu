#include <time.h>
#include <iostream>
#include <fstream>
// stats calculation
#include <algorithm>
#include <numeric>
#include <fstream>
#include <string>

#include <mpi.h>
/********** UPDATE THE INCLUDE PATH FOR LOCAL HEADER FILE HERE ************/
#include "src/qasm_parser.hpp"
#include "src/parser_util.hpp"
/**************************************************************************/

/********** LOAD THE SIMULATION BACKEND ************/
#include "../include/state.hpp"              // base class for state
#include "../include/svsim/sv_nvgpu_mpi.cuh" // derived class for state
/**************************************************************************/

using namespace NWQSim;
ValType pass_threshold = 0.99;

ValType run_brnchmark(IdxType index, IdxType total_shots, bool is_basis, int rank, bool print = false);
Circuit *get_bv_circuit(IdxType n_qubits, IdxType shots);

int main(int argc, char **argv)
{
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    IdxType total_shots = 16384;
    bool run_with_basis = false;

    if (cmdOptionExists(argv, argv + argc, "--shots"))
    {
        const char *shots_str = getCmdOption(argv, argv + argc, "--shots");
        total_shots = stoi(shots_str);
    }

    if (cmdOptionExists(argv, argv + argc, "--basis"))
    {
        run_with_basis = true;
    }

    if (cmdOptionExists(argv, argv + argc, "--bv"))
    {
        const char *bv_qubit_str = getCmdOption(argv, argv + argc, "--bv");
        IdxType bv_qubits = stoi(bv_qubit_str);
        Circuit *circuit = get_bv_circuit(bv_qubits, total_shots);

        QuantumState *state = new SV_NVGPU_MPI(bv_qubits);
        state->sim(circuit);
        auto results = state->get_results();
        auto counts = outcome_to_dict(results, bv_qubits, total_shots);

        if (rank == 0)
        {
            print_counts(counts, total_shots);
        }

        delete state;
        delete counts;
        delete circuit;
    }

    if (cmdOptionExists(argv, argv + argc, "--app"))
    {
        const char *filename = getCmdOption(argv, argv + argc, "--app");

        qasm_parser parser(filename);

        QuantumState *state = new SV_NVGPU_MPI(parser.num_qubits());
        map<string, IdxType> *counts = parser.execute(state, total_shots);

        if (rank == 0)
            print_counts(counts, total_shots);

        delete state;
        delete counts;
    }

    if (cmdOptionExists(argv, argv + argc, "-t"))
    {
        int benchmark_index = stoi(getCmdOption(argv, argv + argc, "-t"));

        auto fidelity = run_brnchmark(benchmark_index, total_shots, run_with_basis, rank);
        if (rank == 0)
            printf("Fidelity between SVSim and Qiskit Execution: %.4f\n", fidelity);
    }

    if (cmdOptionExists(argv, argv + argc, "-a"))
    {
        bool passed = true;
        for (int benchmark_index = 8; benchmark_index < 36; benchmark_index++)
        {
            ValType fidelity = run_brnchmark(benchmark_index, total_shots, run_with_basis, rank);
            if (fidelity < pass_threshold)
            {
                if (rank == 0)
                    cout << "Benchmark " << benchmark_index << " fidelity: " << fidelity << " Failed!" << endl;
                passed = false;
            }
        }
        if (passed)
        {
            if (rank == 0)
                cout << "TEST PASSED!" << endl;
        }
    }

    MPI_Finalize();

    return 0;
}

Circuit *get_bv_circuit(IdxType n_qubits, IdxType shots)
{
    Circuit *circuit = new Circuit(n_qubits);
    n_qubits--;

    // Create an additional qubit for storing the output
    int output_qubit = n_qubits;

    // Apply a Hadamard gate to all qubits
    for (int i = 0; i <= n_qubits; i++)
    {
        circuit->H(i);
    }

    // Apply a Pauli-X gate to the output qubit
    circuit->X(output_qubit);

    // Apply CX gate for each qubit
    for (int i = 0; i < n_qubits; i++)
    {
        circuit->CX(i, output_qubit);
    }

    // Apply a Hadamard gate to all qubits except the last one
    for (int i = 0; i < n_qubits; i++)
    {
        circuit->H(i);
    }

    // Measure all qubits
    circuit->MA(shots);
    return circuit;
}

ValType run_brnchmark(IdxType index, IdxType total_shots, bool is_basis, int rank, bool print)
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

    QuantumState *state = new SV_NVGPU_MPI(parser.num_qubits());

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

    for (const auto &entry : *svsim_counts)
    {
        auto key = entry.first;
        auto val = entry.second;

        ValType out_prob = val / (ValType)total_shots;
        svsim_norm += out_prob * out_prob;
        svsim_probs.insert({key, out_prob});
    }

    ref_norm = sqrt(ref_norm);
    svsim_norm = sqrt(svsim_norm);

    ValType fidelity = 0;

    for (const auto &entry : svsim_probs)
    {
        auto key = entry.first;
        auto val = entry.second;
        if (ref_probs.count(key))
            fidelity += (val / svsim_norm) * (ref_probs.at(key) / ref_norm);
    }

    if (rank == 0 && print)
    {
        print_counts(svsim_counts, total_shots);
        // state->print_res_sv();
    }

    delete state;
    delete svsim_counts;
    return fidelity;
}