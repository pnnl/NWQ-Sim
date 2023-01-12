#include <time.h>
#include <iostream>
#include <fstream>
// stats calculation
#include <algorithm>
#include <numeric>
#include <fstream>
#include <string>

/********** UPDATE THE INCLUDE PATH FOR LOCAL HEADER FILE HERE ************/
#include "../qasm/qasm_parser.hpp"
#include "../qasm/parser_util.hpp"
/**************************************************************************/

using namespace NWQSim;
ValType pass_threshold = 0.99;

ValType run_brnchmark(IdxType index, IdxType total_shots, bool is_basis = false);

int main(int argc, char **argv)
{
    IdxType total_shots = 16384;
    bool run_with_basis = false;

    if (cmdOptionExists(argv, argv + argc, "-shots"))
    {
        const char *shots_str = getCmdOption(argv, argv + argc, "-shots");
        total_shots = stoi(shots_str);
    }

    if (cmdOptionExists(argv, argv + argc, "-basis"))
    {
        run_with_basis = true;
    }
    if (cmdOptionExists(argv, argv + argc, "-app"))
    {
        const char *filename = getCmdOption(argv, argv + argc, "-app");

        qasm_parser parser(filename);

        Simulation sim(parser.num_qubits());
        map<string, IdxType> *counts = parser.execute(sim, total_shots);

        print_counts(counts, total_shots);
    }

    if (cmdOptionExists(argv, argv + argc, "-t"))
    {
        int benchmark_index = stoi(getCmdOption(argv, argv + argc, "-t"));

        printf("Fidelity between SVSim and Qiskit Execution: %.4f\n", run_brnchmark(benchmark_index, total_shots, run_with_basis));
    }

    if (cmdOptionExists(argv, argv + argc, "-a"))
    {
        bool passed = true;
        for (int benchmark_index = 0; benchmark_index < 36; benchmark_index++)
        {
            ValType fidelity = run_brnchmark(benchmark_index, total_shots, run_with_basis);
            if (fidelity < pass_threshold)
            {
                cout << "Benchmark " << benchmark_index << " fidelity: " << fidelity << " Failed!" << endl;
                passed = false;
            }
        }
        if (passed)
            cout << "TEST PASSED!" << endl;
    }
    return 0;
}

ValType run_brnchmark(IdxType index, IdxType total_shots, bool is_basis)
{
    stringstream ss_file, ss_result;

    if (is_basis)
    {
        ss_file << "../../../data/benchmarks_basis/circuits/" << index << ".qasm";
        ss_result << "../../../data/benchmarks_basis/results/" << index << "_result.txt";
    }
    else
    {
        ss_file << "../../../data/benchmarks/circuits/" << index << ".qasm";
        ss_result << "../../../data/benchmarks/results/" << index << "_result.txt";
    }

    ifstream resultFile(ss_result.str().c_str());
    if (!resultFile)
    {
        printf("%s\n", "Could not open result file\n");
        return -1;
    }

    qasm_parser parser(ss_file.str().c_str());

    IdxType num_qubits = parser.num_qubits();
    Simulation sim(num_qubits);
    map<string, IdxType> *svsim_counts = parser.execute(sim, total_shots);
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
    for (auto const &[key, val] : *svsim_counts)
    {
        ValType out_prob = val / (ValType)total_shots;
        svsim_norm += out_prob * out_prob;
        svsim_probs.insert({key, out_prob});
    }
    ref_norm = sqrt(ref_norm);
    svsim_norm = sqrt(svsim_norm);

    ValType fidelity = 0;
    for (auto const &[key, val] : svsim_probs)
        if (ref_probs.count(key))
            fidelity += (val / svsim_norm) * (ref_probs.at(key) / ref_norm);

    return fidelity;
}
