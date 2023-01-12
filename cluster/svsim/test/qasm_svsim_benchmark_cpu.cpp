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
#include "../qasm/qasm_parser.hpp"
#include "../qasm/parser_util.hpp"
/**************************************************************************/

using namespace NWQSim;
double pass_threshold = 0.998;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int i_gpu = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &i_gpu); 
    IdxType total_shots = 16384;

    if (cmdOptionExists(argv, argv + argc, "-shots"))
    {
        const char *shots_str = getCmdOption(argv, argv + argc, "-shots");
        total_shots = stoi(shots_str);
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
        stringstream ss;
        ss << "../../../data/benchmarks/results/" << benchmark_index << "_result.txt";
        ifstream resultFile(ss.str().c_str());
        if (!resultFile)
        {
            printf("%s\n", "Could not open result file\n");
            return -1;
        }
        ss.str(string());
        ss << "../../../data/benchmarks/circuits/" << benchmark_index << ".qasm";

        qasm_parser parser(ss.str().c_str());
        int num_qubits = parser.num_qubits();
        Simulation sim(num_qubits);
        map<string, IdxType> *svsim_counts = parser.execute(sim, total_shots);
        map<string, ValType> ref_probs;
        map<string, ValType> svsim_probs;
        string line;
        double ref_norm = 0, svsim_norm = 0;
        while (getline(resultFile, line))
        {
            vector<string> ref_result = split(line, ' ');
            string outcome = ref_result[0];
            double ref_prob = stod(ref_result[1]);
            ref_probs.insert({outcome, ref_prob});
            ref_norm += ref_prob * ref_prob;
        }
        for (auto const &[key, val] : *svsim_counts)
        {
            double out_prob = val / (double)total_shots;
            svsim_norm += out_prob * out_prob;
            svsim_probs.insert({key, out_prob});
        }
        ref_norm = sqrt(ref_norm);
        svsim_norm = sqrt(svsim_norm);

        double fidelity = 0;
        for (auto const &[key, val] : svsim_probs)
            if (ref_probs.count(key))
                fidelity += (val / svsim_norm) * (ref_probs.at(key) / ref_norm);
        if (i_gpu == 0)
        {
            print_counts(svsim_counts, total_shots);
            printf("Fidelity between SVSim and Qiskit Execution: %.4f\n", fidelity);
        }
    }

    if (cmdOptionExists(argv, argv + argc, "-a"))
    {
        bool passed = true;
        for (int benchmark_index = 16; benchmark_index < 20; benchmark_index++)
        {
            int total_shots = 16384;
            stringstream ss;
            ss << "../../../data/benchmarks/results/" << benchmark_index << "_result.txt";
            ifstream resultFile(ss.str().c_str());
            if (!resultFile)
            {
                printf("%s \n", "Could not open result file\n");
                return -1;
            }
            ss.str(string());
            ss << "../../../data/benchmarks/circuits/" << benchmark_index << ".qasm";
            qasm_parser parser(ss.str().c_str());
            int num_qubits = parser.num_qubits();
            Simulation sim(num_qubits);
            map<string, IdxType> *svsim_counts = parser.execute(sim, total_shots);
            map<string, ValType> ref_probs;
            map<string, ValType> svsim_probs;
            string line;
            double ref_norm = 0, svsim_norm = 0;
            while (getline(resultFile, line))
            {
                vector<string> ref_result = split(line, ' ');

                string outcome = ref_result[0];
                double ref_prob = stod(ref_result[1]);

                ref_probs.insert({outcome, ref_prob});

                ref_norm += ref_prob * ref_prob;
            }

            for (auto const &[key, val] : *svsim_counts)
            {
                double out_prob = val / (double)total_shots;
                svsim_norm += out_prob * out_prob;

                svsim_probs.insert({key, out_prob});
            }

            ref_norm = sqrt(ref_norm);
            svsim_norm = sqrt(svsim_norm);

            double fidelity = 0;
            for (auto const &[key, val] : svsim_probs)
                if (ref_probs.count(key))
                    fidelity += (val / svsim_norm) * (ref_probs.at(key) / ref_norm);

            if (fidelity < pass_threshold && i_gpu == 0)
            {
                cout << "Benchmark " << benchmark_index << " Failed!" << endl;
                passed = false;
            }
        }
        if (passed && i_gpu == 0) cout << "TEST PASSED!" << endl;
    }

    MPI_Finalize();
    return 0;
}
