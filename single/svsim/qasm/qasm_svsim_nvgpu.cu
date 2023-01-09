#include <time.h>
#include <iostream>
#include <fstream>
// stats calculation
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>

#include "qasm_parser.hpp"

using namespace NWQSim;

int main(int argc, char **argv)
{
    const char *filepath;
    IdxType n_shots = 1024;

    // get file name
    if (cmdOptionExists(argv, argv + argc, "-q"))
    {
        filepath = getCmdOption(argv, argv + argc, "-q");
        qasm_parser parser(filepath);

        IdxType n_qubits = parser.num_qubits();
        Simulation sim(n_qubits);
        map<string, IdxType>* res = parser.execute(sim, n_shots);
        //if (sim.i_gpu==0) printf("\n -------- qasm_svsim ---------\n");
        if (sim.i_gpu==0) print_counts(res, n_shots);
        fflush(stdout);
        //sim.print_res_sv();
    }
    else
    {
        cerr << "\n ------------- SVSim for qasm ---------------\n" << "  Use ./qasm_svsim -q input.qasm" << "\n --------------------------------------------\n";
    }
    
    return 0;
}
