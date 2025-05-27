#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <csignal>
#include "src/qasm_parser.hpp"
#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/nwq_util.hpp"
#include <tamm/tamm.hpp>

using namespace NWQSim;
using ValType = double;
namespace fs = std::filesystem;

void segfault_handler(int) {
    std::cerr << "[FATAL] Segmentation fault\n";
    std::exit(1);
}

int main(int argc, char** argv) {
    std::signal(SIGSEGV, segfault_handler);
    tamm::initialize(argc, argv);

    // Defaults
    int n_qubits = 10;
    int bd_min   = 1;
    int bd_max   = 20;
    int CIRCUITS = 5;
    int SHOTS    = 8192;
    std::string qasm_dir = "qv_qasm";

    // Positional overrides:
    //   argv[1] = n_qubits
    //   argv[2] = bd_min
    //   argv[3] = bd_max
    //   argv[4] = CIRCUITS
    //   argv[5] = SHOTS
    //   argv[6] = qasm_dir
    if (argc > 1) n_qubits = std::stoi(argv[1]);
    if (argc > 2) bd_min   = std::stoi(argv[2]);
    if (argc > 3) bd_max   = std::stoi(argv[3]);
    if (argc > 4) CIRCUITS = std::stoi(argv[4]);
    if (argc > 5) SHOTS    = std::stoi(argv[5]);
    if (argc > 6) qasm_dir = argv[6];

    if (bd_min > bd_max) {
        std::cerr << "[ERROR] bd_min (" << bd_min
                  << ") > bd_max (" << bd_max << ")\n";
        return 1;
    }

    std::cerr << "[INFO] n_qubits="  << n_qubits
              << "  bd_min="   << bd_min
              << "  bd_max="   << bd_max
              << "  circuits=" << CIRCUITS
              << "  shots="    << SHOTS
              << "  qasm_dir=" << qasm_dir << "\n";

    const std::string csv_path = "bond_vs_qv.csv";
    std::ofstream csv(csv_path);
    if (!csv) {
        std::cerr << "[ERROR] Cannot open " << csv_path << "\n";
        return 1;
    }
    csv << "bond_dimension,mean_hop,std_error\n";

    size_t dim = size_t(1) << n_qubits;
    std::vector<ValType> p_sv(dim);
    std::vector<bool> is_heavy(dim);

    // Loop over bond dimensions
    for (int bd = bd_min; bd <= bd_max; ++bd) {
        std::cerr << "[INFO] Scanning bond_dim=" << bd << "\n";
        std::vector<ValType> hops;
        hops.reserve(CIRCUITS);

        // Loop over circuits
        for (int c = 0; c < CIRCUITS; ++c) {
            // Build filename QV_<n_qubits>q_<ccc>.qasm
            std::ostringstream fn;
            fn << qasm_dir << "/QV_" << n_qubits << "q_"
               << std::setw(3) << std::setfill('0') << c
               << ".qasm";
            fs::path path = fn.str();
            if (!fs::exists(path)) {
                std::cerr << "[ERROR] Missing file: " << path << "\n";
                return 1;
            }

            // --- Exact SV sample to build p_sv ---
            qasm_parser p_sv_parser;
            p_sv_parser.load_qasm_file(path.string());
            auto state_sv = BackendManager::create_state("CPU", n_qubits, "SV");
            auto counts_sv = p_sv_parser.execute(state_sv, "", "", SHOTS);

            std::fill(p_sv.begin(), p_sv.end(), 0.0);
            for (auto &ent : *counts_sv) {
                int idx = std::stoi(ent.first, nullptr, 2);
                p_sv[idx] = ent.second / ValType(SHOTS);
            }

            // Compute median of p_sv
            std::vector<ValType> sorted = p_sv;
            std::sort(sorted.begin(), sorted.end());
            ValType median = (dim & 1)
                           ? sorted[dim/2]
                           : 0.5*(sorted[dim/2 - 1] + sorted[dim/2]);

            // Mark heavy outputs
            for (size_t i = 0; i < dim; ++i) {
                is_heavy[i] = (p_sv[i] > median);
            }

            // --- TN sample with bond=bd ---
            qasm_parser p_tn_parser;
            p_tn_parser.load_qasm_file(path.string());
            auto state_tn = BackendManager::create_state("CPU_TAMM",
                                                         n_qubits,
                                                         "TN",
                                                         bd);
            auto counts_tn = p_tn_parser.execute(state_tn, "", "", SHOTS);

            int heavy_hits = 0;
            for (auto &ent : *counts_tn) {
                int idx = std::stoi(ent.first, nullptr, 2);
                if (is_heavy[idx]) heavy_hits += ent.second;
            }
            ValType hop = heavy_hits / ValType(SHOTS);
            hops.push_back(hop);
            std::cerr << "  circuit=" << c << "  hop=" << hop << "\n";
        }

        // Compute mean & standard error
        ValType sum = 0, sum_sq = 0;
        for (auto h : hops) {
            sum    += h;
            sum_sq += h*h;
        }
        ValType mean = sum / CIRCUITS;
        ValType var  = (CIRCUITS>1)
                     ? (sum_sq - CIRCUITS*mean*mean)/(CIRCUITS-1)
                     : 0.0;
        ValType stderr = std::sqrt(var / CIRCUITS);

        csv << bd << "," << mean << "," << stderr << "\n";
        std::cerr << "[RESULT] bd=" << bd
                  << "  mean_hop=" << mean
                  << "  stderr="  << stderr << "\n";
    }

    csv.close();
    tamm::finalize();
    return 0;
}

