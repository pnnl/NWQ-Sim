#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <memory>
#include <csignal>
#include "src/qasm_parser.hpp"
#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/nwq_util.hpp"

#include <tamm/tamm.hpp>

using namespace NWQSim;
using ValType = double;
namespace fs = std::filesystem;

static constexpr int MIN_QB       = 9;
static constexpr int MAX_QB       = 10;
static constexpr int CIRCUITS     = 5;
static constexpr int SHOTS        = 10024;
static constexpr double THRESHOLD = 2.0 / 3.0;

void segfault_handler(int) {
    std::exit(1);
}

int main(int argc, char **argv) {
    tamm::initialize(argc, argv);
    std::signal(SIGSEGV, segfault_handler);

    const std::string qasm_dir = "qv_qasm";
    const std::string csv_path = "bond_vs_qv.csv";

    std::ofstream csv(csv_path);
    if (!csv) {
        std::cerr << "[ERROR] Unable to open " << csv_path << "\n";
        return 1;
    }
    csv << "bond_dimension,log2(QV),QV\n";

    for (int bd = 15; bd <= 15; ++bd) {
        std::cerr << "[INFO] Starting bond_dimension=" << bd << "\n";
        int max_pass = 0;

        for (int k = MIN_QB; k <= MAX_QB; ++k) {
            std::cerr << "[INFO] Testing k=" << k << " qubits\n";
            double sum_hop = 0.0;
            size_t dim = size_t(1) << k;

            for (int c = 0; c < CIRCUITS; ++c) {
                std::cerr << "[DEBUG] bd=" << bd
                          << " k=" << k
                          << " circuit=" << c
                          << " -- loading QASM\n";

                std::ostringstream fname;
                fname << qasm_dir << "/QV_" << k << "q_"
                      << std::setfill('0') << std::setw(3) << c
                      << ".qasm";
                fs::path qasm_path = fname.str();
                if (!fs::exists(qasm_path)) {
                    std::cerr << "[ERROR] Missing file: " << qasm_path << "\n";
                    return 1;
                }

                // Ideal statevector sampling to estimate p_sv
                std::cerr << "[DEBUG] Executing SV simulation\n";
                qasm_parser parser_sv;
                parser_sv.load_qasm_file(qasm_path.string());
                auto state_sv = BackendManager::create_state("CPU", k, "SV");
                auto counts_sv = parser_sv.execute(state_sv, "", "", SHOTS);
                std::cerr << "[DEBUG] counts_sv size = " << counts_sv->size() << "\n";

                std::vector<ValType> p_sv(dim, 0.0);
                for (auto &ent : *counts_sv) {
                    int idx = std::stoi(ent.first, nullptr, 2);
                    p_sv[idx] = ent.second / ValType(SHOTS);
                }

                auto p_sorted = p_sv;
                std::sort(p_sorted.begin(), p_sorted.end());
                ValType median = (dim & 1)
                    ? p_sorted[dim/2]
                    : 0.5 * (p_sorted[dim/2 - 1] + p_sorted[dim/2]);
                std::cerr << "[DEBUG] median = " << median << "\n";

                std::vector<bool> is_heavy(dim);
                for (size_t i = 0; i < dim; ++i) {
                    is_heavy[i] = (p_sv[i] > median);
                }

                // Tensor-network sampling
                std::cerr << "[DEBUG] Executing TN simulation with bd=" << bd << "\n";
                qasm_parser parser_tn;
                parser_tn.load_qasm_file(qasm_path.string());
                auto state_tn = BackendManager::create_state("CPU_TAMM", k, "TN", bd);
                auto counts_tn = parser_tn.execute(state_tn, "", "", SHOTS);
                std::cerr << "[DEBUG] counts_tn size = " << counts_tn->size() << "\n";

                int heavy_hits = 0;
                for (auto &ent : *counts_tn) {
                    int idx = std::stoi(ent.first, nullptr, 2);
                    if (is_heavy[idx]) {
                        heavy_hits += ent.second;
                    }
                }

                ValType hop = heavy_hits / ValType(SHOTS);
                sum_hop += hop;
                std::cerr << "[DEBUG] circuit=" << c
                          << " hop=" << hop << "\n";
            }

            ValType avg_hop = sum_hop / CIRCUITS;
            std::cerr << "[INFO] k=" << k
                      << " average_hop=" << avg_hop << "\n";
            if (avg_hop > THRESHOLD) {
                max_pass = k;
                std::cerr << "[INFO] hop > threshold, updating max_pass=" << max_pass << "\n";
            } else {
                std::cerr << "[INFO] hop <= threshold, breaking k-loop\n";
                break;
            }
        }

        int log2_qv = max_pass;
        int qv      = (log2_qv > 0 ? (1 << log2_qv) : 0);
        csv << bd << "," << log2_qv << "," << qv << "\n";
        std::cerr << "[RESULT] bond_dimension=" << bd
                  << " log2(QV)=" << log2_qv
                  << " QV=" << qv << "\n";
    }
    tamm::finalize();

    printf("Finished!");
    return 0;
}


