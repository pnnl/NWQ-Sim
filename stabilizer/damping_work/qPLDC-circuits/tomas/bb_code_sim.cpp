// #define EIGEN
// #define OMP_ENABLED

#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <ranges>
#include <fstream>
#include <limits>

#include <mpi.h>

#include "../../../../include/backendManager.hpp"
#include "../../../../include/state.hpp"
#include "../../../../include/circuit.hpp"
#include "../../../../include/nwq_util.hpp"

#include "../../../../include/stabsim/src/qasm_extraction.hpp"
namespace NWQSim{
void sim_batch(std::shared_ptr<NWQSim::Circuit> circuit, IdxType shots, int n, double &sim_time, const std::string& m_filepath)
{
    // MPI rank/size
    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Partition shots across ranks (contiguous blocks)
    IdxType base = shots / world_size;
    IdxType rem  = shots % world_size;
    IdxType local_count = base + (world_rank < rem ? 1 : 0);
    IdxType local_start = world_rank * base + std::min<IdxType>(world_rank, rem);

    // Local simulation
    cpu_timer sim_timer; sim_timer.start_timer();
    double tmp_timer = 0.0;
    auto sim = BackendManager::create_state("CPU", n, "stab");

    // Each rank simulates and stores results in memory
    std::vector<std::vector<int32_t>> local_results;
    local_results.reserve(static_cast<size_t>(local_count));

    for (IdxType j = 0; j < local_count; ++j) {
        IdxType g = local_start + j; // global shot index
        sim->set_seed(static_cast<IdxType>(g));
        sim->sim(circuit, tmp_timer);
        local_results.push_back(sim->get_measurement_results());
        sim->reset_state();
    }

    sim_timer.stop_timer();
    double local_time = sim_timer.measure();
    double max_time = 0.0;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        sim_time = max_time;
    }

    // Ranks write to a single file sequentially to avoid race conditions
    for (int i = 0; i < world_size; ++i) {
        if (world_rank == i) {
            // Rank 0 opens in write mode (truncates), others append
            std::ofstream out_file(m_filepath, (i == 0) ? std::ios::out : std::ios::app);
            if (!out_file.is_open()) {
                std::cerr << "Rank " << world_rank << " error opening " << m_filepath << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            for (const auto& measurements : local_results) {
                for (size_t k = 0; k < measurements.size(); ++k) {
                    out_file << measurements[k] << (k == measurements.size() - 1 ? "" : " ");
                }
                out_file << "\n";
            }
            out_file.close();
        }
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize before next rank writes
    }
}
}

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc < 5) {
        if (world_rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <qubits> <shots> <qasm filepath> <measurement filepath>\n";
        }
        MPI_Finalize();
        return 1;
    }

    int n_qubits = std::atoi(argv[1]);
    int iters = std::atoi(argv[2]);

    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

    std::string qasmfile = argv[3];
    appendQASMToCircuit(circuit, qasmfile, n_qubits);
    std::string m_filepath = argv[4];

    double total_time = 0.0;

    // Distribute simulation across MPI processes
    NWQSim::sim_batch(circuit, iters, n_qubits, total_time, m_filepath);

    // Rank 0 reports time
    if (world_rank == 0) {
        std::cout << "Total C++ simulation time: " << total_time / 1000 << "s" << std::endl;
    }

    MPI_Finalize();
    return 0;
}