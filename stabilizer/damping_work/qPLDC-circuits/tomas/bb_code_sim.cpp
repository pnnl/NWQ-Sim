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

// ...existing code...
    sim_timer.stop_timer();
    double local_time = sim_timer.measure();
    double max_time = 0.0;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        sim_time = max_time;
    }

    // --- Parallel I/O Strategy using MPI-IO ---

    // 1. Convert local results (vector of vectors) to a single string buffer for writing.
    std::stringstream ss;
    for (const auto& measurements : local_results) {
        for (size_t k = 0; k < measurements.size(); ++k) {
            ss << measurements[k] << (k == measurements.size() - 1 ? "" : " ");
        }
        ss << "\n";
    }
    std::string local_data_str = ss.str();
    long long local_size = local_data_str.length();

    // 2. Use MPI_Exscan to calculate the file offset for each rank.
    // Each rank's offset is the sum of the data sizes of all preceding ranks.
    long long offset = 0;
    MPI_Exscan(&local_size, &offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    // 3. Open the file in parallel.
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, m_filepath.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    // 4. Each rank writes its data block to its calculated offset in the file.
    // This is a collective operation; all ranks participate.
    MPI_File_write_at_all(fh, offset, local_data_str.c_str(), static_cast<int>(local_size), MPI_CHAR, MPI_STATUS_IGNORE);

    // 5. Close the file.
    MPI_File_close(&fh);
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