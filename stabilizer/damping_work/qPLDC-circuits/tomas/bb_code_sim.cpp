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
void sim_batch(std::shared_ptr<NWQSim::Circuit> circuit, IdxType shots, std::vector<std::vector<int32_t>>& all_results, int n, double &sim_time)
{
    // std::cerr << "In sim batch!" << std::endl;
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

    // Determine measurement length; ensure consistent across ranks
    int mlen_local = local_results.empty() ? -1 : static_cast<int>(local_results[0].size());
    int mlen_max = 0;
    MPI_Allreduce(&mlen_local, &mlen_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    int big = std::numeric_limits<int>::max();
    int mlen_min_local = local_results.empty() ? big : mlen_local;
    int mlen_min = 0;
    MPI_Allreduce(&mlen_min_local, &mlen_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (mlen_min == big) mlen_min = 0; // all empty (shouldn't happen if shots>0)

    if (mlen_min != mlen_max && mlen_min != 0) {
        if (world_rank == 0) {
            std::cerr << "Error: measurement length mismatch across MPI ranks." << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int mlen = mlen_max;

    // Flatten local results for gather
    std::vector<int32_t> local_flat;
    local_flat.resize(static_cast<size_t>(local_count) * mlen);
    for (IdxType j = 0; j < local_count; ++j) {
        const auto &v = local_results[static_cast<size_t>(j)];
        for (int k = 0; k < mlen; ++k) {
            local_flat[static_cast<size_t>(j) * mlen + k] = v[static_cast<size_t>(k)];
        }
    }

    // Gather counts
    int local_count_i = static_cast<int>(local_count);
    std::vector<int> recv_counts;
    if (world_rank == 0) recv_counts.resize(world_size, 0);
    MPI_Gather(&local_count_i, 1, MPI_INT,
               world_rank == 0 ? recv_counts.data() : nullptr, 1, MPI_INT,
               0, MPI_COMM_WORLD);

    // Prepare Gatherv parameters on root
    std::vector<int> displs_elems; // element displacements (in ints)
    std::vector<int> recv_counts_elems;
    std::vector<int32_t> recv_flat; // flattened buffer on root
    if (world_rank == 0) {
        displs_elems.resize(world_size, 0);
        recv_counts_elems.resize(world_size, 0);
        int offset_elems = 0;
        for (int r = 0; r < world_size; ++r) {
            displs_elems[r] = offset_elems;
            recv_counts_elems[r] = recv_counts[r] * mlen;
            offset_elems += recv_counts_elems[r];
        }
        recv_flat.resize(static_cast<size_t>(offset_elems));
    }

    // Gather all results to root
    int send_elems = local_count_i * mlen;
    MPI_Gatherv(local_flat.data(), send_elems, MPI_INT,
                world_rank == 0 ? recv_flat.data() : nullptr,
                world_rank == 0 ? recv_counts_elems.data() : nullptr,
                world_rank == 0 ? displs_elems.data() : nullptr,
                MPI_INT, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        // Reconstruct all_results in global shot order
        all_results.clear();
        all_results.resize(shots);
        for (IdxType g = 0; g < shots; ++g) {
            all_results[static_cast<size_t>(g)].assign(
                recv_flat.begin() + static_cast<IdxType>(g) * mlen,
                recv_flat.begin() + static_cast<IdxType>(g + 1) * mlen
            );
        }
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

    std::ofstream out_file;
    if (world_rank == 0) {
        out_file.open(m_filepath);
        if (!out_file.is_open()) {
            std::cerr << "Error opening measurement_results.txt" << std::endl;
            MPI_Finalize();
            return 1;
        }
    }

    std::vector<std::vector<int32_t>> all_measurements; // only used on rank 0
    double total_time = 0.0;

    // Distribute simulation across MPI processes
    all_measurements.resize(world_rank == 0 ? iters : 0);
    // std::cerr << "Starting sim batch!!" << std::endl;
    NWQSim::sim_batch(circuit, iters, all_measurements, n_qubits, total_time);

    // Rank 0 writes results
    if (world_rank == 0) {
        for (const auto& measurements : all_measurements) {
            for (size_t j = 0; j < measurements.size(); ++j) {
                out_file << measurements[j] << (j == measurements.size() - 1 ? "" : " ");
            }
            out_file << "\n";
        }
        out_file.close();
        std::cout << "Total C++ simulation time: " << total_time / 1000 << "s" << std::endl;
    }

    MPI_Finalize();
    return 0;
}