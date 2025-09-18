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

#include "include/backendManager.hpp"
#include "include/state.hpp"
#include "include/circuit.hpp"
#include "include/nwq_util.hpp"

#include "include/stabsim/src/qasm_extraction.hpp"
// ...existing includes...
#include <cstdint>
#include <cstring>
#include <memory>
#include <cstdlib>
// ...existing code...

using IdxType = long long int;
namespace NWQSim{
void sim_batch(std::shared_ptr<NWQSim::Circuit> circuit,
               IdxType shots,
               const std::string& out_path,
               int n,
               double &sim_time)
{
    // MPI rank/size
    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Partition shots across ranks (contiguous blocks)
    IdxType base = shots / world_size;
    IdxType rem  = shots % world_size;
    IdxType local_count = base + (world_rank < rem ? 1 : 0);
    IdxType local_start = static_cast<IdxType>(world_rank) * base + static_cast<IdxType>(std::min<IdxType>(world_rank, rem));

    // Simulator
    cpu_timer sim_timer; sim_timer.start_timer();
    double tmp_timer = 0.0;

    // Probe first local shot to learn measurement length
    std::vector<int32_t> first_meas;
    long long mlen_local_ll = -1LL;
    if (local_count > 0) {
        auto sim = BackendManager::create_state("CPU", n, "stab");
        IdxType g0 = local_start;
        sim->set_seed(static_cast<IdxType>(g0));
        sim->sim(circuit, tmp_timer);
        first_meas = sim->get_measurement_results();
        mlen_local_ll = static_cast<long long>(first_meas.size());
        // sim is destroyed here, ensuring clean state
    }

    // Agree on measurement length across ranks
    long long mlen_max_ll = 0;
    MPI_Allreduce(&mlen_local_ll, &mlen_max_ll, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);

    long long big_ll = std::numeric_limits<long long>::max();
    long long mlen_min_local_ll = (mlen_local_ll >= 0) ? mlen_local_ll : big_ll;
    long long mlen_min_ll = 0;
    MPI_Allreduce(&mlen_min_local_ll, &mlen_min_ll, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    if (mlen_min_ll == big_ll) mlen_min_ll = 0; // all empty

    if (mlen_min_ll != mlen_max_ll && mlen_min_ll != 0) {
        if (world_rank == 0) {
            std::cerr << "Error: measurement length mismatch across MPI ranks." << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    IdxType mlen = static_cast<IdxType>(mlen_max_ll);

    // Open output file collectively and pre-size it (binary + header)
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD,
                  const_cast<char*>(out_path.c_str()),
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);

    // Header: magic 'NWQB' (u32), version (u32), mlen (u64), shots (u64)
    const uint32_t MAGIC = 0x4E575142u; // 'NWQB'
    const uint32_t VERSION = 1u;
    const MPI_Offset header_bytes = static_cast<MPI_Offset>(sizeof(uint32_t) * 2 + sizeof(uint64_t) * 2); // 24 bytes

    const IdxType bytes_per_shot_idx = (mlen + 7) / 8;
    const MPI_Offset bytes_per_shot = static_cast<MPI_Offset>(bytes_per_shot_idx);
    const MPI_Offset total_bytes = header_bytes + static_cast<MPI_Offset>(shots) * bytes_per_shot;

    if (world_rank == 0) {
        uint8_t hdr[24];
        std::memset(hdr, 0, sizeof(hdr));
        uint64_t mlen64  = static_cast<uint64_t>(mlen);
        uint64_t shots64 = static_cast<uint64_t>(shots);
        std::memcpy(hdr + 0,  &MAGIC,   sizeof(uint32_t));
        std::memcpy(hdr + 4,  &VERSION, sizeof(uint32_t));
        std::memcpy(hdr + 8,  &mlen64,  sizeof(uint64_t));
        std::memcpy(hdr + 16, &shots64, sizeof(uint64_t));
        MPI_Status st;
        MPI_File_write_at(fh, 0, hdr, static_cast<int>(sizeof(hdr)), MPI_BYTE, &st);
        // Log header info for debugging alignment with reader
        std::cout << "NWQB header v" << VERSION
                  << ": mlen=" << mlen
                  << ", shots=" << shots
                  << ", bytes_per_shot=" << static_cast<unsigned long long>((mlen + 7) / 8)
                  << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_set_size(fh, total_bytes);
    MPI_Barrier(MPI_COMM_WORLD);

    // Pack one shot’s bits into little-endian bytes
    auto pack_shot_le = [mlen](const std::vector<int32_t>& v, uint8_t* out) {
        const size_t B = static_cast<size_t>((mlen + 7) / 8);
        std::memset(out, 0, B);
        const size_t avail = v.size();
        for (IdxType k = 0; k < mlen; ++k) {
            if (static_cast<size_t>(k) < avail && v[static_cast<size_t>(k)] != 0) {
                out[static_cast<size_t>(k >> 3)] |= static_cast<uint8_t>(1u << (k & 7));
            }
        }
    };

    // Chunked write to reduce MPI calls
    IdxType chunk_shots = 8192; // configurable: export NWQSIM_IO_CHUNK_SHOTS=16384
    if (const char* env_chunk = std::getenv("NWQSIM_IO_CHUNK_SHOTS")) {
        long long v = std::atoll(env_chunk);
        if (v > 0) chunk_shots = static_cast<IdxType>(v);
    }

    // Main simulation loop
    IdxType j0 = 0;
    if (local_count > 0) {
        auto sim = BackendManager::create_state("CPU", n, "stab");
        while (j0 < local_count) {
            IdxType blk = std::min(chunk_shots, local_count - j0);
            std::vector<uint8_t> buf(static_cast<size_t>(blk) * static_cast<size_t>(bytes_per_shot));

            for (IdxType t = 0; t < blk; ++t) {
                std::vector<int32_t> meas;
                if (j0 == 0 && t == 0) {
                    meas.swap(first_meas);
                } else {
                    IdxType g = local_start + j0 + t;
                    sim->set_seed(static_cast<IdxType>(g));
                    sim->sim(circuit, tmp_timer);
                    meas = sim->get_measurement_results();
                    sim->reset_state();
                }
                uint8_t* dst = buf.data() + static_cast<size_t>(t) * static_cast<size_t>(bytes_per_shot);
                pack_shot_le(meas, dst);
            }

            MPI_Offset off = header_bytes + (static_cast<MPI_Offset>(local_start + j0)) * bytes_per_shot;
            MPI_Status st;
            MPI_File_write_at_all(fh, off, buf.data(), static_cast<int>(buf.size()), MPI_BYTE, &st);

            j0 += blk;
        }
    }

    MPI_File_close(&fh);

    sim_timer.stop_timer();
    double local_time = sim_timer.measure();
    double max_time = 0.0;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        sim_time = max_time;
        std::cout << "Total C++ simulation time: " << sim_time << "s" << std::endl;
    }
}
} // namespace NWQSim

// ...existing main...
int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc < 5) {
        if (world_rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <n_qubits> <shots> <qasm_file> <out_path>\n";
        }
        MPI_Finalize();
        return 1;
    }

    int n_qubits = std::atoi(argv[1]);
    IdxType iters = static_cast<IdxType>(std::atoll(argv[2]));
    std::string qasmfile = argv[3];
    std::string m_filepath = argv[4];

    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
    appendQASMToCircuit(circuit, qasmfile, n_qubits);

    double total_time = 0.0;

    // Parallel simulation and direct MPI-IO writing (no root gather)
    NWQSim::sim_batch(circuit, iters, m_filepath, n_qubits, total_time);

    MPI_Finalize();
    return 0;
}