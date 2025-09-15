// #define EIGEN
// #define OMP_ENABLED

#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <ranges>
#include <fstream>
#include <cstdint>
#include <cstring>
#include <algorithm>

#include <mpi.h>

#include "../../include/backendManager.hpp"
#include "../../include/state.hpp"
#include "../../include/circuit.hpp"
#include "../../include/nwq_util.hpp"

#include "../../include/stabsim/src/qasm_extraction.hpp"
namespace NWQSim{
void sim_batch(std::shared_ptr<NWQSim::Circuit> circuit, IdxType shots, std::vector<std::vector<int32_t>>& all_results, int n, double &sim_time)
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
    int mlen_local = local_results.empty() ? 0 : static_cast<int>(local_results[0].size());
    int mlen_min = 0, mlen_max = 0;
    MPI_Allreduce(&mlen_local, &mlen_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&mlen_local, &mlen_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (mlen_min != mlen_max) {
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

// New scalable streaming simulation with MPI-IO (supports billions of shots)
// If output filename ends with .bin a compact bit-packed binary format is used:
//   Header (64 bytes): magic "NWQBIN01", n_qubits (u64), shots (u64), mlen (u64), format=1
//   Followed by shots * ceil(mlen/8) bytes (LSB-first within each byte)
// Otherwise a fixed-width text format identical to the legacy output is used:
//   Each line: mlen bits separated by single spaces + newline (length = 2*mlen bytes)
void sim_stream_write(std::shared_ptr<NWQSim::Circuit> circuit,
                      uint64_t shots,
                      int n_qubits,
                      const std::string &outfile,
                      double &sim_time_max)
{
    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    uint64_t base = shots / world_size;
    uint64_t rem  = shots % world_size;
    uint64_t local_count = base + (static_cast<uint64_t>(world_rank) < rem ? 1 : 0);
    uint64_t local_start = static_cast<uint64_t>(world_rank) * base + std::min<uint64_t>(world_rank, rem);

    bool binary_bitpacked = false;   // format=1
    bool binary_rawbytes  = false;   // format=2 (NEW: 1 byte per measurement bit)
    bool text_legacy      = false;   // no header

    // Determine mode from extension first
    if (outfile.size() >= 4) {
        std::string ext = outfile.substr(outfile.size()-4);
        if (ext == ".bin") binary_bitpacked = true; // default .bin -> bitpacked
    }
    // Environment override
    const char *env_fmt = std::getenv("NWQSIM_OUTPUT_FORMAT");
    if (env_fmt) {
        std::string fmt(env_fmt);
        for (auto &c: fmt) c = (char)std::tolower(c);
        if (fmt == "bin" || fmt == "bit" || fmt == "bitpacked") {
            binary_bitpacked = true; binary_rawbytes = false; text_legacy = false;
        } else if (fmt == "raw" || fmt == "rawbytes") {
            binary_bitpacked = false; binary_rawbytes = true; text_legacy = false;
        } else if (fmt == "txt" || fmt == "text") {
            binary_bitpacked = false; binary_rawbytes = false; text_legacy = true;
        }
    }
    if (!binary_bitpacked && !binary_rawbytes && !text_legacy) {
        // default if nothing chosen
        binary_bitpacked = true; // keep previous behavior by default
    }

    if (world_rank == 0) {
        if (text_legacy) std::cerr << "[info] Output format = TEXT legacy (large).\n";
        else if (binary_bitpacked) std::cerr << "[info] Output format = BIN bit-packed (format=1).\n";
        else if (binary_rawbytes) std::cerr << "[info] Output format = BIN RAW bytes (format=2, debugging friendly).\n";
    }

    // Open MPI file collectively
    MPI_File fh; MPI_File_open(MPI_COMM_WORLD, outfile.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    auto sim = BackendManager::create_state("CPU", n_qubits, "stab");
    double tmp_timer = 0.0;
    cpu_timer wall_timer; wall_timer.start_timer();

    size_t mlen = 0; // measurement length
    size_t stride = 0; // per-shot bytes
    const size_t header_size = (text_legacy ? 0 : 64);
    std::vector<int32_t> meas_first;

    if (local_count > 0) {
        sim->set_seed(local_start);
        sim->sim(circuit, tmp_timer);
        meas_first = sim->get_measurement_results();
        mlen = meas_first.size();
    }
    uint64_t mlen_local = mlen, mlen_min = 0, mlen_max = 0;
    MPI_Allreduce(&mlen_local, &mlen_min, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&mlen_local, &mlen_max, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
    if (mlen_min != mlen_max) {
        if (world_rank == 0) std::cerr << "[error] Measurement length mismatch across ranks." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
    mlen = mlen_max;

    uint64_t fmt_flag = 0; // 0 unused, 1 = bitpacked, 2 = rawbytes
    if (binary_bitpacked) {
        stride = (mlen + 7) / 8; fmt_flag = 1;
    } else if (binary_rawbytes) {
        stride = mlen; fmt_flag = 2; // 1 byte per bit
    } else { // text
        stride = 2 * mlen; // digits + spaces
    }

    if (!text_legacy && world_rank == 0) {
        std::vector<unsigned char> header(header_size, 0);
        const char magic[8] = {'N','W','Q','B','I','N','0','1'}; // unchanged
        std::memcpy(header.data(), magic, 8);
        std::memcpy(header.data()+8,  &n_qubits, sizeof(uint64_t));
        std::memcpy(header.data()+16, &shots,    sizeof(uint64_t));
        std::memcpy(header.data()+24, &mlen,     sizeof(uint64_t));
        std::memcpy(header.data()+32, &fmt_flag, sizeof(uint64_t));
        MPI_File_write_at(fh, 0, header.data(), (int)header.size(), MPI_BYTE, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<unsigned char> bin_buf((binary_bitpacked||binary_rawbytes)? stride : 0);
    std::vector<char> txt_buf(text_legacy ? stride : 0, '\n');

    auto write_measurement = [&](uint64_t g_shot, const std::vector<int32_t>& mv){
        MPI_Offset base_off = (MPI_Offset)header_size + (MPI_Offset)g_shot * (MPI_Offset)stride;
        if (binary_bitpacked) {
            std::fill(bin_buf.begin(), bin_buf.end(), 0);
            for (size_t k = 0; k < mv.size(); ++k) if (mv[k]) bin_buf[k >> 3] |= (unsigned char)(1u << (k & 7));
            MPI_File_write_at(fh, base_off, bin_buf.data(), (int)bin_buf.size(), MPI_BYTE, MPI_STATUS_IGNORE);
        } else if (binary_rawbytes) {
            for (size_t k = 0; k < mv.size(); ++k) bin_buf[k] = (unsigned char)(mv[k] ? 1 : 0);
            MPI_File_write_at(fh, base_off, bin_buf.data(), (int)mv.size(), MPI_BYTE, MPI_STATUS_IGNORE);
        } else { // text
            for (size_t k = 0; k < mlen; ++k) { txt_buf[2*k] = mv[k] ? '1' : '0'; if (k < mlen - 1) txt_buf[2*k+1] = ' '; }
            txt_buf[stride-1] = '\n';
            MPI_File_write_at(fh, base_off, txt_buf.data(), (int)txt_buf.size(), MPI_CHAR, MPI_STATUS_IGNORE);
        }
    };

    if (local_count > 0) { write_measurement(local_start, meas_first); sim->reset_state(); }
    for (uint64_t j = (local_count > 0 ? 1 : 0); j < local_count; ++j) {
        uint64_t g = local_start + j;
        sim->set_seed(g); sim->sim(circuit, tmp_timer);
        const auto &mv = sim->get_measurement_results();
        if (mv.size() != mlen) { std::cerr << "[error] Measurement length changed." << std::endl; MPI_Abort(MPI_COMM_WORLD, 3); }
        write_measurement(g, mv); sim->reset_state();
    }

    wall_timer.stop_timer();
    double local_time_ms = wall_timer.measure();
    double max_time_ms = 0.0;
    MPI_Reduce(&local_time_ms, &max_time_ms, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        sim_time_max = max_time_ms;
        std::cout << "Total C++ simulation time: " << (sim_time_max / 1000.0) << "s" << std::endl;
        long double est_bytes = (long double)shots * (long double)stride + header_size;
        long double gib = est_bytes / (1024.0L*1024.0L*1024.0L);
        std::cerr << "[info] Approx output size: " << std::fixed << std::setprecision(2) << gib << " GiB" << (text_legacy?" (text)":" (binary)") << std::endl;
    }

    MPI_File_close(&fh);
}
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int world_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (argc < 5) {
        if (world_rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <qubits> <shots> <qasm filepath> <measurement filepath[.txt|.bin]>\n";
        }
        MPI_Finalize();
        return 1;
    }

    int n_qubits = std::atoi(argv[1]);
    uint64_t shots = std::strtoull(argv[2], nullptr, 10);
    std::string qasmfile = argv[3];
    std::string m_filepath = argv[4];

    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
    appendQASMToCircuit(circuit, qasmfile, n_qubits);

    double total_time_ms = 0.0;
    // Scalable streaming write (no massive gather)
    NWQSim::sim_stream_write(circuit, shots, n_qubits, m_filepath, total_time_ms);

    MPI_Finalize();
    return 0;
}