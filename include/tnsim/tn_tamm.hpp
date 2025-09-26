#pragma once

#include <mpi.h>
#include <numeric> 

#include "../state.hpp"

#include "../nwq_util.hpp"
#include "../gate.hpp"
#include "../circuit.hpp"
#include "../config.hpp"
#include "private/exp_gate_declarations_host.hpp"

#include "../circuit_pass/fusion.hpp"
#include "../private/macros.hpp"
#include "../private/sim_gate.hpp"

#include <random>
#include <vector>
#include <string>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <sstream>

#include <tamm/errors.hpp>
#include <tamm/symbol.hpp>
#include <tamm/tensor.hpp>
#include <tamm/tiled_index_space.hpp>
#include <tamm/index_space.hpp>
#include <tamm/execution_context.hpp>
#include <tamm/scheduler.hpp>
#include <tamm/tamm_io.hpp>

#include <tamm/rmm_memory_pool.hpp>
#include <tamm/gpu_streams.hpp>

#include <tamm/tamm.hpp>

#include <gsl/span>
#include <iostream>

#include <complex>
#include <map>
#include <cstring>

#include <fstream>
#include <iomanip>

#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <Eigen/Dense>

namespace NWQSim
{

    using Cplx = std::complex<ValType>;

    class TN_TAMM;

    struct LocalGateResult {
        bool is_valid = false;
        IdxType q0, q1;
        IdxType new_bond_dim;
        std::vector<Cplx> new_T0_data;
        std::vector<Cplx> new_T1_data;
        int original_rank;
    };

    struct GateUpdateMetadata {
        bool is_valid = false;
        IdxType q0, q1;
        IdxType new_bond_dim;
        int original_rank;
    };

    struct CuCtx {
        cusolverDnHandle_t solver = nullptr;
        cudaStream_t stream = nullptr;
        gesvdjInfo_t jp = nullptr;
        int lwork_jac = 0;
        cuDoubleComplex* d_work_jac = nullptr;

        CuCtx() {
            cusolverDnCreate(&solver);
            cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
            cusolverDnSetStream(solver, stream);
            cusolverDnCreateGesvdjInfo(&jp);
        }
        ~CuCtx() {
            if(d_work_jac) cudaFree(d_work_jac);
            if(jp) cusolverDnDestroyGesvdjInfo(jp);
            if(solver) cusolverDnDestroy(solver);
            if(stream) cudaStreamDestroy(stream);
        }
    };


    using Eigen::Index;
    class TN_TAMM : public QuantumState
    {
    public:
        TN_TAMM(IdxType n_qubits,
                IdxType max_bond_dim = 100,
                double sv_cutoff = 0.0,
                std::string backend = "TN_TAMM_CPU")
        : QuantumState(SimType::TN),
            n_qubits(n_qubits),
            block_size(1024),
            max_bond_dim(max_bond_dim),
            sv_cutoff(sv_cutoff),
            pg(init_pg()),
            ec(pg, tamm::DistributionKind::dense, tamm::MemoryManagerKind::ga),
            pg_local_(tamm::ProcGroup::create_self()),
            ec_local_(pg_local_, tamm::DistributionKind::dense, 
                      tamm::MemoryManagerKind::local),
            sch_local_(ec_local_)                     
        {
            i_proc = pg.rank().value();
            
            // print the tamm execution context, this can be commented out if desired
            if(ec.print()) {
                auto current_time   = std::chrono::system_clock::now();
                auto current_time_t = std::chrono::system_clock::to_time_t(current_time);
                auto cur_local_time = localtime(&current_time_t);
                std::cout << std::endl << "date: " << std::put_time(cur_local_time, "%c") << std::endl;
                std::cout << "nnodes: " << ec.nnodes() << ", ";
                std::cout << "nproc_per_node: " << ec.ppn() << ", ";
                std::cout << "nproc_total: " << ec.nnodes() * ec.ppn() << ", ";
                if(ec.has_gpu()) {
                  std::cout << "ngpus_per_node: " << ec.gpn() << ", ";
                  std::cout << "ngpus_total: " << ec.nnodes() * ec.gpn() << std::endl;
                }
                std::cout << std::endl;
                ec.print_mem_info();
                std::cout << std::endl;
            }

            if (backend == "TN_TAMM_CPU")
            {
                exec_hw = tamm::ExecutionHW::CPU;
            }
            else if(backend == "TN_TAMM_GPU")
            {
                exec_hw = tamm::ExecutionHW::GPU;
            }

            // initialize bond index spaces
            bond_tis.resize(n_qubits + 1);
            bond_dims.resize(n_qubits + 1);
            for (IdxType i = 0; i <= n_qubits; ++i)
            {
                bond_dims[i] = 1;
                tamm::IndexSpace is{ tamm::range(1) };
                bond_tis[i] = tamm::TiledIndexSpace(is, block_size);
            }
        
            // initialize physical index spaces
            phys_tis.resize(n_qubits);
            phys_dims.resize(n_qubits);
            for (IdxType i = 0; i < n_qubits; ++i)
            {
                phys_dims[i] = 2;
                tamm::IndexSpace is{ tamm::range(2) };
                phys_tis[i] = tamm::TiledIndexSpace(is, 1);
            }
        
            // allocate and initialize MPS tensors
            mps_tensors.reserve(n_qubits);
            for (IdxType i = 0; i < n_qubits; ++i)
            {
                tamm::Tensor<Cplx> T({ bond_tis[i], phys_tis[i], bond_tis[i + 1] });
                T.set_dense();
                T.allocate(&ec);
        
                for (const auto & blockid : T.loop_nest())
                {
                    const size_t block_size = T.block_size(blockid);
                    std::vector<Cplx> hostbuf(block_size);
        
                    auto dims = T.block_dims(blockid);
                    auto offsets = T.block_offsets(blockid);
        
                    for (size_t idx = 0; idx < block_size; ++idx)
                    {
                        size_t rem = idx;
                        size_t l_loc = rem % dims[0];
                        rem /= dims[0];
                        size_t p_loc = rem % dims[1];
                        rem /= dims[1];
                        size_t r_loc = rem % dims[2];
        
                        IdxType l = offsets[0] + l_loc;
                        IdxType p = offsets[1] + p_loc;
                        IdxType r = offsets[2] + r_loc;
        
                        hostbuf[idx] = (l == 0 && p == 0 && r == 0)
                                       ? Cplx{1.0, 0.0}
                                       : Cplx{0.0, 0.0};
                    }
        
                    T.put(blockid, hostbuf);
                }
        
                mps_tensors.push_back(std::move(T));
            }
        }

        ~TN_TAMM() noexcept override 
        {
            SAFE_FREE_HOST(results);
        }

        void reset_state() override {
            for(IdxType i = 0; i < n_qubits; ++i)
            {
                auto& T = mps_tensors[i];
                T.loop_nest().iterate([&](auto const& idxs)
                {
                    Cplx v = (idxs[0] == 0 && idxs[1] == 0 && idxs[2] == 0)
                             ? Cplx(1.0,0.0)
                             : Cplx(0.0,0.0);
                    T.put(idxs, gsl::span<Cplx>(&v,1));
                });
            }
        }

        void set_seed(IdxType seed) override
        {
            throw std::runtime_error("TN_TAMM does not use RNG seed, not accessible form cutensornet API\n");
        }

        void set_initial(std::string fpath, std::string format) override
        {
            std::cout << "set function was called" << std::endl;
        }

        void dump_res_state(std::string outpath) override
        {
            std::cout << "dump function was called" << std::endl;
        }

        void sim(std::shared_ptr<NWQSim::Circuit> circuit) override
        {
            IdxType original_gate_count = circuit->num_gates();
            std::vector<SVGate> gates = fuse_circuit_sv(circuit);
            IdxType fused_gate_count = gates.size();
            assert(circuit->num_qubits() == n_qubits);
        
            pg.barrier();
            simulation_kernel(gates);
            pg.barrier();
        
        }

        IdxType* get_results() override
        {
            return results;
        }

        IdxType measure(IdxType qubit) override
        {
            throw std::runtime_error("TN_TAMM::measure not implemented");
        }

        IdxType* measure_all(IdxType repetition) override 
        {
            MA_GATE(repetition);
            return results; 
        }

        ValType* get_real() const override
        {
            throw std::runtime_error("TN_TAMM::get_real not implemented");
        }

        ValType* get_imag() const override
        {
            throw std::runtime_error("TN_TAMM::get_imag not implemented");
        }

        ValType get_exp_z() override
        {
            throw std::runtime_error("TN_TAMM::get_exp_z() not implemented");
        }

        ValType get_exp_z(const std::vector<size_t>& in_bits) override
        {
            throw std::runtime_error("TN_TAMM::get_exp_z(bits) not implemented");
        }

        void print_res_state() override
        {
            throw std::runtime_error("TN_TAMM::print_res_state not implemented");
        }

        static SVGate make_swap_sv(int a, int b)
        {
            SVGate s(OP::C2, b, a);
        
            static const ValType real[16] = {
                1,0,0,0,
                0,0,1,0,
                0,1,0,0,
                0,0,0,1
            };
            static const ValType imag[16] = {0};
            memcpy(s.gm_real, real, 16 * sizeof(ValType));
            memcpy(s.gm_imag, imag, 16 * sizeof(ValType));
            return s;
        }
        
        static SVGate make_local_c2_sv(const SVGate& g, int left, int right)
        {
            SVGate t(g);
            t.ctrl = left;
            t.qubit = right;
            return t;
        }
        
        bool has_conflict(int qubit, const std::vector<SVGate>& layer)
        {
            for (const auto& gate_in_layer : layer) {
                if (gate_in_layer.op_name == OP::C1) {
                    if (gate_in_layer.qubit == qubit) return true;
                } else if (gate_in_layer.op_name == OP::C2) {
                    if (gate_in_layer.qubit == qubit || gate_in_layer.ctrl == qubit) return true;
                }
            }
            return false;
        }
        
        bool has_conflict(int qubit1, int qubit2, const std::vector<SVGate>& layer)
        {
            return has_conflict(qubit1, layer) || has_conflict(qubit2, layer);
        }
        
        void place_c1(
            const SVGate& s,
            std::vector<std::vector<SVGate>>& layers,
            std::map<int,int>& last_layer_map)
        {
            // determine earliest layer
            int q = s.qubit;
            int L = last_layer_map[q] + 1;
        
            // find first conflict-free layer
            while (true) {
                if (L > layers.size()) layers.resize(L);
                if (!has_conflict(q, layers[L - 1])) {
                    layers[L - 1].push_back(s);
                    last_layer_map[q] = L;
                    return;
                }
                L++;
            }
        }
        
        void place_c2(
            const SVGate& t, int a, int b,
            std::vector<std::vector<SVGate>>& layers,
            std::map<int,int>& last_layer_map)
        {
            // determine earliest layer
            int L = 1 + std::max(last_layer_map[a], last_layer_map[b]);
        
            SVGate x = t;
            x.ctrl = a;
            x.qubit = b;
        
            // find first conflict-free layer
            while (true) {
                if (L > layers.size()) layers.resize(L);
                if (!has_conflict(a, b, layers[L - 1])) {
                    layers[L - 1].push_back(x);
                    last_layer_map[a] = L;
                    last_layer_map[b] = L;
                    return;
                }
                L++;
            }
        }

    protected:
        IdxType n_qubits;
        IdxType* results = NULL;
        IdxType max_bond_dim;
        int block_size;
        double sv_cutoff;
        tamm::ExecutionHW exec_hw;

        tamm::ProcGroup pg_local_;
        tamm::ExecutionContext ec_local_;
        tamm::Scheduler sch_local_;

        tamm::ProcGroup pg;
        tamm::ExecutionContext ec;
        std::vector<IdxType> bond_dims;
        std::vector<IdxType> phys_dims;
        std::vector<tamm::TiledIndexSpace> bond_tis;
        std::vector<tamm::TiledIndexSpace> phys_tis;
        std::vector<tamm::Tensor<Cplx>> mps_tensors;
        IdxType* result = nullptr;
        CuCtx cu_ctx_;

        virtual void simulation_kernel(const std::vector<SVGate>& gates)
        {
            // separate parallel and sequential gates
            std::vector<SVGate> parallel_gates;
            std::vector<SVGate> sequential_gates;
            for (const auto& g : gates)
            {
                if (g.op_name == OP::C1 || g.op_name == OP::C2)
                {
                    parallel_gates.push_back(g);
                }
                else if (g.op_name == OP::M || g.op_name == OP::MA || g.op_name == OP::RESET)
                {
                    sequential_gates.push_back(g);
                }
            }
        
            pg.barrier();
        
            // process parallel gates
            if (!parallel_gates.empty())
            {
                std::vector<SVGate> flat_gates;
                flat_gates.reserve(parallel_gates.size() * 2);
        
                for (const auto& g : parallel_gates)
                {
                    if (g.op_name == OP::C1)
                    {
                        flat_gates.push_back(g);
                    }
                    else
                    {
                        int a = g.ctrl;
                        int b = g.qubit;
        
                        if (std::abs(b - a) > 1)
                        {
                            int start = std::min(a, b);
                            int end = std::max(a, b);
                            for (int k = start; k < end - 1; ++k)
                            {
                                flat_gates.push_back(make_swap_sv(k, k + 1));
                            }
                            flat_gates.push_back(make_local_c2_sv(g, end - 1, end));
                            for (int k = end - 2; k >= start; --k)
                            {
                                flat_gates.push_back(make_swap_sv(k, k + 1));
                            }
                        }
                        else
                        {
                            flat_gates.push_back(g);
                        }
                    }
                }
        
                std::vector<std::vector<SVGate>> layers;
                layers.reserve(flat_gates.size());
                std::map<int, int> last_layer_map;
        
                for (const auto& g : flat_gates)
                {
                    if (g.op_name == OP::C1)
                    {
                        place_c1(g, layers, last_layer_map);
                    }
                    else
                    {
                        place_c2(g, g.ctrl, g.qubit, layers, last_layer_map);
                    }
                }
        
                pg.barrier();
        
                // execute layers
                for (int layer_idx = 0; layer_idx < layers.size(); ++layer_idx)
                {
                    const auto& layer = layers[layer_idx];
                    if (layer.empty()) continue;
        
                    pg.barrier();
                    auto local_update_results = run_gates_parallel(layer);
        
                    pg.barrier();
                    apply_collective_updates(local_update_results);
                    pg.barrier();
                }
            }
        
            pg.barrier();
        
            // process sequential gates
            for (const auto& g : sequential_gates)
            {
                if (g.op_name == OP::RESET)
                {
                    RESET_GATE(g.qubit);
                }
                else if (g.op_name == OP::M)
                {
                    M_GATE(g.qubit);
                }
                else if (g.op_name == OP::MA)
                {
                    MA_GATE(g.qubit);
                }
            }
        
            pg.barrier();
        }

        std::vector<LocalGateResult> run_gates_parallel(const std::vector<SVGate>& batch)
        {
            // initialize atomic counter
            tamm::AtomicCounterGA gate_counter(pg, 1);
            gate_counter.allocate(0);
            pg.barrier();
        
            std::vector<LocalGateResult> local_results;
        
            // process gates in parallel
            while (true)
            {
                long long gate_idx = gate_counter.fetch_add(0, 1);
                if (gate_idx >= static_cast<long long>(batch.size())) break;
        
                const SVGate& g = batch[gate_idx];
        
                if (g.op_name == OP::C1)
                {
                    LocalGateResult result = C1_GATE_COMPUTE(g);
                    if (result.is_valid)
                    {
                        local_results.push_back(std::move(result));
                    }
                }
                else if (g.op_name == OP::C2)
                {
                    std::array<Cplx, 16> U4;
                    for (int i = 0; i < 16; ++i) U4[i] = Cplx(g.gm_real[i], g.gm_imag[i]);
        
                    LocalGateResult result = C2_GATE_COMPUTE(U4, g.ctrl, g.qubit);
                    if (result.is_valid)
                    {
                        local_results.push_back(std::move(result));
                    }
                }
            }
        
            // finalize
            pg.barrier();
            gate_counter.deallocate();
        
            return local_results;
        }

       std::vector<GateUpdateMetadata> allgather_metadata(const std::vector<LocalGateResult>& local_results)
        {
            // build local metadata
            std::vector<GateUpdateMetadata> local_metadata;
            local_metadata.reserve(local_results.size());
            for (const auto& res : local_results)
            {
                if (res.is_valid)
                {
                    local_metadata.push_back({
                        true,
                        res.q0,
                        res.q1,
                        res.new_bond_dim,
                        res.original_rank
                    });
                }
            }
        
            // gather sizes
            int local_size_bytes = local_metadata.size() * sizeof(GateUpdateMetadata);
            std::vector<int> all_sizes_bytes(pg.size().value());
            pg.allgather(&local_size_bytes, 1, all_sizes_bytes.data(), 1);
        
            // compute displacements
            std::vector<int> displacements_bytes(pg.size().value(), 0);
            int total_size_bytes = 0;
            for (size_t i = 0; i < all_sizes_bytes.size(); ++i)
            {
                if (i > 0)
                {
                    displacements_bytes[i] = displacements_bytes[i - 1] + all_sizes_bytes[i - 1];
                }
                total_size_bytes += all_sizes_bytes[i];
            }
        
            // gather metadata
            std::vector<GateUpdateMetadata> all_metadata;
            if (total_size_bytes > 0)
            {
                all_metadata.resize(total_size_bytes / sizeof(GateUpdateMetadata));
                MPI_Allgatherv(local_metadata.data(),
                               local_size_bytes,
                               MPI_BYTE,
                               all_metadata.data(),
                               all_sizes_bytes.data(),
                               displacements_bytes.data(),
                               MPI_BYTE,
                               pg.comm());
            }
        
            return all_metadata;
        }
 
        void populate_tensor_from_local_data(
            tamm::Tensor<Cplx>& global_tensor,
            const std::vector<Cplx>& local_data,
            const std::vector<size_t>& full_dims)
        {
            // iterate over tensor blocks
            for (const auto& blockid : global_tensor.loop_nest())
            {
                auto block_dims = global_tensor.block_dims(blockid);
                auto block_offsets = global_tensor.block_offsets(blockid);
                size_t block_size = global_tensor.block_size(blockid);
        
                std::vector<Cplx> block_buf(block_size);
        
                size_t c = 0;
                for (size_t i = block_offsets[0]; i < block_offsets[0] + block_dims[0]; ++i)
                {
                    for (size_t j = block_offsets[1]; j < block_offsets[1] + block_dims[1]; ++j)
                    {
                        for (size_t k = block_offsets[2]; k < block_offsets[2] + block_dims[2]; ++k, ++c)
                        {
                            size_t source_idx = i * (full_dims[1] * full_dims[2]) + j * full_dims[2] + k;
                            if (source_idx < local_data.size())
                            {
                                block_buf[c] = local_data[source_idx];
                            }
                        }
                    }
                }
        
                global_tensor.put(blockid, block_buf);
            }
        }

        void apply_collective_updates(std::vector<LocalGateResult>& local_results)
        {
            int rank = pg.rank().value();
        
            // gather metadata
            auto all_metadata = allgather_metadata(local_results);
            tamm::Scheduler sch_global{ec};
        
            // deallocate old tensors and update bond dimensions
            std::set<IdxType> deallocated_sites;
            for (const auto& meta : all_metadata)
            {
                if (!meta.is_valid) continue;
        
                if (deallocated_sites.find(meta.q0) == deallocated_sites.end())
                {
                    sch_global.deallocate(mps_tensors[meta.q0]);
                    deallocated_sites.insert(meta.q0);
                }
        
                if (meta.q1 != -1 && deallocated_sites.find(meta.q1) == deallocated_sites.end())
                {
                    sch_global.deallocate(mps_tensors[meta.q1]);
                    deallocated_sites.insert(meta.q1);
                }
        
                if (meta.q1 != -1)
                {
                    bond_dims[meta.q0 + 1] = meta.new_bond_dim;
                    tamm::IndexSpace is_new_bond{tamm::range(meta.new_bond_dim)};
                    bond_tis[meta.q0 + 1] = tamm::TiledIndexSpace(is_new_bond, block_size);
                }
            }
        
            // allocate new tensors
            std::map<IdxType, tamm::Tensor<Cplx>> site_to_new_tensor;
            for (const auto& site : deallocated_sites)
            {
                site_to_new_tensor.emplace(
                    site,
                    tamm::Tensor<Cplx>{bond_tis[site], phys_tis[site], bond_tis[site + 1]}
                );
                site_to_new_tensor.at(site).set_dense();
                sch_global.allocate(site_to_new_tensor.at(site));
            }
        
            sch_global.execute(exec_hw);
        
            // populate tensors with results
            int local_result_idx = 0;
            for (const auto& meta : all_metadata)
            {
                if (!meta.is_valid) continue;
        
                if (rank == meta.original_rank)
                {
                    auto& result_data = local_results[local_result_idx++];
                    auto& new_T0_ref = site_to_new_tensor.at(meta.q0);
        
                    std::vector<size_t> t0_full_dims = {
                        (size_t)bond_dims[meta.q0],
                        (size_t)phys_dims[meta.q0],
                        (size_t)bond_dims[meta.q0 + 1]
                    };
                    populate_tensor_from_local_data(new_T0_ref, result_data.new_T0_data, t0_full_dims);
        
                    if (meta.q1 != -1)
                    {
                        auto& new_T1_ref = site_to_new_tensor.at(meta.q1);
                        std::vector<size_t> t1_full_dims = {
                            (size_t)bond_dims[meta.q1],
                            (size_t)phys_dims[meta.q1],
                            (size_t)bond_dims[meta.q1 + 1]
                        };
                        populate_tensor_from_local_data(new_T1_ref, result_data.new_T1_data, t1_full_dims);
                    }
                }
            }
        
            // synchronize ranks
            pg.barrier();
        
            // update MPS state vector
            for (auto const& [site, new_tensor] : site_to_new_tensor)
            {
                mps_tensors[site] = new_tensor;
            }
        }


        LocalGateResult C2_GATE_COMPUTE(const std::array<Cplx, 16>& U4, IdxType q0, IdxType q1)
        {
            int rank = pg.rank().value();
        
            // create local tensors
            tamm::Tensor<Cplx> T0_local({bond_tis[q0], phys_tis[q0], bond_tis[q0 + 1]});
            tamm::Tensor<Cplx> T1_local({bond_tis[q1], phys_tis[q1], bond_tis[q1 + 1]});
            tamm::Tensor<Cplx> M_local({bond_tis[q0], phys_tis[q0], phys_tis[q1], bond_tis[q1 + 1]});
            tamm::Tensor<Cplx> G4_local({phys_tis[q0], phys_tis[q1], phys_tis[q0], phys_tis[q1]});
            tamm::Tensor<Cplx> M2_local({bond_tis[q0], phys_tis[q0], phys_tis[q1], bond_tis[q1 + 1]});
        
            T0_local.set_dense();
            T1_local.set_dense();
            M_local.set_dense();
            G4_local.set_dense();
            M2_local.set_dense();
        
            sch_local_.allocate(T0_local, T1_local, M_local, G4_local, M2_local).execute(exec_hw);
        
            // copy global tensors into local
            sch_local_
                (T0_local() = mps_tensors[q0]())
                (T1_local() = mps_tensors[q1]())
                .execute(exec_hw);
        
            // set gate matrix
            Cplx* g4_buf = G4_local.access_local_buf();
            for (size_t i = 0; i < 16; ++i)
            {
                g4_buf[i] = U4[i];
            }
        
            // perform contractions
            sch_local_
                (M_local("l","p0","p1","r") = T0_local("l","p0","b") * T1_local("b","p1","r"))
                .execute(exec_hw);
        
            sch_local_
                (M2_local("l","p0p","p1p","r") = G4_local("p0p","p1p","p0","p1") * M_local("l","p0","p1","r"))
                .execute(exec_hw);
        
            // perform SVD and reconstruct
            std::vector<Cplx> Ti_new_data, Tj_new_data;
            IdxType new_bond_dim = local_svd_and_reconstruct_data(M2_local, Ti_new_data, Tj_new_data, q0, q1);
        
            sch_local_.deallocate(T0_local, T1_local, M_local, G4_local, M2_local).execute(exec_hw);
        
            // package result
            LocalGateResult result;
            result.is_valid = true;
            result.q0 = q0;
            result.q1 = q1;
            result.new_bond_dim = new_bond_dim;
            result.new_T0_data = std::move(Ti_new_data);
            result.new_T1_data = std::move(Tj_new_data);
            result.original_rank = rank;
        
            return result;
        }

        LocalGateResult C1_GATE_COMPUTE(const SVGate& g)
        {
            int rank = pg.rank().value();
            IdxType q_idx = g.qubit;
        
            // extract target tensor spaces
            auto& target_tensor_global = mps_tensors[q_idx];
            auto tis_l = target_tensor_global.tiled_index_spaces()[0];
            auto tis_p = target_tensor_global.tiled_index_spaces()[1];
            auto tis_r = target_tensor_global.tiled_index_spaces()[2];
        
            // define local tensors
            tamm::Tensor<Cplx> T_in_local({tis_l, tis_p, tis_r});
            tamm::Tensor<Cplx> G_local({tis_p, tis_p});
            tamm::Tensor<Cplx> T_new_local({tis_l, tis_p, tis_r});
        
            T_in_local.set_dense();
            G_local.set_dense();
            T_new_local.set_dense();
        
            sch_local_.allocate(T_in_local, G_local, T_new_local).execute(exec_hw);
        
            // copy global tensor into local
            sch_local_(T_in_local() = target_tensor_global()).execute(exec_hw);
        
            // set gate matrix
            std::array<Cplx, 4> U;
            for (int i = 0; i < 4; ++i) U[i] = Cplx(g.gm_real[i], g.gm_imag[i]);
            Cplx* g_local_buf = G_local.access_local_buf();
            std::memcpy(g_local_buf, U.data(), 4 * sizeof(Cplx));
        
            // perform contraction
            sch_local_(T_new_local("l","p_prime","r") = G_local("p_prime","p") * T_in_local("l","p","r")).execute(exec_hw);
        
            // extract result tensor data
            std::vector<Cplx> t_out_data(T_new_local.size());
            Cplx* result_buffer_ptr = T_new_local.access_local_buf();
            std::memcpy(t_out_data.data(), result_buffer_ptr, t_out_data.size() * sizeof(Cplx));
        
            sch_local_.deallocate(T_in_local, T_new_local, G_local).execute(exec_hw);
        
            // package result
            LocalGateResult result;
            result.is_valid = true;
            result.q0 = q_idx;
            result.q1 = -1;
            result.new_bond_dim = bond_dims[q_idx + 1];
            result.new_T0_data = std::move(t_out_data);
            result.original_rank = rank;
        
            return result;
        }

        void gpu_svd_jacobi(
            const Cplx* A_h, int m, int n,
            std::vector<double>& S,
            std::vector<Cplx>& U_row,
            std::vector<Cplx>& VT_row)
        {
            int rank = pg.rank().value();
        
            // configure solver parameters
            cusolverDnXgesvdjSetTolerance(cu_ctx_.jp, 1e-3);
            cusolverDnXgesvdjSetMaxSweeps(cu_ctx_.jp, 5);
        
            int lda = m;
            int ldu = m;
            int ldv = n;
            int econ = 1;
            int k = std::min(m, n);
        
            // allocate device memory
            cuDoubleComplex* d_A = nullptr;
            double* d_S = nullptr;
            cuDoubleComplex* d_U = nullptr;
            cuDoubleComplex* d_V = nullptr;
            int* d_info = nullptr;
        
            cudaMalloc((void**)&d_A, sizeof(cuDoubleComplex) * (size_t)lda * (size_t)n);
            cudaMalloc((void**)&d_S, sizeof(double) * (size_t)k);
            cudaMalloc((void**)&d_U, sizeof(cuDoubleComplex) * (size_t)ldu * (size_t)k);
            cudaMalloc((void**)&d_V, sizeof(cuDoubleComplex) * (size_t)ldv * (size_t)k);
            cudaMalloc((void**)&d_info, sizeof(int));
        
            // copy input matrix to device
            cudaMemcpyAsync(d_A, reinterpret_cast<const cuDoubleComplex*>(A_h),
                            sizeof(cuDoubleComplex) * (size_t)lda * (size_t)n,
                            cudaMemcpyHostToDevice, cu_ctx_.stream);
        
            // query workspace size
            int lwork_req = 0;
            cusolverDnZgesvdj_bufferSize(cu_ctx_.solver, CUSOLVER_EIG_MODE_VECTOR, econ,
                                         m, n, d_A, lda, d_S, d_U, ldu, d_V, ldv,
                                         &lwork_req, cu_ctx_.jp);
        
            if (lwork_req > cu_ctx_.lwork_jac)
            {
                if (cu_ctx_.d_work_jac) cudaFree(cu_ctx_.d_work_jac);
                cu_ctx_.lwork_jac = lwork_req;
                cudaMalloc((void**)&cu_ctx_.d_work_jac, sizeof(cuDoubleComplex) * (size_t)cu_ctx_.lwork_jac);
            }
        
            // perform SVD
            cusolverDnZgesvdj(cu_ctx_.solver, CUSOLVER_EIG_MODE_VECTOR, econ,
                              m, n, d_A, lda, d_S, d_U, ldu, d_V, ldv,
                              reinterpret_cast<cuDoubleComplex*>(cu_ctx_.d_work_jac),
                              cu_ctx_.lwork_jac, d_info, cu_ctx_.jp);
        
            cudaStreamSynchronize(cu_ctx_.stream);
        
            // copy results to host
            S.resize((size_t)k);
            std::vector<Cplx> U_col((size_t)ldu * (size_t)k);
            std::vector<Cplx> V_col((size_t)ldv * (size_t)k);
        
            cudaMemcpy(S.data(), d_S, sizeof(double) * (size_t)k, cudaMemcpyDeviceToHost);
            cudaMemcpy(U_col.data(), d_U, sizeof(Cplx) * (size_t)ldu * (size_t)k, cudaMemcpyDeviceToHost);
            cudaMemcpy(V_col.data(), d_V, sizeof(Cplx) * (size_t)ldv * (size_t)k, cudaMemcpyDeviceToHost);
        
            // convert U to row-major
            U_row.resize((size_t)m * (size_t)k);
            for (int i = 0; i < m; ++i)
                for (int j = 0; j < k; ++j)
                    U_row[(size_t)i * (size_t)k + (size_t)j] = U_col[(size_t)i + (size_t)j * (size_t)ldu];
        
            // convert V to row-major with conjugation
            VT_row.resize((size_t)k * (size_t)n);
            for (int i = 0; i < k; ++i)
                for (int j = 0; j < n; ++j)
                    VT_row[(size_t)i * (size_t)n + (size_t)j] = std::conj(V_col[(size_t)j + (size_t)i * (size_t)ldv]);
        
            // free device memory
            cudaFree(d_info);
            cudaFree(d_V);
            cudaFree(d_U);
            cudaFree(d_S);
            cudaFree(d_A);
        }

        IdxType local_svd_and_reconstruct_data(
            tamm::Tensor<Cplx>& M2_local,
            std::vector<Cplx>& Ti_new_data,
            std::vector<Cplx>& Tj_new_data,
            IdxType q0, IdxType q1)
        {
            int rank = pg.rank().value();
        
            // extract dimensions
            const IdxType phys_dim = 2;
            IdxType Dl = M2_local.tiled_index_spaces()[0].index_space().num_indices();
            IdxType Dr = M2_local.tiled_index_spaces()[3].index_space().num_indices();
        
            int m = Dl * phys_dim;
            int n = phys_dim * Dr;
        
            // access local TAMM buffer
            Cplx* M2_hostbuf_rowmajor = M2_local.access_local_buf();
        
            // convert to column-major layout
            std::vector<Cplx> M2_col_major(m * n);
            size_t c = 0;
            for (size_t l = 0; l < Dl; ++l)
            {
                for (size_t p0 = 0; p0 < phys_dim; ++p0)
                {
                    for (size_t p1 = 0; p1 < phys_dim; ++p1)
                    {
                        for (size_t r = 0; r < Dr; ++r, ++c)
                        {
                            size_t row = l * phys_dim + p0;
                            size_t col = p1 * Dr + r;
                            M2_col_major[col * m + row] = M2_hostbuf_rowmajor[c];
                        }
                    }
                }
            }
        
            // perform GPU SVD
            std::vector<double> S_vals;
            std::vector<Cplx> U_mat_rowmajor, VT_mat_rowmajor;
            gpu_svd_jacobi(M2_col_major.data(), m, n, S_vals, U_mat_rowmajor, VT_mat_rowmajor);
        
            // determine truncation
            std::vector<IdxType> keep_indices;
            keep_indices.reserve(S_vals.size());
            for (size_t i = 0; i < S_vals.size(); ++i)
            {
                if (S_vals[i] >= sv_cutoff)
                {
                    keep_indices.push_back(i);
                }
            }
            IdxType chi = std::min<IdxType>(max_bond_dim, IdxType(keep_indices.size()));
            if (chi == 0 && !S_vals.empty())
            {
                chi = 1;
            }
        
            // build left tensor data
            Ti_new_data.resize(Dl * phys_dim * chi);
            c = 0;
            for (size_t l = 0; l < Dl; ++l)
            {
                for (size_t p0 = 0; p0 < phys_dim; ++p0)
                {
                    for (size_t b = 0; b < chi; ++b, ++c)
                    {
                        size_t row = l * phys_dim + p0;
                        size_t col = keep_indices[b];
                        Ti_new_data[c] = U_mat_rowmajor[row * S_vals.size() + col];
                    }
                }
            }
        
            // build right tensor data
            Tj_new_data.resize(chi * phys_dim * Dr);
            c = 0;
            for (size_t b = 0; b < chi; ++b)
            {
                for (size_t p1 = 0; p1 < phys_dim; ++p1)
                {
                    for (size_t r = 0; r < Dr; ++r, ++c)
                    {
                        size_t row = keep_indices[b];
                        size_t col = p1 * Dr + r;
                        Tj_new_data[c] = Cplx(S_vals[row], 0.0) * VT_mat_rowmajor[row * n + col];
                    }
                }
            }
        
            return chi;
        }

        void right_canonicalize(std::vector<tamm::Tensor<Cplx>> &MPS)
        {
            // canonicalize MPS from right end toward left
            for (IdxType i = n_qubits - 1; i > 0; --i)
            {
                // extract dimensions and assemble Mmat for SVD
                IdxType Dl_old = bond_dims[i];
                IdxType Dr     = bond_dims[i + 1];
                IdxType d      = phys_dims[i];
                Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Mmat(Dl_old, d * Dr);
                {
                    auto &T = MPS[i];
                    for (const auto &blockid : T.loop_nest())
                    {
                        const size_t bs = T.block_size(blockid);
                        std::vector<Cplx> hostbuf(bs);
                        T.get(blockid, hostbuf);
                        auto dims = T.block_dims(blockid);
                        auto offs = T.block_offsets(blockid);
                        size_t idx = 0;
                        for (size_t ll = offs[0]; ll < offs[0] + dims[0]; ++ll)
                        {
                            for (size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
                            {
                                for (size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx)
                                {
                                    Mmat(ll, pp * Dr + rr) = hostbuf[idx];
                                }
                            }
                        }
                    }
                }
        
                // perform SVD and truncate to bond dimension
                Eigen::BDCSVD<decltype(Mmat)> svd(Mmat, Eigen::ComputeThinU | Eigen::ComputeThinV);
                auto svals = svd.singularValues();
                IdxType chi = std::min<IdxType>(IdxType(svals.size()), max_bond_dim);
                auto Umat  = svd.matrixU().leftCols(chi);
                auto Vh    = svd.matrixV().leftCols(chi).adjoint();
                bond_dims[i] = chi;
                {
                    tamm::IndexSpace is_new{ tamm::range(chi) };
                    bond_tis[i] = tamm::TiledIndexSpace(is_new, block_size);
                }
        
                // build updated right tensor via Vh
                tamm::Tensor<Cplx> Tnew({ bond_tis[i], phys_tis[i], bond_tis[i + 1] });
                Tnew.set_dense();
                Tnew.allocate(&ec);
                {
                    auto &T = Tnew;
                    for (const auto &blockid : T.loop_nest())
                    {
                        const size_t bs = T.block_size(blockid);
                        std::vector<Cplx> hostbuf(bs);
                        auto dims = T.block_dims(blockid);
                        auto offs = T.block_offsets(blockid);
                        size_t idx = 0;
                        for (size_t ll = offs[0]; ll < offs[0] + dims[0]; ++ll)
                        {
                            for (size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
                            {
                                for (size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx)
                                {
                                    hostbuf[idx] = Vh(ll, pp * Dr + rr);
                                }
                            }
                        }
                        T.put(blockid, hostbuf);
                    }
                }
        
                // update left neighbor via Umat * S
                IdxType Dl_prev = bond_dims[i - 1];
                IdxType d_prev  = phys_dims[i - 1];
                Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Mprev(Dl_prev * d_prev, Dl_old);
                {
                    auto &Told = MPS[i - 1];
                    for (const auto &blockid : Told.loop_nest())
                    {
                        const size_t bs = Told.block_size(blockid);
                        std::vector<Cplx> hostbuf(bs);
                        Told.get(blockid, hostbuf);
                        auto dims = Told.block_dims(blockid);
                        auto offs = Told.block_offsets(blockid);
                        size_t idx = 0;
                        for (size_t ll = offs[0]; ll < offs[0] + dims[0]; ++ll)
                        {
                            for (size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
                            {
                                for (size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx)
                                {
                                    Mprev(ll * d_prev + pp, rr) = hostbuf[idx];
                                }
                            }
                        }
                    }
                }
                auto US     = Umat * svals.head(chi).asDiagonal();
                auto Mprev2 = Mprev * US;
        
                // build updated left tensor via Mprev2
                tamm::Tensor<Cplx> Tprev({ bond_tis[i - 1], phys_tis[i - 1], bond_tis[i] });
                Tprev.set_dense();
                Tprev.allocate(&ec);
                {
                    auto &T = Tprev;
                    for (const auto &blockid : T.loop_nest())
                    {
                        const size_t bs = T.block_size(blockid);
                        std::vector<Cplx> hostbuf(bs);
                        auto dims = T.block_dims(blockid);
                        auto offs = T.block_offsets(blockid);
                        size_t idx = 0;
                        for (size_t ll = offs[0]; ll < offs[0] + dims[0]; ++ll)
                        {
                            for (size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
                            {
                                for (size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx)
                                {
                                    hostbuf[idx] = Mprev2(ll * d_prev + pp, rr);
                                }
                            }
                        }
                        T.put(blockid, hostbuf);
                    }
                }
        
                // replace tensors in MPS
                MPS[i].deallocate();
                MPS[i - 1].deallocate();
                MPS[i]     = std::move(Tnew);
                MPS[i - 1] = std::move(Tprev);
            }
        }

        void left_canonicalize(std::vector<tamm::Tensor<Cplx>>& MPS)
        {
            for (IdxType i = 0; i < n_qubits - 1; ++i)
            {
                // assemble matrix from MPS[i]
                IdxType Dl     = bond_dims[i];
                IdxType Dr_old = bond_dims[i + 1];
                IdxType d      = phys_dims[i];
                Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Mmat(Dl * d, Dr_old);
                {
                    auto& T = MPS[i];
                    for (const auto& blockid : T.loop_nest())
                    {
                        size_t bs = T.block_size(blockid);
                        std::vector<Cplx> hostbuf(bs);
                        T.get(blockid, hostbuf);
        
                        auto dims = T.block_dims(blockid);
                        auto offs = T.block_offsets(blockid);
                        size_t idx = 0;
                        for (size_t ll = offs[0]; ll < offs[0] + dims[0]; ++ll)
                        {
                            for (size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
                            {
                                for (size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx)
                                {
                                    Mmat(ll * d + pp, rr) = hostbuf[idx];
                                }
                            }
                        }
                    }
                }
        
                // compute truncated SVD of Mmat
                Eigen::BDCSVD<decltype(Mmat)> svd(Mmat, Eigen::ComputeThinU | Eigen::ComputeThinV);
                auto svals = svd.singularValues();
                IdxType chi = std::min<IdxType>(IdxType(svals.size()), max_bond_dim);
                auto Umat  = svd.matrixU().leftCols(chi);
                auto Sdiag  = svals.head(chi).asDiagonal();
                auto Vh     = svd.matrixV().leftCols(chi).adjoint();
        
                // update bond dimension and index space
                bond_dims[i + 1] = chi;
                {
                    tamm::IndexSpace is_new{ tamm::range(chi) };
                    bond_tis[i + 1] = tamm::TiledIndexSpace(is_new, block_size);
                }
        
                // build new left tensor from Umat
                tamm::Tensor<Cplx> Tleft({ bond_tis[i], phys_tis[i], bond_tis[i + 1] });
                Tleft.set_dense();
                Tleft.allocate(&ec);
                for (const auto& blockid : Tleft.loop_nest())
                {
                    size_t bs = Tleft.block_size(blockid);
                    std::vector<Cplx> hostbuf(bs);
        
                    auto dims = Tleft.block_dims(blockid);
                    auto offs = Tleft.block_offsets(blockid);
                    size_t idx = 0;
                    for (size_t ll = offs[0]; ll < offs[0] + dims[0]; ++ll)
                    {
                        for (size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
                        {
                            for (size_t bb = offs[2]; bb < offs[2] + dims[2]; ++bb, ++idx)
                            {
                                hostbuf[idx] = Umat(ll * d + pp, bb);
                            }
                        }
                    }
        
                    Tleft.put(blockid, hostbuf);
                }
        
                // assemble matrix from MPS[i+1]
                IdxType d_next   = phys_dims[i + 1];
                IdxType Dr_right = bond_dims[i + 2];
                Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Tnext_mat(chi, d_next * Dr_right);
                {
                    auto& Tnext_old = MPS[i + 1];
                    for (const auto& blockid : Tnext_old.loop_nest())
                    {
                        size_t bs = Tnext_old.block_size(blockid);
                        std::vector<Cplx> hostbuf(bs);
                        Tnext_old.get(blockid, hostbuf);
        
                        auto dims = Tnext_old.block_dims(blockid);
                        auto offs = Tnext_old.block_offsets(blockid);
                        size_t idx = 0;
                        for (size_t bb = offs[0]; bb < offs[0] + dims[0]; ++bb)
                        {
                            for (size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
                            {
                                for (size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx)
                                {
                                    Tnext_mat(bb, pp * Dr_right + rr) = hostbuf[idx];
                                }
                            }
                        }
                    }
                }
        
                // apply spectrum to right neighbor tensor
                auto M2_eig    = Sdiag * Vh;
                auto Tnext_new = M2_eig * Tnext_mat;
                tamm::Tensor<Cplx> Tnext({ bond_tis[i + 1], phys_tis[i + 1], bond_tis[i + 2] });
                Tnext.set_dense();
                Tnext.allocate(&ec);
                for (const auto& blockid : Tnext.loop_nest())
                {
                    size_t bs = Tnext.block_size(blockid);
                    std::vector<Cplx> hostbuf(bs);
        
                    auto dims = Tnext.block_dims(blockid);
                    auto offs = Tnext.block_offsets(blockid);
                    size_t idx = 0;
                    for (size_t bb = offs[0]; bb < offs[0] + dims[0]; ++bb)
                    {
                        for (size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
                        {
                            for (size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx)
                            {
                                hostbuf[idx] = Tnext_new(bb, pp * Dr_right + rr);
                            }
                        }
                    }
        
                    Tnext.put(blockid, hostbuf);
                }
        
                // replace tensors in MPS
                MPS[i].deallocate();
                MPS[i + 1].deallocate();
                MPS[i]     = std::move(Tleft);
                MPS[i + 1] = std::move(Tnext);
            }
        }

        /*MA Gate
        *1) Free and reallocate a buffer to store repetition number of measurement results.

        *2) Canonicalize the MPS from both left and right to stabilize subsequent contractions.
        
        *3) Initialize a random number generator for sampling measurement outcomes.
        
        *4) For each repetition:
        *a. Initialize the environment vector with amplitude 1.
        *b. Loop over all qubit sites from left to right:
        *i. Extract the MPS tensor at the current site.
        *ii. Contract the current environment with the tensor to produce two branch environments for qubit outcomes 0 and 1.
        *iii. Compute squared norms of both branches to determine measurement probabilities.
        *iv. Sample a measurement outcome based on these probabilities.
        *v. Record the sampled bit in the packed result.
        *vi. Renormalize the chosen branch and set it as the new environment.
        
        *5) Store the final bitstring for each repetition in the results buffer. */
        virtual void MA_GATE(const IdxType repetition)
        {
            // allocate result buffer
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
            std::memset(results, 0, sizeof(IdxType) * repetition);
        
            // bring MPS into canonical form
            left_canonicalize(mps_tensors);
            right_canonicalize(mps_tensors);
        
            // perform measurement assignment repetitions
            std::mt19937_64 rng{std::random_device{}()};
            std::uniform_real_distribution<double> dist(0.0, 1.0);
            for (IdxType rep = 0; rep < repetition; ++rep)
            {
                IdxType packed = 0;
                std::vector<Cplx> env(1, Cplx{1.0, 0.0});
        
                // sweep through all qubits
                for (IdxType site = 0; site < n_qubits; ++site)
                {
                    auto& T = mps_tensors[site];
                    IdxType Dl = bond_dims[site];
                    IdxType Dr = bond_dims[site + 1];
        
                    std::vector<Cplx> env0(Dr, Cplx{0.0, 0.0});
                    std::vector<Cplx> env1(Dr, Cplx{0.0, 0.0});
        
                    // unpack tensor into branch environments
                    for (const auto& blockid : T.loop_nest())
                    {
                        size_t bs = T.block_size(blockid);
                        std::vector<Cplx> hostbuf(bs);
                        T.get(blockid, hostbuf);
        
                        auto dims = T.block_dims(blockid);
                        auto offs = T.block_offsets(blockid);
                        size_t idx = 0;
                        for (size_t jj = 0; jj < dims[0]; ++jj)
                        {
                            for (size_t kk = 0; kk < dims[1]; ++kk)
                            {
                                for (size_t ll = 0; ll < dims[2]; ++ll, ++idx)
                                {
                                    IdxType l = offs[0] + jj;
                                    IdxType s = offs[1] + kk;
                                    IdxType r = offs[2] + ll;
                                    Cplx prod = env[l] * hostbuf[idx];
                                    if (s == 0)
                                    {
                                        env0[r] += prod;
                                    }
                                    else
                                    {
                                        env1[r] += prod;
                                    }
                                }
                            }
                        }
                    }
        
                    // compute probabilities and sample outcome
                    long double w0 = 0.0L;
                    long double w1 = 0.0L;
                    for (IdxType r = 0; r < Dr; ++r)
                    {
                        w0 += std::norm(env0[r]);
                        w1 += std::norm(env1[r]);
                    }
                    long double sumw = w0 + w1;
                    long double p0 = sumw > 0.0L ? (w0 / sumw) : 0.0L;
                    bool outcome1 = (dist(rng) >= static_cast<double>(p0));
                    if (outcome1)
                    {
                        packed |= (IdxType(1) << site);
                    }
        
                    auto& chosen = outcome1 ? env1 : env0;
                    long double norm_branch = outcome1 ? w1 : w0;
                    long double invnorm = norm_branch > 0.0L
                        ? (1.0L / std::sqrt(norm_branch))
                        : 0.0L;
                    env.assign(Dr, Cplx{0.0, 0.0});
                    for (IdxType r = 0; r < Dr; ++r)
                    {
                        env[r] = chosen[r] * static_cast<Cplx>(invnorm);
                    }
                }
        
                results[rep] = packed;
            }
        }

        void local_left_step(IdxType i)
        {
            // extract dimensions for sites i and i+1
            const IdxType Dl       = bond_dims[i];
            const IdxType Dr       = bond_dims[i+1];
            const IdxType d        = phys_dims[i];
            const IdxType d_next   = phys_dims[i+1];
            const IdxType Dr_next  = bond_dims[i+2];
        
            // reshape MPS[i] into matrix Mmat of size (Dl*d) × Dr
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Mmat(Dl * d, Dr);
            auto& Ti = mps_tensors[i];
            Ti.loop_nest().iterate([&](auto const& idxs)
            {
                IdxType l = idxs[0], p = idxs[1], r = idxs[2];
                Cplx v;
                Ti.get(idxs, gsl::span<Cplx>(&v,1));
                Mmat(l * d + p, r) = v;
            });
        
            // compute SVD of Mmat
            Eigen::BDCSVD<decltype(Mmat)> svd(
                Mmat, Eigen::ComputeThinU | Eigen::ComputeThinV);
            auto Umat  = svd.matrixU();
            auto Sdiag = svd.singularValues().asDiagonal();
            auto Vh    = svd.matrixV().adjoint();
        
            // update bond dimension at i+1
            const IdxType chi = Umat.cols();
            bond_dims[i+1] = chi;
            {
                tamm::IndexSpace is_new{ tamm::range(chi) };
                bond_tis[i+1] = tamm::TiledIndexSpace(is_new, 1);
            }
        
            // write back new left tensor Ti_new
            tamm::Tensor<Cplx> Ti_new({ bond_tis[i], phys_tis[i], bond_tis[i+1] });
            Ti_new.set_dense();
            Ti_new.allocate(&ec);
            Ti_new.loop_nest().iterate([&](auto const& idxs)
            {
                IdxType l = idxs[0], p = idxs[1], b = idxs[2];
                Cplx val = Umat(l * d + p, b);
                Ti_new.put(idxs, gsl::span<Cplx>(&val,1));
            });
        
            // absorb spectrum S·Vh into MPS[i+1]
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> A = Sdiag * Vh;
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Bmat(
                Dr, d_next * Dr_next);
            auto& Tj = mps_tensors[i+1];
            Tj.loop_nest().iterate([&](auto const& idxs)
            {
                IdxType b = idxs[0], p = idxs[1], r = idxs[2];
                Cplx v;
                Tj.get(idxs, gsl::span<Cplx>(&v,1));
                Bmat(b, p * Dr_next + r) = v;
            });
            auto Bnew = A * Bmat;
        
            // write back new right tensor Tj_new
            tamm::Tensor<Cplx> Tj_new({ bond_tis[i+1], phys_tis[i+1], bond_tis[i+2] });
            Tj_new.set_dense();
            Tj_new.allocate(&ec);
            Tj_new.loop_nest().iterate([&](auto const& idxs)
            {
                IdxType b = idxs[0], p = idxs[1], r = idxs[2];
                Cplx val = Bnew(b, p * Dr_next + r);
                Tj_new.put(idxs, gsl::span<Cplx>(&val,1));
            });
        
            // replace old tensors
            Ti.deallocate();
            Tj.deallocate();
            mps_tensors[i]   = std::move(Ti_new);
            mps_tensors[i+1] = std::move(Tj_new);
        }
        
        /* Sets the MPS in a mixed gauge centered around a specified position */
        void position(IdxType site)
        {
            assert(site < n_qubits);
        
            // bring orthogonality center to qubit 0
            right_canonicalize(mps_tensors);
        
            // sweep center from 0 to target site
            for (IdxType i = 0; i < site; ++i)
            {
                local_left_step(i);
            }
        }

        virtual void M_GATE(const IdxType qubit)
        {
            throw std::runtime_error("Not implemented");
        }

        virtual void EXPECT_GATE(ObservableList *o)
        {
            throw std::runtime_error("Not implemented");
        }

        virtual void RESET_GATE(const IdxType qubit)
        {
            throw std::runtime_error("Not implemented");
        }

        static tamm::ProcGroup init_pg()
        {
            int argc = 0; char** argv = nullptr;
            return tamm::ProcGroup::create_world_coll();
        }
    };
} // namespace NWQS
