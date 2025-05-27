#pragma once

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

#include <Eigen/Dense>
#include <gsl/span>
#include <iostream>

#include <complex>

namespace NWQSim
{
    using Cplx = std::complex<ValType>;
    using Eigen::Index;
    class TN_CUDA : public QuantumState
    {
    public:
        TN_CUDA(IdxType _n_qubits,
               IdxType _max_bond_dim = 100)
        : QuantumState(SimType::TN),
          n_qubits(_n_qubits),
          block_size(2048),
          max_bond_dim(_max_bond_dim),
          pg(init_pg()),
          ec(pg, tamm::DistributionKind::dense, tamm::MemoryManagerKind::ga)
        {
          // 1) build bond index spaces, each tiled by block_size
          bond_tis.resize(n_qubits + 1);
          bond_dims.resize(n_qubits + 1);
          for(IdxType i = 0; i <= n_qubits; ++i) {
            bond_dims[i] = 1;
            tamm::IndexSpace is{tamm::range(1)};
            bond_tis[i] = tamm::TiledIndexSpace(is, block_size);
          }
        
          // 2) build physical index spaces (always small, tile-size=1)
          phys_tis.resize(n_qubits);
          phys_dims.resize(n_qubits);
          for(IdxType i = 0; i < n_qubits; ++i) {
            phys_dims[i] = 2;
            tamm::IndexSpace is{tamm::range(2)};
            phys_tis[i] = tamm::TiledIndexSpace(is, /*phys tilesize=*/1);
          }
        
          // 3) allocate and initialize MPS tensors via bulk puts
          mps_tensors.reserve(n_qubits);
          for(IdxType i = 0; i < n_qubits; ++i) {
            tamm::Tensor<Cplx> T({ bond_tis[i],
                                   phys_tis[i],
                                   bond_tis[i+1] });
            T.set_dense();
            T.allocate(&ec);
        
            // initialize each tile in one call
            for(const auto& blockid : T.loop_nest()) {
              const size_t bs = T.block_size(blockid);
              std::vector<Cplx> hostbuf(bs);
        
              auto dims    = T.block_dims(blockid);    // e.g. [d_l, d_p, d_r]
              auto offsets = T.block_offsets(blockid); // e.g. [o_l, o_p, o_r]
        
              for(size_t idx = 0; idx < bs; ++idx) {
                size_t rem = idx;
        
                // local coords within this block
                size_t l_loc = rem % dims[0]; rem /= dims[0];
                size_t p_loc = rem % dims[1]; rem /= dims[1];
                size_t r_loc = rem % dims[2];
        
                // global indices
                IdxType l = offsets[0] + l_loc;
                IdxType p = offsets[1] + p_loc;
                IdxType r = offsets[2] + r_loc;
        
                // set the |0…0> amplitude to 1, rest 0
                hostbuf[idx] = (l == 0 && p == 0 && r == 0)
                               ? Cplx{1.0,0.0}
                               : Cplx{0.0,0.0};
              }
        
              T.put(blockid, hostbuf);
            }
        
            mps_tensors.push_back(std::move(T));
          }
        }

         //Virtual destructor inherits from QuantumState
        ~TN_CUDA() noexcept override 
        {
            SAFE_FREE_HOST(results);
            //GA_Terminate();
            //MPI_Finalize();
            //tamm::finalize();
        }

        void reset_state() override {
            //printf("Inside reset gate\n");
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
            throw std::runtime_error("TN_CUDA does not use RNG seed, not accessible form cutensornet API\n");
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
            //printf("Inside sim\n");
            IdxType origional_gates = circuit->num_gates();
            std::vector<SVGate> gates = fuse_circuit_sv(circuit);
            IdxType n_gates = gates.size();
            assert(circuit->num_qubits() == n_qubits);
            double sim_time;
            cpu_timer sim_timer;
            sim_timer.start_timer();
	    
            //printf("Going inside sim kernel\n");
            simulation_kernel(gates);


            sim_timer.stop_timer();
            sim_time = sim_timer.measure();

            if (Config::PRINT_SIM_TRACE)
            {
                printf("\n============== TN-Sim ===============\n");
                printf("n_qubits:%lld, n_gates:%lld, sim_gates:%lld, ncpus:%lld, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_cpu:%.3lf MB\n",
                       n_qubits, origional_gates, n_gates, sim_time, 0.,
                       sim_time);
                printf("=====================================\n");
            }
        }


        IdxType* get_results() override
        {
            return results;
        }

        IdxType measure(IdxType qubit) override
        {
            throw std::runtime_error("TN_CUDA::measure not implemented");
        }

        IdxType* measure_all(IdxType repetition) override 
        {
            MA_GATE(repetition);
	        std::cout<<"measure_all was called"<<std::endl;
            return results; 
        }

        ValType* get_real() const override
        {
            throw std::runtime_error("TN_CUDA::get_real not implemented");
        }

        ValType* get_imag() const override
        {
            throw std::runtime_error("TN_CUDA::get_imag not implemented");
        }

        ValType get_exp_z() override
        {
            throw std::runtime_error("TN_CUDA::get_exp_z() not implemented");
        }

        ValType get_exp_z(const std::vector<size_t>& in_bits) override
        {
            throw std::runtime_error("TN_CUDA::get_exp_z(bits) not implemented");
        }

        void print_res_state() override
        {
            throw std::runtime_error("TN_CUDA::print_res_state not implemented");
        }

    protected:
        IdxType n_qubits;
        IdxType* results = NULL;
        IdxType max_bond_dim;
        int block_size;

        tamm::ProcGroup pg;
        tamm::ExecutionContext ec;
        std::vector<IdxType> bond_dims;
        std::vector<IdxType> phys_dims;
        std::vector<tamm::TiledIndexSpace> bond_tis;
        std::vector<tamm::TiledIndexSpace> phys_tis;
        std::vector<tamm::Tensor<Cplx>> mps_tensors;
        IdxType* result = nullptr;


        virtual void simulation_kernel(const std::vector<SVGate> &gates)
        {
            using namespace std::chrono;
        
            duration<double> c1_total(0);
            duration<double> c2_local_total(0);
            duration<double> c2_nonlocal_total(0);
            duration<double> ma_total(0);
        
            int c1_count = 0;
            int c2_local_count = 0;
            int c2_nonlocal_count = 0;
            int ma_count = 0;
        
            int n_gates = gates.size();
            for (int i = 0; i < n_gates; i++)
            {
                auto g = gates[i];
                auto op_start = steady_clock::now();
        
                if (g.op_name == OP::C1)
                {
                    std::array<Cplx,4> U;
                    for (int idx = 0; idx < 4; ++idx)
                    {
                        U[idx] = Cplx(g.gm_real[idx], g.gm_imag[idx]);
                    }
                    C1_GATE(U, g.qubit);
        
                    auto op_end = steady_clock::now();
                    c1_total += duration_cast<duration<double>>(op_end - op_start);
                    ++c1_count;
                }
                else if (g.op_name == OP::C2)
                {
                    std::array<Cplx,16> U4;
                    for (int idx = 0; idx < 16; ++idx)
                    {
                        U4[idx] = Cplx(g.gm_real[idx], g.gm_imag[idx]);
                    }
                    C2_GATE(U4, g.ctrl, g.qubit);
        
                    auto op_end = steady_clock::now();
                    int delta = std::abs(g.ctrl - g.qubit);
                    if (delta == 1)
                    {
                        c2_local_total += duration_cast<duration<double>>(op_end - op_start);
                        ++c2_local_count;
                    }
                    else
                    {
                        c2_nonlocal_total += duration_cast<duration<double>>(op_end - op_start);
                        ++c2_nonlocal_count;
                    }
                }
                else if (g.op_name == OP::RESET)
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
        
                    auto op_end = steady_clock::now();
                    ma_total += duration_cast<duration<double>>(op_end - op_start);
                    ++ma_count;
                }
                else if (g.op_name == OP::EXPECT)
                {
                    ObservableList *o = static_cast<ObservableList *>(g.data);
                    EXPECT_GATE(o);
                }
                else
                {
                    std::cout << "Unrecognized gates" << std::endl
                              << OP_NAMES[g.op_name] << std::endl;
                    throw std::logic_error("Invalid gate type");
                }
            }
        
            std::cout << "Gate Timing Report (seconds):" << std::endl;
        
            std::cout << "C1 (1-qubit) total:   " << c1_total.count()
                      << "    avg: " << (c1_count ? c1_total.count() / c1_count : 0) << std::endl;
        
            std::cout << "C2 local total:        " << c2_local_total.count()
                      << "    avg: " << (c2_local_count ? c2_local_total.count() / c2_local_count : 0) << std::endl;
        
            std::cout << "C2 non-local total:    " << c2_nonlocal_total.count()
                      << "    avg: " << (c2_nonlocal_count ? c2_nonlocal_total.count() / c2_nonlocal_count : 0) << std::endl;
        
            std::cout << "MA total:              " << ma_total.count()
                      << "    avg: " << (ma_count ? ma_total.count() / ma_count : 0) << std::endl;
        }

       void C1_GATE(const std::array<Cplx,4>& U, IdxType site) {
          using Clock = std::chrono::high_resolution_clock;
          using MS    = std::chrono::milliseconds;
        
          // 1) Build the 2×2 gate tensor G by bulk‐putting each tile
          auto t0 = Clock::now();
          tamm::Tensor<Cplx> G({ phys_tis[site], phys_tis[site] });
          G.set_dense();
          G.allocate(&ec);
        
          for(const auto& blockid : G.loop_nest()) {
            size_t bs = G.block_size(blockid);
            std::vector<Cplx> hostbuf(bs);
        
            auto dims    = G.block_dims(blockid);
            auto offsets = G.block_offsets(blockid);
        
            // each element in the tile corresponds to some (p_out,p_in)
            for(size_t j = 0; j < bs; ++j) {
              size_t rem = j;
              size_t pout_loc = rem % dims[0]; rem /= dims[0];
              size_t pin_loc  = rem % dims[1];
        
              // global indices
              size_t pout = offsets[0] + pout_loc;
              size_t pin  = offsets[1] + pin_loc;
        
              hostbuf[j] = U[pout*2 + pin];
            }
        
            G.put(blockid, hostbuf);
          }
          auto t1 = Clock::now();
          std::cout << "[C1_GATE] Build G —      "
                    << std::chrono::duration_cast<MS>(t1 - t0).count()
                    << " ms\n";
        
          // 2) Apply it
          auto t2 = Clock::now();
          auto& T = mps_tensors[site];
          tamm::Tensor<Cplx> Tnew({ bond_tis[site], phys_tis[site], bond_tis[site+1] });
          Tnew.set_dense();
          Tnew.allocate(&ec);
        
          tamm::Scheduler sch{ec};
          sch(Tnew("l","p'","r") = G("p'","p") * T("l","p","r"),
              "apply_one_qubit", tamm::ExecutionHW::GPU);
          sch.execute(tamm::ExecutionHW::GPU);
          auto t3 = Clock::now();
          std::cout << "[C1_GATE] Apply —        "
                    << std::chrono::duration_cast<MS>(t3 - t2).count()
                    << " ms\n";
        
          // 3) Cleanup
          T.deallocate();
          mps_tensors[site] = std::move(Tnew);
          G.deallocate();
        
          auto t4 = Clock::now();
          std::cout << "[C1_GATE] Total —        "
                    << std::chrono::duration_cast<MS>(t4 - t0).count()
                    << " ms\n";
        }
 
        std::vector<Cplx> flatten_mps_state() {
            const size_t N = mps_tensors.size();
            // 1) Initialize with tensor 0: shape = [d0, χ1]
            const IdxType d0   = phys_dims[0];
            const IdxType chi1 = bond_dims[1];
            std::vector<Cplx> state(d0 * chi1, Cplx{0,0});
            mps_tensors[0].loop_nest().iterate([&](auto const& idxs){
                // idxs = [l=0, p0∈[0,d0), b1∈[0,χ1))
                Cplx v; 
                mps_tensors[0].get(idxs, gsl::span<Cplx>(&v,1));
                state[idxs[1] * chi1 + idxs[2]] = v;
            });
        
            // 2) For each subsequent tensor i=1..N-1, absorb into 'state'
            for(size_t i = 1; i < N; ++i) {
                const auto& T      = mps_tensors[i];
                const IdxType di   = phys_dims[i];
                const IdxType chi_i   = bond_dims[i];
                const IdxType chi_ip1 = bond_dims[i+1];
        
                const size_t rows = state.size() / chi_i;
                // new_rows = rows * di, new_cols = χ_{i+1}
                std::vector<Cplx> new_state(rows * di * chi_ip1, Cplx{0,0});
        
                // Pre-copy T into host‐linear Tdata[(bi*di + pi)*χ_{i+1} + bip1]
                std::vector<Cplx> Tdata(chi_i * di * chi_ip1);
                T.loop_nest().iterate([&](auto const& idxs){
                    Cplx v;
                    T.get(idxs, gsl::span<Cplx>(&v,1));
                    const auto bi   = idxs[0];
                    const auto pi   = idxs[1];
                    const auto bip1 = idxs[2];
                    Tdata[(bi * di + pi) * chi_ip1 + bip1] = v;
                });
        
                // Contract: new_state[(r*di + pi)*χ_{i+1} + bip1] +=
                //     state[r*χ_i + bi] * Tdata[(bi*di + pi)*χ_{i+1} + bip1]
                for(size_t r = 0; r < rows; ++r) {
                    for(IdxType bi = 0; bi < chi_i; ++bi) {
                        const Cplx α = state[r * chi_i + bi];
                        if(α == Cplx{0,0}) continue;
                        for(IdxType pi = 0; pi < di; ++pi) {
                            for(IdxType bip1 = 0; bip1 < chi_ip1; ++bip1) {
                                new_state[(pi * rows + r) * chi_ip1 + bip1]
                                  += α * Tdata[(bi*di + pi) * chi_ip1 + bip1];
                            }
                        }
                    }
                }
        
                state.swap(new_state);
            }
        
            // At the end χ_N == 1, so state.size() == ∏ₖ dₖ
            return state;
        }

        void dump_state(const std::string& tag) {
            auto psi = flatten_mps_state();
            const size_t N          = phys_dims.size();      // number of qubits
            const size_t dim        = psi.size();            // should be 2^N
            printf("[DUMP %s] dim = %zu\n", tag.c_str(), dim);
        
            // For each basis index i = 0..2^N-1:
            for(size_t i = 0; i < dim; ++i) {
                // 1) Build the binary label "b_{N-1}...b_1b_0"
                //    where b_k = (i >> k) & 1, and we print k=N-1 down to 0.
                std::string bits;
                bits.reserve(N);
                for(int q = int(N) - 1; q >= 0; --q) {
                    bits.push_back(char('0' + ((i >> q) & 1)));
                }
        
                // 2) Print:    000:(re,im)  001:(re,im)  010:(re,im) ...
                printf("%s:(%.6f,%.6f)  ",
                       bits.c_str(),
                       std::real(psi[i]), std::imag(psi[i]));
            }
            printf("\n");
        }

        virtual void C2_GATE_L(const std::array<Cplx,16>& U4,
                               IdxType q0,
                               IdxType q1)
        {
            using Clock = std::chrono::high_resolution_clock;
        
            auto t_start = Clock::now();
        
            std::array<Cplx,16> U4_eff = U4;
        
            IdxType Dl = bond_dims[q0];
            IdxType Dr = bond_dims[q1+1];
        
            tamm::Tensor<Cplx> M({bond_tis[q0], phys_tis[q0],
                                  phys_tis[q1], bond_tis[q1+1]});
            M.set_dense();
            M.allocate(&ec);
        
            tamm::Scheduler sch{ec};
            sch(M("l","p0","p1","r") =
                  mps_tensors[q0]("l","p0","b") *
                  mps_tensors[q1]("b","p1","r"),
                "merge_two", tamm::ExecutionHW::GPU);
            sch.execute(tamm::ExecutionHW::GPU);
        
            tamm::Tensor<Cplx> G4({ phys_tis[q0], phys_tis[q1],
                                    phys_tis[q0], phys_tis[q1] });
            G4.set_dense();
            G4.allocate(&ec);
        
            for(const auto& blockid : G4.loop_nest()) {
                size_t bs = G4.block_size(blockid);
                std::vector<Cplx> hostbuf(bs);
                auto dims = G4.block_dims(blockid);
                auto offs = G4.block_offsets(blockid);
                size_t c = 0;
                for(size_t p0p = offs[0]; p0p < offs[0] + dims[0]; ++p0p) {
                    for(size_t p1p = offs[1]; p1p < offs[1] + dims[1]; ++p1p) {
                        for(size_t p0 = offs[2]; p0 < offs[2] + dims[2]; ++p0) {
                            for(size_t p1 = offs[3]; p1 < offs[3] + dims[3]; ++p1, ++c) {
                                int row = int(p0p*2 + p1p);
                                int col = int(p0*2 + p1);
                                hostbuf[c] = U4_eff[row*4 + col];
                            }
                        }
                    }
                }
                G4.put(blockid, hostbuf);
            }
        
            tamm::Tensor<Cplx> M2({bond_tis[q0], phys_tis[q0],
                                   phys_tis[q1], bond_tis[q1+1]});
            M2.set_dense();
            M2.allocate(&ec);
        
            {
                tamm::Scheduler sch2{ec};
                sch2(M2("l","p0p","p1p","r") =
                       G4("p0p","p1p","p0","p1") * M("l","p0","p1","r"),
                     "apply_two", tamm::ExecutionHW::GPU);
                sch2.execute(tamm::ExecutionHW::GPU);
            }
        
            M.deallocate();
            G4.deallocate();
        
            Eigen::Index rows = Dl * 2;
            Eigen::Index cols = 2 * Dr;
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> mat(rows, cols);
        
            for(const auto& blockid : M2.loop_nest()) {
                size_t bs = M2.block_size(blockid);
                std::vector<Cplx> hostbuf(bs);
                M2.get(blockid, hostbuf);
                auto dims = M2.block_dims(blockid);
                auto offs = M2.block_offsets(blockid);
                size_t c = 0;
                for(size_t l = offs[0]; l < offs[0] + dims[0]; ++l) {
                    for(size_t p0 = offs[1]; p0 < offs[1] + dims[1]; ++p0) {
                        for(size_t p1 = offs[2]; p1 < offs[2] + dims[2]; ++p1) {
                            for(size_t r = offs[3]; r < offs[3] + dims[3]; ++r, ++c) {
                                mat(l*2 + p0, p1*Dr + r) = hostbuf[c];
                            }
                        }
                    }
                }
            }
        
            Eigen::BDCSVD<decltype(mat)> svd(mat,
                Eigen::ComputeThinU | Eigen::ComputeThinV);
            auto svals = svd.singularValues();
            IdxType chi = std::min<IdxType>(max_bond_dim,
                                             IdxType(svals.size()));
            auto Umat  = svd.matrixU().leftCols(chi);
            auto Sdiag = svals.head(chi).asDiagonal();
            auto Vh    = svd.matrixV().leftCols(chi).adjoint();
        
            bond_dims[q0+1] = chi;
            {
                tamm::IndexSpace is_new{ tamm::range(chi) };
                bond_tis[q0+1] = tamm::TiledIndexSpace(is_new, block_size);
            }
        
            tamm::Tensor<Cplx> Ti_new({
                bond_tis[q0], phys_tis[q0],
                tamm::TiledIndexSpace(tamm::range(chi), block_size)
            });
            Ti_new.set_dense();
            Ti_new.allocate(&ec);
        
            for(const auto& blockid : Ti_new.loop_nest()) {
                size_t bs = Ti_new.block_size(blockid);
                std::vector<Cplx> hostbuf(bs);
                auto dims = Ti_new.block_dims(blockid);
                auto offs = Ti_new.block_offsets(blockid);
                size_t c = 0;
                for(size_t l = offs[0]; l < offs[0] + dims[0]; ++l) {
                    for(size_t p0 = offs[1]; p0 < offs[1] + dims[1]; ++p0) {
                        for(size_t b = offs[2]; b < offs[2] + dims[2]; ++b, ++c) {
                            hostbuf[c] = Umat(l*2 + p0, b);
                        }
                    }
                }
                Ti_new.put(blockid, hostbuf);
            }
        
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> SV = Sdiag * Vh;
        
            tamm::Tensor<Cplx> Tj_new({
                tamm::TiledIndexSpace(tamm::range(chi), block_size),
                phys_tis[q1],
                bond_tis[q1+1]
            });
            Tj_new.set_dense();
            Tj_new.allocate(&ec);
        
            for(const auto& blockid : Tj_new.loop_nest()) {
                size_t bs = Tj_new.block_size(blockid);
                std::vector<Cplx> hostbuf(bs);
                auto dims = Tj_new.block_dims(blockid);
                auto offs = Tj_new.block_offsets(blockid);
                size_t c = 0;
                for(size_t b = offs[0]; b < offs[0] + dims[0]; ++b) {
                    for(size_t p1 = offs[1]; p1 < offs[1] + dims[1]; ++p1) {
                        for(size_t r = offs[2]; r < offs[2] + dims[2]; ++r, ++c) {
                            hostbuf[c] = SV(b, p1*Dr + r);
                        }
                    }
                }
                Tj_new.put(blockid, hostbuf);
            }
        
            mps_tensors[q0].deallocate();
            mps_tensors[q1].deallocate();
            mps_tensors[q0] = std::move(Ti_new);
            mps_tensors[q1] = std::move(Tj_new);
        
            M2.deallocate();
        }

        virtual void C2_GATE_NL(const std::array<Cplx,16>& U4, IdxType q0,IdxType q1)
        {
            // 0) Center MPS at control qubit
            position(q0);
        
            // 1) Gate SVD: 4×4 → Ugate·Sgate·Vhgate
            Eigen::Matrix<Cplx,4,4> Gmat;
            for(int r = 0; r < 4; ++r)
                for(int c = 0; c < 4; ++c)
                    Gmat(r,c) = U4[r*4 + c];
            Eigen::BDCSVD<decltype(Gmat)> gate_svd(
                Gmat, Eigen::ComputeThinU | Eigen::ComputeThinV);
            auto svals   = gate_svd.singularValues();                  // length‐4 vector
            IdxType chiG = std::min<IdxType>(max_bond_dim,
                                              static_cast<IdxType>(svals.size()));
            auto Ugate   = gate_svd.matrixU().leftCols(chiG);          // (4×chiG)
            auto Vhgate  = gate_svd.matrixV().leftCols(chiG).adjoint(); // (chiG×4)
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> SVfull =
                svals.head(chiG).asDiagonal() * Vhgate;
        
            // 2) Merge control site with Ugate + local SVD
            IdxType Dl     = bond_dims[q0];
            IdxType Dr_old = bond_dims[q0+1];
            tamm::TiledIndexSpace mid_ti  { tamm::IndexSpace(tamm::range(chiG)), 1 };
            tamm::TiledIndexSpace r_old_ti{ tamm::IndexSpace(tamm::range(Dr_old)),1 };
        
            // Build Uten[p_out,p_in,a]
            tamm::Tensor<Cplx> Uten({ phys_tis[q0], phys_tis[q0], mid_ti });
            Uten.set_dense(); Uten.allocate(&ec);
            Uten.loop_nest().iterate([&](auto const& idxs){
                auto pout = idxs[0], pin = idxs[1], a = idxs[2];
                Cplx v = Ugate(pout*2 + pin, a);
                Uten.put(idxs, gsl::span<Cplx>(&v,1));
            });
        
            // Contract into M0[l,pout,a,r_old]
            tamm::Tensor<Cplx> M0({ bond_tis[q0], phys_tis[q0], mid_ti, r_old_ti });
            M0.set_dense(); M0.allocate(&ec);
            {
                tamm::Scheduler sch{ec};
                sch(M0("l","pout","a","r") =
                      Uten("pout","pin","a") * mps_tensors[q0]("l","pin","r"),
                    "merge_control", tamm::ExecutionHW::CPU);
                sch.execute(tamm::ExecutionHW::CPU);
            }
            Uten.deallocate();
        
            // Flatten and SVD M0 → U0,S0,V0h
            Eigen::Index rows0 = Dl * 2 * chiG;
            Eigen::Index cols0 = Dr_old;
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> mat0(rows0, cols0);
            M0.loop_nest().iterate([&](auto const& idxs){
                auto l    = idxs[0];
                auto pout = idxs[1];
                auto a    = idxs[2];
                auto r    = idxs[3];
                Cplx v; M0.get(idxs, gsl::span<Cplx>(&v,1));
                mat0(l*2*chiG + pout*chiG + a, r) = v;
            });
            auto svd0   = Eigen::BDCSVD<decltype(mat0)>(
                             mat0,
                             Eigen::ComputeThinU | Eigen::ComputeThinV);
            auto s0     = svd0.singularValues();                      // length‐≤chiG
            IdxType chi0= std::min<IdxType>(max_bond_dim,
                                             static_cast<IdxType>(s0.size()));
            auto U0mat  = svd0.matrixU().leftCols(chi0);               // (rows0×chi0)
            auto V0h    = svd0.matrixV().leftCols(chi0).adjoint();     // (chi0×Dr_old)
        
            // Update left tensor at q0
            bond_dims[q0+1] = chi0;
            {
                tamm::IndexSpace is_new{ tamm::range(chi0) };
                bond_tis[q0+1] = tamm::TiledIndexSpace(is_new,1);
            }
            tamm::Tensor<Cplx> T0new({ bond_tis[q0], phys_tis[q0], bond_tis[q0+1] });
            T0new.set_dense(); T0new.allocate(&ec);
            T0new.loop_nest().iterate([&](auto const& idxs){
                auto l    = idxs[0];
                auto pout = idxs[1];
                auto b    = idxs[2];
                Cplx v = U0mat(l*2*chiG + pout*chiG + b, b);
                T0new.put(idxs, gsl::span<Cplx>(&v,1));
            });
            mps_tensors[q0].deallocate();
            mps_tensors[q0] = std::move(T0new);
            M0.deallocate();
        
            // 3) Build initial propagation bond prop(a,l) = s0(a) * V0h(a,l)
            tamm::Tensor<Cplx> prop({ tamm::TiledIndexSpace{ tamm::IndexSpace{tamm::range(chi0)},1 },
                                      r_old_ti });
            prop.set_dense(); prop.allocate(&ec);
            prop.loop_nest().iterate([&](auto const& idxs){
                auto a = idxs[0], l = idxs[1];
                Cplx v = s0(a) * V0h(a,l);
                prop.put(idxs, gsl::span<Cplx>(&v,1));
            });
            IdxType currChi = chi0;
        
            // Propagate through each intermediate site
            for(IdxType site = q0+1; site < q1; ++site) {
                IdxType Dr_i = bond_dims[site+1];
                tamm::Tensor<Cplx> M1({ tamm::TiledIndexSpace{ tamm::IndexSpace{tamm::range(currChi)},1 },
                                          phys_tis[site],
                                          tamm::TiledIndexSpace{ tamm::IndexSpace{tamm::range(Dr_i)},1 } });
                M1.set_dense(); M1.allocate(&ec);
                {
                    tamm::Scheduler sch{ec};
                    sch(M1("a","p","r") =
                          prop("a","l") * mps_tensors[site]("l","p","r"),
                        "merge_prop", tamm::ExecutionHW::CPU);
                    sch.execute(tamm::ExecutionHW::CPU);
                }
        
                // SVD at site
                Eigen::Index rows1 = currChi * 2;
                Eigen::Index cols1 = Dr_i;
                Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> mat1(rows1, cols1);
                M1.loop_nest().iterate([&](auto const& idxs){
                    auto a = idxs[0], p = idxs[1], r = idxs[2];
                    Cplx v; M1.get(idxs, gsl::span<Cplx>(&v,1));
                    mat1(a*2 + p, r) = v;
                });
                auto svd1   = Eigen::BDCSVD<decltype(mat1)>(
                                  mat1,
                                  Eigen::ComputeThinU | Eigen::ComputeThinV);
                auto s1     = svd1.singularValues();                  // length‐≤currChi*2
                IdxType chi1= std::min<IdxType>(max_bond_dim,
                                                static_cast<IdxType>(s1.size()));
                auto U1mat  = svd1.matrixU().leftCols(chi1);
                auto V1h    = svd1.matrixV().leftCols(chi1).adjoint();
        
                // Write back updated tensor
                bond_dims[site+1] = chi1;
                bond_tis[site+1]  = tamm::TiledIndexSpace{ tamm::IndexSpace{tamm::range(chi1)},1 };
                tamm::Tensor<Cplx> T1new({ bond_tis[site],
                                           phys_tis[site],
                                           bond_tis[site+1] });
                T1new.set_dense(); T1new.allocate(&ec);
                T1new.loop_nest().iterate([&](auto const& idxs){
                    auto a = idxs[0], p = idxs[1], b = idxs[2];
                    Cplx v = U1mat(a*2 + p, b);
                    T1new.put(idxs, gsl::span<Cplx>(&v,1));
                });
                mps_tensors[site].deallocate();
                mps_tensors[site] = std::move(T1new);
                M1.deallocate();
        
                // Rebuild prop = s1(b) * V1h(b,r)
                tamm::Tensor<Cplx> propnew({ tamm::TiledIndexSpace{ tamm::IndexSpace{tamm::range(chi1)},1 },
                                              tamm::TiledIndexSpace{ tamm::IndexSpace{tamm::range(bond_dims[site+1])},1 } });
                propnew.set_dense(); propnew.allocate(&ec);
                propnew.loop_nest().iterate([&](auto const& idxs){
                    auto b = idxs[0], r = idxs[1];
                    Cplx v = s1(b) * V1h(b,r);
                    propnew.put(idxs, gsl::span<Cplx>(&v,1));
                });
                prop.deallocate();
                prop = std::move(propnew);
                currChi = chi1;
            }
        
            // 4a) Propagate into the target site
            tamm::Tensor<Cplx> Tmid({ tamm::TiledIndexSpace{ tamm::IndexSpace{tamm::range(currChi)},1 },
                                      phys_tis[q1],
                                      tamm::TiledIndexSpace{ tamm::IndexSpace{tamm::range(bond_dims[q1+1])},1 } });
            Tmid.set_dense(); Tmid.allocate(&ec);
            {
                tamm::Scheduler sch{ec};
                sch(Tmid("a","pin","r") =
                      prop("a","l") * mps_tensors[q1]("l","pin","r"),
                    "propagate_to_target", tamm::ExecutionHW::CPU);
                sch.execute(tamm::ExecutionHW::CPU);
            }
            mps_tensors[q1].deallocate();
        
            // 4b) Absorb the gate spectrum into Tmid → Tfin
            tamm::Tensor<Cplx> Gten({ tamm::TiledIndexSpace{ tamm::IndexSpace{tamm::range(chiG)},1 },
                                      phys_tis[q1],
                                      phys_tis[q1] });
            Gten.set_dense(); Gten.allocate(&ec);
            Gten.loop_nest().iterate([&](auto const& idxs){
                auto a    = idxs[0];
                auto pout = idxs[1];
                auto pin  = idxs[2];
                Cplx v = SVfull(a, pout*2 + pin);
                Gten.put(idxs, gsl::span<Cplx>(&v,1));
            });
        
            tamm::Tensor<Cplx> Tfin({ bond_tis[q1], phys_tis[q1], bond_tis[q1+1] });
            Tfin.set_dense(); Tfin.allocate(&ec);
            {
                tamm::Scheduler sch{ec};
                sch(Tfin("a","p","r") =
                      Gten("a","p","pin") * Tmid("a","pin","r"),
                    "absorb_gate_spectrum", tamm::ExecutionHW::CPU);
                sch.execute(tamm::ExecutionHW::CPU);
            }
        
            // Cleanup and install final tensor
            Tmid.deallocate();
            prop.deallocate();
            Gten.deallocate();
            mps_tensors[q1] = std::move(Tfin);
        }

        // The 4×4 SWAP matrix in row‐major order, in the computational basis {00,01,10,11}
        static constexpr std::array<Cplx,16> SWAP_U4 = {
            // row 00 → col 00
            Cplx(1,0), Cplx(0,0), Cplx(0,0), Cplx(0,0),
            // row 01 → col 10
            Cplx(0,0), Cplx(0,0), Cplx(1,0), Cplx(0,0),
            // row 10 → col 01
            Cplx(0,0), Cplx(1,0), Cplx(0,0), Cplx(0,0),
            // row 11 → col 11
            Cplx(0,0), Cplx(0,0), Cplx(0,0), Cplx(1,0)
        };
        
        // Apply an arbitrary two‐qubit gate U4 on qubits q0, q1 (not necessarily adjacent)
        void C2_GATE_NL_SWAP(const std::array<Cplx,16>& U4, IdxType q0, IdxType q1) {
            // 1) Determine ordering and (if necessary) permute U4
            bool reversed = (q0 > q1);
            IdxType i = std::min(q0, q1), j = std::max(q0, q1);
            std::array<Cplx,16> U4_eff;
            if(!reversed) {
                U4_eff = U4;
            } else {
                // swap the role of the two qubits in U4:
                // (p0',p1') ↔ (p1',p0'), (p0,p1) ↔ (p1,p0)
                for(int r0=0; r0<2; ++r0)
                for(int r1=0; r1<2; ++r1)
                for(int c0=0; c0<2; ++c0)
                for(int c1=0; c1<2; ++c1) {
                    int src = (r0*2 + r1)*4 + (c0*2 + c1);
                    int dst = (r1*2 + r0)*4 + (c1*2 + c0);
                    U4_eff[dst] = U4[src];
                }
            }
        
            // 2) Bubble qubit 'i' forward until it sits at position j−1
            for(IdxType k = i; k < j-1; ++k) {
                C2_GATE_L(SWAP_U4, k, k+1);
            }
        
            // 3) Apply the target gate on the now‐adjacent pair (j−1, j)
            C2_GATE_L(U4_eff, j-1, j);
        
            // 4) Undo the SWAP cascade
            for(IdxType k = j-1; k > i; --k) {
                C2_GATE_L(SWAP_U4, k-1, k);
            }
        }
        
        // Finally, modify your dispatcher in C2_GATE:
        virtual void C2_GATE(const std::array<Cplx,16>& U4, IdxType q0, IdxType q1) {
            //dump_state("Before C2");
            //position(q0);  // center the MPS on q0 first
            if(std::abs(q0 - q1) == 1) {
                C2_GATE_L(U4, q0, q1);
            } else {
                C2_GATE_NL_SWAP(U4, q0, q1);
            }
            //dump_state("After C2");
        }

        //TODO: Rewrite these canonialization functions using QR

        void right_canonicalize(std::vector<tamm::Tensor<Cplx>>& MPS) {
          for(IdxType i = n_qubits - 1; i > 0; --i) {
            // 0) dimensions
            IdxType Dl_old = bond_dims[i];
            IdxType Dr     = bond_dims[i+1];
            IdxType d      = phys_dims[i];
            IdxType Dl_prev= bond_dims[i-1];
            IdxType d_prev = phys_dims[i-1];
        
            //
            // 1) Bulk-get MPS[i] into Eigen Mmat of shape (Dl_old × (d*Dr))
            //
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Mmat(Dl_old, d * Dr);
            {
              auto& T = MPS[i];
              for(const auto& blockid : T.loop_nest()) {
                // allocate host buffer for entire tile
                const size_t bs = T.block_size(blockid);
                std::vector<Cplx> hostbuf(bs);
                T.get(blockid, hostbuf);
        
                // unravel each element into its (l,p,r) and write into Mmat
                auto dims    = T.block_dims(blockid);
                auto offs    = T.block_offsets(blockid);
                size_t idx = 0;
                for(size_t ll = offs[0]; ll < offs[0] + dims[0]; ++ll)
                for(size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
                for(size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx) {
                  Mmat(ll, pp * Dr + rr) = hostbuf[idx];
                }
              }
            }
        
            //
            // 2) SVD: Mmat = U · S · Vh
            //
            Eigen::BDCSVD<decltype(Mmat)> svd(
              Mmat, Eigen::ComputeThinU | Eigen::ComputeThinV
            );
            auto svals      = svd.singularValues();
            IdxType full    = IdxType(svals.size());
            IdxType chi     = std::min<IdxType>(full, max_bond_dim);
            auto Umat       = svd.matrixU().leftCols(chi);            // (Dl_old × chi)
            auto Sdiag      = svals.head(chi).asDiagonal();           // (chi × chi)
            auto Vh         = svd.matrixV().leftCols(chi).adjoint();  // (chi × (d*Dr))
        
            //
            // 3) Update bond_dims and bond_tis
            //
            bond_dims[i] = chi;
            {
              tamm::IndexSpace is_new{tamm::range(chi)};
              bond_tis[i] = tamm::TiledIndexSpace(is_new, block_size);
            }
        
            //
            // 4) Build new right tensor Tnew via bulk-puts of Vh
            //
            tamm::Tensor<Cplx> Tnew({bond_tis[i], phys_tis[i], bond_tis[i+1]});
            Tnew.set_dense();
            Tnew.allocate(&ec);
            for(const auto& blockid : Tnew.loop_nest()) {
              const size_t bs = Tnew.block_size(blockid);
              std::vector<Cplx> hostbuf(bs);
        
              auto dims = Tnew.block_dims(blockid);
              auto offs = Tnew.block_offsets(blockid);
              size_t idx = 0;
              for(size_t ll = offs[0]; ll < offs[0] + dims[0]; ++ll)
              for(size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
              for(size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx) {
                hostbuf[idx] = Vh(ll, pp * Dr + rr);
              }
        
              Tnew.put(blockid, hostbuf);
            }
        
            //
            // 5) Bulk-get MPS[i-1] into Eigen Mprev of shape ((Dl_prev*d_prev) × Dl_old)
            //
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Mprev(Dl_prev * d_prev, Dl_old);
            {
              auto& Told = MPS[i-1];
              for(const auto& blockid : Told.loop_nest()) {
                const size_t bs = Told.block_size(blockid);
                std::vector<Cplx> hostbuf(bs);
                Told.get(blockid, hostbuf);
        
                auto dims = Told.block_dims(blockid);
                auto offs = Told.block_offsets(blockid);
                size_t idx = 0;
                for(size_t ll = offs[0]; ll < offs[0] + dims[0]; ++ll)
                for(size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
                for(size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx) {
                  Mprev(ll * d_prev + pp, rr) = hostbuf[idx];
                }
              }
            }
        
            //
            // 6) Form Mprev2 = Mprev * (Umat * Sdiag)
            //
            auto US     = Umat * Sdiag;    //   (Dl_old × chi)
            auto Mprev2 = Mprev * US;      // (Dl_prev*d_prev × chi)
        
            //
            // 7) Build new left tensor Tprev via bulk-puts of Mprev2
            //
            tamm::Tensor<Cplx> Tprev({bond_tis[i-1], phys_tis[i-1], bond_tis[i]});
            Tprev.set_dense();
            Tprev.allocate(&ec);
            for(const auto& blockid : Tprev.loop_nest()) {
              const size_t bs = Tprev.block_size(blockid);
              std::vector<Cplx> hostbuf(bs);
        
              auto dims = Tprev.block_dims(blockid);
              auto offs = Tprev.block_offsets(blockid);
              size_t idx = 0;
              for(size_t ll = offs[0]; ll < offs[0] + dims[0]; ++ll)
              for(size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
              for(size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx) {
                hostbuf[idx] = Mprev2(ll * d_prev + pp, rr);
              }
        
              Tprev.put(blockid, hostbuf);
            }
        
            //
            // 8) Swap into MPS and cleanup
            //
            MPS[i].deallocate();
            MPS[i-1].deallocate();
            MPS[i]     = std::move(Tnew);
            MPS[i-1]   = std::move(Tprev);
          }
        }


        void left_canonicalize(std::vector<tamm::Tensor<Cplx>>& MPS) {
          for(IdxType i = 0; i < n_qubits - 1; ++i) {
            // 0) dimensions
            IdxType Dl       = bond_dims[i];
            IdxType Dr_old   = bond_dims[i+1];
            IdxType d        = phys_dims[i];
            IdxType d_next   = phys_dims[i+1];
            IdxType Dr_right = bond_dims[i+2];
        
            //
            // 1) Bulk-get MPS[i] into Eigen Mmat of shape (Dl*d) × Dr_old
            //
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Mmat(Dl * d, Dr_old);
            {
              auto& T = MPS[i];
              for(const auto& blockid : T.loop_nest()) {
                const size_t bs = T.block_size(blockid);
                std::vector<Cplx> hostbuf(bs);
                T.get(blockid, hostbuf);
        
                auto dims = T.block_dims(blockid);
                auto offs = T.block_offsets(blockid);
                size_t idx = 0;
                // unravel tile
                for(size_t ll = offs[0]; ll < offs[0] + dims[0]; ++ll)
                for(size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
                for(size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx) {
                  Mmat(ll * d + pp, rr) = hostbuf[idx];
                }
              }
            }
        
            //
            // 2) SVD: Mmat = U · S · Vh
            //
            Eigen::BDCSVD<decltype(Mmat)> svd(
              Mmat, Eigen::ComputeThinU | Eigen::ComputeThinV
            );
            auto svals      = svd.singularValues();
            IdxType full    = IdxType(svals.size());
            IdxType chi     = std::min<IdxType>(full, max_bond_dim);
            auto Umat       = svd.matrixU().leftCols(chi);            // (Dl*d) × chi
            auto Sdiag      = svals.head(chi).asDiagonal();           // chi × chi
            auto Vh         = svd.matrixV().leftCols(chi).adjoint();  // chi × Dr_old
        
            //
            // 3) Update bond_dims and bond_tis
            //
            bond_dims[i+1] = chi;
            {
              tamm::IndexSpace is_new{tamm::range(chi)};
              bond_tis[i+1] = tamm::TiledIndexSpace(is_new, block_size);
            }
        
            //
            // 4) Build new left tensor Tleft via bulk-puts of U
            //
            tamm::Tensor<Cplx> Tleft({bond_tis[i], phys_tis[i], bond_tis[i+1]});
            Tleft.set_dense();
            Tleft.allocate(&ec);
            for(const auto& blockid : Tleft.loop_nest()) {
              const size_t bs = Tleft.block_size(blockid);
              std::vector<Cplx> hostbuf(bs);
        
              auto dims = Tleft.block_dims(blockid);
              auto offs = Tleft.block_offsets(blockid);
              size_t idx = 0;
              for(size_t ll = offs[0]; ll < offs[0] + dims[0]; ++ll)
              for(size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
              for(size_t bb = offs[2]; bb < offs[2] + dims[2]; ++bb, ++idx) {
                // Umat maps (l*d + p) → b
                hostbuf[idx] = Umat(ll * d + pp, bb);
              }
        
              Tleft.put(blockid, hostbuf);
            }
        
            //
            // 5) Bulk-get MPS[i+1] into Eigen Tnext_mat of shape Dr_old × (d_next * Dr_right)
            //
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Tnext_mat(Dr_old, d_next * Dr_right);
            {
              auto& Tnext_old = MPS[i+1];
              for(const auto& blockid : Tnext_old.loop_nest()) {
                const size_t bs = Tnext_old.block_size(blockid);
                std::vector<Cplx> hostbuf(bs);
                Tnext_old.get(blockid, hostbuf);
        
                auto dims = Tnext_old.block_dims(blockid);
                auto offs = Tnext_old.block_offsets(blockid);
                size_t idx = 0;
                for(size_t bb = offs[0]; bb < offs[0] + dims[0]; ++bb)
                for(size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
                for(size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx) {
                  Tnext_mat(bb, pp * Dr_right + rr) = hostbuf[idx];
                }
              }
            }
        
            //
            // 6) Multiply through the spectrum: M2_eig = Sdiag · Vh  (chi × Dr_old)
            //                  then Tnext_mat2 = M2_eig · Tnext_mat
            //
            auto M2_eig     = Sdiag * Vh;         // chi × Dr_old
            auto Tnext_mat2 = M2_eig * Tnext_mat; // chi × (d_next*Dr_right)
        
            //
            // 7) Build new right tensor Tnext via bulk-puts of Tnext_mat2
            //
            tamm::Tensor<Cplx> Tnext({bond_tis[i+1], phys_tis[i+1], bond_tis[i+2]});
            Tnext.set_dense();
            Tnext.allocate(&ec);
            for(const auto& blockid : Tnext.loop_nest()) {
              const size_t bs = Tnext.block_size(blockid);
              std::vector<Cplx> hostbuf(bs);
        
              auto dims = Tnext.block_dims(blockid);
              auto offs = Tnext.block_offsets(blockid);
              size_t idx = 0;
              for(size_t bb = offs[0]; bb < offs[0] + dims[0]; ++bb)
              for(size_t pp = offs[1]; pp < offs[1] + dims[1]; ++pp)
              for(size_t rr = offs[2]; rr < offs[2] + dims[2]; ++rr, ++idx) {
                hostbuf[idx] = Tnext_mat2(bb, pp * Dr_right + rr);
              }
        
              Tnext.put(blockid, hostbuf);
            }
        
            //
            // 8) Swap in and cleanup
            //
            MPS[i].deallocate();
            MPS[i+1].deallocate();
            MPS[i]   = std::move(Tleft);
            MPS[i+1] = std::move(Tnext);
          }
        }


        virtual void MA_GATE(const IdxType repetition) {
          using Clock = std::chrono::steady_clock;
          using MS    = std::chrono::duration<double, std::milli>;
          std::mt19937_64 rng{ std::random_device{}() };
          std::uniform_real_distribution<double> dist(0.0,1.0);

        
          // 1) Allocate result buffer
          auto t_alloc_start = Clock::now();
          SAFE_FREE_HOST(results);
          SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
          std::memset(results, 0, sizeof(IdxType) * repetition);
          auto t_alloc_end = Clock::now();
          std::cout << "[MA_GATE] Result‐buffer alloc: "
                    << MS(t_alloc_end - t_alloc_start).count() << " ms\n";
        
          // RNG setup (negligible, so not timed)
        
          // 2) Canonicalize
          auto t_can_start = Clock::now();
          left_canonicalize(mps_tensors);
          right_canonicalize(mps_tensors);
          auto t_can_end = Clock::now();
          std::cout << "[MA_GATE] Canonicalization: "
                    << MS(t_can_end - t_can_start).count() << " ms\n";
        
          // 3) Loop over repetitions
          for(IdxType rep = 0; rep < repetition; ++rep) {
            auto t_rep_start = Clock::now();
        
            IdxType packed = 0;
            std::vector<Cplx> env(1, Cplx{1.0,0.0});
        
            // 4) Sweep through each site
            for(IdxType site = 0; site < n_qubits; ++site) {
              auto t_site_start = Clock::now();
        
              auto& T = mps_tensors[site];
              IdxType Dl = bond_dims[site];
              IdxType Dr = bond_dims[site+1];
        
              // Branch accumulators
              std::vector<Cplx> env0(Dr, Cplx{0.0,0.0});
              std::vector<Cplx> env1(Dr, Cplx{0.0,0.0});
        
              // 4a) Bulk‐get + unpack each tile
              auto t_tiles_start = Clock::now();
              for(const auto& blockid : T.loop_nest()) {
                size_t bs = T.block_size(blockid);
                std::vector<Cplx> hostbuf(bs);
                T.get(blockid, hostbuf);
        
                auto dims    = T.block_dims(blockid);
                auto offsets = T.block_offsets(blockid);
        
                // unpack tile
                size_t idx = 0;
                for(size_t jj = 0; jj < dims[0]; ++jj)
                for(size_t kk = 0; kk < dims[1]; ++kk)
                for(size_t ll = 0; ll < dims[2]; ++ll, ++idx) {
                  IdxType l = offsets[0] + jj;
                  IdxType s = offsets[1] + kk;
                  IdxType r = offsets[2] + ll;
                  Cplx v = hostbuf[idx];
                  Cplx prod = env[l] * v;
                  if(s == 0) env0[r] += prod;
                  else       env1[r] += prod;
                }
              }
              auto t_tiles_end = Clock::now();
        
              // 4b) Compute weights & sample
              auto t_sample_start = Clock::now();
              long double w0 = 0, w1 = 0;
              for(IdxType r = 0; r < Dr; ++r) {
                w0 += std::norm(env0[r]);
                w1 += std::norm(env1[r]);
              }
              long double sumw = w0 + w1;
              long double p0   = sumw > 0 ? (w0 / sumw) : 0.0L;
              bool outcome1    = (dist(rng) >= static_cast<double>(p0));
              if(outcome1) packed |= (IdxType(1) << site);
              auto& chosen    = outcome1 ? env1 : env0;
              long double norm_branch = outcome1 ? w1 : w0;
              long double invnorm      = norm_branch > 0
                                        ? (1.0L / std::sqrt(norm_branch))
                                        : 0.0L;
              env.assign(Dr, Cplx{0.0,0.0});
              for(IdxType r = 0; r < Dr; ++r)
                env[r] = chosen[r] * static_cast<Cplx>(invnorm);
              auto t_sample_end = Clock::now();
        
              auto t_site_end = Clock::now();
        
              // Print per‐site breakdown (for first 3 reps only)
              if(rep < 3) {
                std::cout << "[MA_GATE] rep=" << rep
                          << " site=" << site
                          << " tile‐unpack: " << MS(t_tiles_end - t_tiles_start).count() << " ms, "
                          << "sample: "       << MS(t_sample_end - t_sample_start).count() << " ms, "
                          << "total-site: "  << MS(t_site_end - t_site_start).count() << " ms\n";
              }
            }
        
            results[rep] = packed;
            auto t_rep_end = Clock::now();
        
            if(rep < 3) {
              std::cout << "[MA_GATE] rep=" << rep
                        << " total: " << MS(t_rep_end - t_rep_start).count() << " ms\n";
            }
          }
        
          std::cout << "[MA_GATE] done (" << repetition << " reps, rep‐details printed for first 3)\n";
        }

        void local_left_step(IdxType i) {
            // dimensions
            const IdxType Dl = bond_dims[i];
            const IdxType Dr = bond_dims[i+1];
            const IdxType d  = phys_dims[i];
            const IdxType d_next = phys_dims[i+1];
            const IdxType Dr_next = bond_dims[i+2];
        
            // 1) Build matrix M = reshape MPS[i] into (Dl * d) × Dr
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Mmat(Dl * d, Dr);
            auto& Ti = mps_tensors[i];
            Ti.loop_nest().iterate([&](auto const& idxs){
                IdxType l = idxs[0], p = idxs[1], r = idxs[2];
                Cplx v; Ti.get(idxs, gsl::span<Cplx>(&v,1));
                Mmat(l*d + p, r) = v;
            });
        
            // 2) SVD: M = U · S · Vh
            Eigen::BDCSVD<decltype(Mmat)> svd(
              Mmat, Eigen::ComputeThinU | Eigen::ComputeThinV);
            auto Umat  = svd.matrixU();
            auto Sdiag = svd.singularValues().asDiagonal();
            auto Vh    = svd.matrixV().adjoint();
        
            // 3) Truncate (if desired) and update bond dimension
            const IdxType chi = Umat.cols();
            bond_dims[i+1] = chi;
            {
              tamm::IndexSpace is_new{tamm::range(chi)};
              bond_tis[i+1] = tamm::TiledIndexSpace(is_new, 1);
            }
        
            // 4) Write back new MPS[i] ← reshape U into (Dl × d × chi)
            tamm::Tensor<Cplx> Ti_new({bond_tis[i], phys_tis[i], bond_tis[i+1]});
            Ti_new.set_dense(); Ti_new.allocate(&ec);
            Ti_new.loop_nest().iterate([&](auto const& idxs){
                IdxType l = idxs[0], p = idxs[1], b = idxs[2];
                Cplx val = Umat(l*d + p, b);
                Ti_new.put(idxs, gsl::span<Cplx>(&val,1));
            });
        
            // 5) Absorb S·Vh into MPS[i+1]:
            //    build A = S · Vh  (size chi × Dr)
            //    then contract A into MPS[i+1] on its left bond
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> A = Sdiag * Vh;
        
            // flatten old MPS[i+1] into (Dr × d_next*Dr_next)
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Bmat(
              Dr, d_next * Dr_next);
            auto& Tj = mps_tensors[i+1];
            Tj.loop_nest().iterate([&](auto const& idxs){
                IdxType b = idxs[0], p = idxs[1], r = idxs[2];
                Cplx v; Tj.get(idxs, gsl::span<Cplx>(&v,1));
                Bmat(b, p*Dr_next + r) = v;
            });
        
            // form new right tensor matrix: (chi × d_next*Dr_next)
            auto Bnew = A * Bmat;
        
            // write back MPS[i+1]
            tamm::Tensor<Cplx> Tj_new(
              {bond_tis[i+1], phys_tis[i+1], bond_tis[i+2]});
            Tj_new.set_dense(); Tj_new.allocate(&ec);
            Tj_new.loop_nest().iterate([&](auto const& idxs){
                IdxType b = idxs[0], p = idxs[1], r = idxs[2];
                Cplx val = Bnew(b, p*Dr_next + r);
                Tj_new.put(idxs, gsl::span<Cplx>(&val,1));
            });
        
            // replace old tensors
            Ti.deallocate();
            Tj.deallocate();
            mps_tensors[i]   = std::move(Ti_new);
            mps_tensors[i+1] = std::move(Tj_new);
        }
        
        /// Move the orthogonality center to `site`:
        void position(IdxType site) {
            assert(site < n_qubits);
            // 1) fully right‐canonicalize to put center at 0
            right_canonicalize(mps_tensors);
            // 2) sweep center rightward from 0 → site
            for(IdxType i = 0; i < site; ++i){
                local_left_step(i);
            }
            // now the MPS is in mixed‐canonical form with center exactly at `site`
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
            //MPI_Init(&argc,&argv);
            //GA_Initialize();
            //tamm::initialize(argc, argv);
            //tamm::ProcGroup::self_ga_pgroup(true);
            return tamm::ProcGroup::create_world_coll();
        }
    };
} // namespace NWQS
