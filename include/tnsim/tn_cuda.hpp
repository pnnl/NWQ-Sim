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
        TN_CUDA(IdxType _n_qubits)
        : QuantumState(SimType::TN),
          n_qubits(_n_qubits),
          max_bond_dim(10),
          pg(init_pg()),
          ec(pg, tamm::DistributionKind::dense, tamm::MemoryManagerKind::ga)
        {
            bond_tis.resize(n_qubits + 1);
            bond_dims.resize(n_qubits + 1);
            for(IdxType i = 0; i <= n_qubits; ++i)
            {
                IdxType dim = 1;
                bond_dims[i] = dim; 
                tamm::IndexSpace is{tamm::range(dim)};
                bond_tis[i] = tamm::TiledIndexSpace(is, 1);
            }
        
            phys_tis.resize(n_qubits);
            phys_dims.resize(n_qubits);
            for(IdxType i = 0; i < n_qubits; ++i)
            {
                phys_dims[i] = 2;
                tamm::IndexSpace is{tamm::range(2)};
                phys_tis[i] = tamm::TiledIndexSpace(is, 1);
            }
        
            mps_tensors.reserve(n_qubits);
            for(IdxType i = 0; i < n_qubits; ++i)
            {
                tamm::Tensor<Cplx> T{bond_tis[i], phys_tis[i], bond_tis[i+1]};
                T.set_dense();
                T.allocate(&ec);
                T.loop_nest().iterate([&](auto const& idxs)
                {
                    Cplx v = (idxs[0] == 0 && idxs[1] == 0 && idxs[2] == 0)
                             ? Cplx(1.0,0.0)
                             : Cplx(0.0,0.0);
                    T.put(idxs, gsl::span<Cplx>(&v,1));
                });
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
            throw std::runtime_error("TN_CUDA does not use RNG seed, not accessible form cutensornet API");
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
            IdxType origional_gates = circuit->num_gates();
            std::vector<SVGate> gates = fuse_circuit_sv(circuit);
            IdxType n_gates = gates.size();
            assert(circuit->num_qubits() == n_qubits);
            double sim_time;
            cpu_timer sim_timer;
            sim_timer.start_timer();
	    
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
            auto start = std::chrono::steady_clock::now();
            int n_gates = gates.size();
            for (int i = 0; i < n_gates; i++)
            {
                auto g = gates[i];

                if (g.op_name == OP::C1)
                {
                    std::array<Cplx,4> U;
                    for(int idx=0; idx<4; ++idx)
                    {
                        U[idx] = Cplx(g.gm_real[idx], g.gm_imag[idx]);
                    }

                    C1_GATE(U, g.qubit);
                }
                else if (g.op_name == OP::C2)
                {
                    std::array<Cplx,16> U4;
                    for(int idx=0; idx<16; ++idx)
                        U4[idx] = Cplx(g.gm_real[idx], g.gm_imag[idx]);

                    C2_GATE(U4, g.ctrl, g.qubit);
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
                }
                else if (g.op_name == OP::EXPECT)
                {
                    ObservableList *o = (ObservableList *)(g.data);
                    EXPECT_GATE(o);
                }
                else
                {
                    std::cout << "Unrecognized gates" << std::endl
                              << OP_NAMES[g.op_name] << std::endl;
                    std::logic_error("Invalid gate type");
                }
            }
        }

        void C1_GATE(const std::array<Cplx,4>& U, IdxType site) 
        {
            printf("In One qubit gaate!");
            dump_state("Before");
            tamm::Tensor<Cplx> G({phys_tis[site], phys_tis[site]});
            G.set_dense();
            G.allocate(&ec);
            G.loop_nest().iterate([&](auto const& idxs) 
            {
                Cplx v = U[idxs[0]*2 + idxs[1]];
                G.put(idxs, gsl::span<Cplx>(&v,1));
            });
        
            auto& T = mps_tensors[site];
            tamm::Tensor<Cplx> Tnew({bond_tis[site], phys_tis[site], bond_tis[site+1]});
            Tnew.set_dense();
            Tnew.allocate(&ec);
        
            tamm::Scheduler sch{ec};
            sch(Tnew("l","p'","r") = G("p'","p") * T("l","p","r"),
                "apply_one_qubit", tamm::ExecutionHW::GPU);
            sch.execute(tamm::ExecutionHW::GPU);
        
            T.deallocate();
            mps_tensors[site] = std::move(Tnew);
            dump_state("Before");
        }
        
        // Inside your TN_CUDA class
        std::array<Cplx,4> flatten_two_qubit(const tamm::Tensor<Cplx>& A,
                                             const tamm::Tensor<Cplx>& B)
        {
            // χ = current bond dimension between qubit 0 and 1
            const IdxType chi = bond_dims[1];
            const IdxType d0  = phys_dims[0];  // == 2
            const IdxType d1  = phys_dims[1];  // == 2
        
            // 1) Copy A(1×2×χ) into host vector Adata[p0*χ + b]
            std::vector<Cplx> Adata(d0 * chi);
            A.loop_nest().iterate([&](auto const& idxs){
                IdxType l  = idxs[0];       // == 0
                IdxType p0 = idxs[1];       // 0 or 1
                IdxType b  = idxs[2];       // 0…χ-1
                Cplx    v;
                A.get(idxs, gsl::span<Cplx>(&v,1));
                Adata[p0 * chi + b] = v;
            });
        
            // 2) Copy B(χ×2×1) into host vector Bdata[b*2 + p1]
            std::vector<Cplx> Bdata(chi * d1);
            B.loop_nest().iterate([&](auto const& idxs){
                IdxType b  = idxs[0];       // 0…χ-1
                IdxType p1 = idxs[1];       // 0 or 1
                IdxType r  = idxs[2];       // == 0
                Cplx    v;
                B.get(idxs, gsl::span<Cplx>(&v,1));
                Bdata[b * d1 + p1] = v;
            });
        
            // 3) Build ψ_{p0,p1} = sum_b A(0,p0,b) * B(b,p1,0)
            std::array<Cplx,4> psi{};
            for(IdxType p0 = 0; p0 < d0; ++p0)
            {
                for(IdxType p1 = 0; p1 < d1; ++p1)
                {
                    Cplx sum = Cplx(0.0,0.0);
                    for(IdxType b = 0; b < chi; ++b)
                    {
                        sum += Adata[p0 * chi + b] * Bdata[b * d1 + p1];
                    }
                    psi[p1 * d1 + p0] = sum;
                }
            }
        
            return psi;
        }
 
        void dump_state(const std::string& tag)
        {
          auto psi = flatten_two_qubit(mps_tensors[0], mps_tensors[1]);
          printf("[DUMP %s] ψ = [", tag.c_str());
          for(int i=0; i<4; ++i)
            printf(" (%.6f,%.6f)", std::real(psi[i]), std::imag(psi[i]));
          printf(" ]\n");
        }

        virtual void C2_GATE_L(const std::array<Cplx,16>& U4, IdxType q0, IdxType q1)
        {
            std::array<Cplx,16> U4_eff = U4;
            IdxType Dl = bond_dims[q0];
            IdxType Dr = bond_dims[q1+1];

            tamm::Tensor<Cplx> M({bond_tis[q0], phys_tis[q0], phys_tis[q1], bond_tis[q1+1]});
            M.set_dense(); M.allocate(&ec);

            tamm::Scheduler sch{ec};
            sch(M("l","p0","p1","r") =
                  mps_tensors[q0]("l","p0","b") *
                  mps_tensors[q1]("b","p1","r"),
                "merge_two", tamm::ExecutionHW::GPU);
            sch.execute(tamm::ExecutionHW::GPU);
        
            tamm::Tensor<Cplx> G4({
                phys_tis[q0],  // axis 0 = p0'
                phys_tis[q1],  // axis 1 = p1'
                phys_tis[q0],  // axis 2 = p0
                phys_tis[q1]   // axis 3 = p1
            });

            G4.set_dense(); G4.allocate(&ec);
            G4.loop_nest().iterate([&](auto const& idxs){
                int p0p = idxs[0];
                int p1p = idxs[1];
                int p0  = idxs[2];
                int p1  = idxs[3];
                // Compute row/col in the flattened 4×4 U4_eff
                int row = p0p * 2 + p1p;  
                int col = p0  * 2 + p1;   
                G4.put(idxs, gsl::span<Cplx>(&U4_eff[row * 4 + col], 1));
            });
        
            tamm::Tensor<Cplx> M2({bond_tis[q0], phys_tis[q0], phys_tis[q1], bond_tis[q1+1]});
            M2.set_dense(); M2.allocate(&ec);
            {
                tamm::Scheduler sch{ec};
                sch(M2("l","p0p","p1p","r") =
                      G4("p0p","p1p","p0","p1") * M("l","p0","p1","r"),
                    "apply_two", tamm::ExecutionHW::GPU);
                sch.execute(tamm::ExecutionHW::GPU);
            }
            M.deallocate(); G4.deallocate();
        
            Eigen::Index rows = Dl * 2;
            Eigen::Index cols = 2  * Dr;
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> mat(rows, cols);
            M2.loop_nest().iterate([&](auto const& idxs){
                IdxType l  = idxs[0];
                IdxType p0 = idxs[1];
                IdxType p1 = idxs[2];
                IdxType r  = idxs[3];
                Cplx v; M2.get(idxs, gsl::span<Cplx>(&v,1));
                mat(l*2 + p0, p1*Dr + r) = v;
            });
        
            Eigen::JacobiSVD<decltype(mat)> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
            auto svals = svd.singularValues();
            IdxType chi = std::min<IdxType>(max_bond_dim, (IdxType)svals.size());
        
            auto Umat  = svd.matrixU().leftCols(chi);
            auto Sdiag = svals.head(chi).asDiagonal();
            auto Vh    = svd.matrixV().leftCols(chi).adjoint();
        
            bond_dims[q0+1] = chi;
            {
                tamm::IndexSpace is_new{tamm::range(chi)};
                bond_tis[q0+1] = tamm::TiledIndexSpace(is_new,1);
            }
        
            tamm::Tensor<Cplx> Ti_new({
                bond_tis[q0], phys_tis[q0], tamm::TiledIndexSpace(tamm::range(chi),1)
            });
            Ti_new.set_dense(); Ti_new.allocate(&ec);
            Ti_new.loop_nest().iterate([&](auto const& idxs){
                IdxType l = idxs[0], p0 = idxs[1], b = idxs[2];
                Cplx val = Umat(l*2 + p0, b);
                Ti_new.put(idxs, gsl::span<Cplx>(&val,1));
            });
        
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> SV = Sdiag * Vh;
            tamm::Tensor<Cplx> Tj_new({
                tamm::TiledIndexSpace(tamm::range(chi),1),
                phys_tis[q1],
                bond_tis[q1+1]
            });
            Tj_new.set_dense(); Tj_new.allocate(&ec);
            Tj_new.loop_nest().iterate([&](auto const& idxs){
                IdxType b = idxs[0], p1 = idxs[1], r = idxs[2];
                Cplx val = SV(b, p1*Dr + r);
                Tj_new.put(idxs, gsl::span<Cplx>(&val,1));
            });
            
            mps_tensors[q0].deallocate();
            mps_tensors[q1].deallocate();
            mps_tensors[q0] = std::move(Ti_new);
            mps_tensors[q1] = std::move(Tj_new);
            M2.deallocate();
        }

        virtual void C2_GATE_NL(const std::array<Cplx,16>& U4, IdxType q0, IdxType q1)
        {
            throw std::runtime_error("Non Local 2 Qubit Gates Not Implemnted!");
        }

        virtual void C2_GATE(const std::array<Cplx,16>& U4, IdxType q0, IdxType q1) 
        {
            //dump_state("Before");

            if (std::abs(q0 - q1) != 1)
            {
                C2_GATE_NL(U4, q0, q1); 
            }
            else
            {
                C2_GATE_L(U4, q0, q1);
            }
            //dump_state("after");
        }


        //TODO: Rewrite these canonialization functions using QR

        void right_canonicalize(std::vector<tamm::Tensor<Cplx>>& MPS)
        {
            for (IdxType i = n_qubits - 1; i > 0; --i)
            {
                IdxType Dl_old = bond_dims[i];
                IdxType Dr     = bond_dims[i+1];
                IdxType d      = phys_dims[i];
        
                Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Mmat(Dl_old, d * Dr);
                MPS[i].loop_nest().iterate([&](auto const& idxs)
                {
                    IdxType l = idxs[0], p = idxs[1], r = idxs[2];
                    Cplx v;
                    MPS[i].get(idxs, gsl::span<Cplx>(&v,1));
                    Mmat(l, p * Dr + r) = v;
                });
        
                Eigen::JacobiSVD<decltype(Mmat)> svd(
                    Mmat, Eigen::ComputeThinU | Eigen::ComputeThinV
                );
                auto svals      = svd.singularValues();
                Index full_rank = static_cast<Index>(svals.size());
                Index chi       = std::min<Index>(static_cast<Index>(max_bond_dim), full_rank);
        
                auto Umat  = svd.matrixU().leftCols(chi);
                auto Sdiag = svals.head(chi).asDiagonal();
                auto Vh    = svd.matrixV().leftCols(chi).adjoint();
        
                bond_dims[i] = chi;
                {
                    tamm::IndexSpace is_new{tamm::range(chi)};
                    bond_tis[i] = tamm::TiledIndexSpace(is_new, 1);
                }
        
                tamm::Tensor<Cplx> Tnew({bond_tis[i], phys_tis[i], bond_tis[i+1]});
                Tnew.set_dense();
                Tnew.allocate(&ec);
                Tnew.loop_nest().iterate([&](auto const& idxs)
                {
                    IdxType l = idxs[0], p = idxs[1], r = idxs[2];
                    Cplx val = Vh(l, p * Dr + r);
                    Tnew.put(idxs, gsl::span<Cplx>(&val,1));
                });
        
                IdxType Dl_prev = bond_dims[i-1];
                IdxType d_prev  = phys_dims[i-1];
                Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Mprev(Dl_prev * d_prev, Dl_old);
                MPS[i-1].loop_nest().iterate([&](auto const& idxs)
                {
                    IdxType l = idxs[0], p = idxs[1], r = idxs[2];
                    Cplx v;
                    MPS[i-1].get(idxs, gsl::span<Cplx>(&v,1));
                    Mprev(l * d_prev + p, r) = v;
                });
        
                auto US     = Umat * Sdiag;
                auto Mprev2 = Mprev * US;
        
                tamm::Tensor<Cplx> Tprev({bond_tis[i-1], phys_tis[i-1], bond_tis[i]});
                Tprev.set_dense();
                Tprev.allocate(&ec);
                Tprev.loop_nest().iterate([&](auto const& idxs)
                {
                    IdxType l = idxs[0], p = idxs[1], r = idxs[2];
                    Cplx val = Mprev2(l * d_prev + p, r);
                    Tprev.put(idxs, gsl::span<Cplx>(&val,1));
                });
        
                MPS[i].deallocate();
                MPS[i-1].deallocate();
                MPS[i]   = std::move(Tnew);
                MPS[i-1] = std::move(Tprev);
            }
        }


        void left_canonicalize(std::vector<tamm::Tensor<Cplx>>& MPS)
        {
            for (IdxType i = 0; i < n_qubits - 1; ++i)
            {
                IdxType Dl      = bond_dims[i];
                IdxType Dr_old  = bond_dims[i+1];
                IdxType d       = phys_dims[i];
        
                Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Mmat(Dl * d, Dr_old);
                MPS[i].loop_nest().iterate([&](auto const& idxs)
                {
                    IdxType l = idxs[0], p = idxs[1], r = idxs[2];
                    Cplx v;
                    MPS[i].get(idxs, gsl::span<Cplx>(&v,1));
                    Mmat(l * d + p, r) = v;
                });
        
                Eigen::JacobiSVD<decltype(Mmat)> svd(
                    Mmat, Eigen::ComputeThinU | Eigen::ComputeThinV
                );
                auto svals      = svd.singularValues();
                Index full_rank = static_cast<Index>(svals.size());
                Index chi       = std::min<Index>(static_cast<Index>(max_bond_dim), full_rank);
        
                auto Umat  = svd.matrixU().leftCols(chi);
                auto Sdiag = svals.head(chi).asDiagonal();
                auto Vh    = svd.matrixV().leftCols(chi).adjoint();
        
                bond_dims[i+1] = chi;
                {
                    tamm::IndexSpace is_new{tamm::range(chi)};
                    bond_tis[i+1] = tamm::TiledIndexSpace(is_new, 1);
                }
        
                tamm::Tensor<Cplx> Tleft({bond_tis[i], phys_tis[i], bond_tis[i+1]});
                Tleft.set_dense();
                Tleft.allocate(&ec);
                Tleft.loop_nest().iterate([&](auto const& idxs)
                {
                    IdxType l = idxs[0], p = idxs[1], b = idxs[2];
                    Cplx v = Umat(l * d + p, b);
                    Tleft.put(idxs, gsl::span<Cplx>(&v,1));
                });
        
                auto M2 = Sdiag * Vh;
        
                auto& Tnext_old = MPS[i+1];
                IdxType d_next   = phys_dims[i+1];
                IdxType Dr_right = bond_dims[i+2];
        
                Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Tnext_mat(
                    Dr_old, d_next * Dr_right
                );
                Tnext_old.loop_nest().iterate([&](auto const& idxs)
                {
                    IdxType b = idxs[0], p = idxs[1], r = idxs[2];
                    Cplx v;
                    Tnext_old.get(idxs, gsl::span<Cplx>(&v,1));
                    Tnext_mat(b, p * Dr_right + r) = v;
                });
        
                auto Tnext_mat2 = M2 * Tnext_mat;
        
                tamm::Tensor<Cplx> Tnext({bond_tis[i+1], phys_tis[i+1], bond_tis[i+2]});
                Tnext.set_dense();
                Tnext.allocate(&ec);
                Tnext.loop_nest().iterate([&](auto const& idxs)
                {
                    IdxType b = idxs[0], p = idxs[1], r = idxs[2];
                    Cplx v = Tnext_mat2(b, p * Dr_right + r);
                    Tnext.put(idxs, gsl::span<Cplx>(&v,1));
                });
        
                MPS[i].deallocate();
                Tnext_old.deallocate();
                MPS[i]     = std::move(Tleft);
                MPS[i+1]   = std::move(Tnext);
            }
        }


        virtual void MA_GATE(const IdxType repetition)
        {
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
            std::memset(results, 0, sizeof(IdxType) * repetition);
        
            std::mt19937_64 rng{ std::random_device{}() };
            std::uniform_real_distribution<double> dist(0.0, 1.0);
        
            left_canonicalize(mps_tensors);
            right_canonicalize(mps_tensors);
        
            for (IdxType rep = 0; rep < repetition; ++rep)
            {
                IdxType packed = 0;
                std::vector<Cplx> env(1, Cplx{1.0, 0.0});
        
                for (IdxType site = 0; site < n_qubits; ++site)
                {
                    auto &T = mps_tensors[site];
                    IdxType Dl = bond_dims[site];
                    IdxType Dr = bond_dims[site + 1];
        
                    std::vector<Cplx> env0(Dr, Cplx{0.0, 0.0});
                    std::vector<Cplx> env1(Dr, Cplx{0.0, 0.0});
        
                    T.loop_nest().iterate([&](auto const &idxs)
                    {
                        IdxType l = idxs[0], s = idxs[1], r = idxs[2];
                        Cplx v;
                        T.get(idxs, gsl::span<Cplx>(&v, 1));
                        Cplx prod = env[l] * v;
                        if (s == 0)
                            env0[r] += prod;
                        else
                            env1[r] += prod;
                    });
        
                    long double w0 = 0, w1 = 0;
                    for (IdxType r = 0; r < Dr; ++r)
                    {
                        w0 += std::norm(env0[r]);
                        w1 += std::norm(env1[r]);
                    }
                    long double sumw = w0 + w1;
        
                    long double p0 = (sumw > 0) ? (w0 / sumw) : 0.0L;
                    bool outcome1 = (dist(rng) >= static_cast<double>(p0));
                    if (outcome1)
                        packed |= (IdxType(1) << site);
        
                    auto &chosen = outcome1 ? env1 : env0;
                    long double norm_branch = outcome1 ? w1 : w0;
                    long double invnorm = (norm_branch > 0)
                                             ? (1.0L / std::sqrt(norm_branch))
                                             : 0.0L;
        
                    env.assign(Dr, Cplx{0.0, 0.0});
                    for (IdxType r = 0; r < Dr; ++r)
                        env[r] = chosen[r] * static_cast<Cplx>(invnorm);
                }
        
                results[rep] = packed;
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
            //MPI_Init(&argc,&argv);
            //GA_Initialize();
            //tamm::initialize(argc, argv);
            //tamm::ProcGroup::self_ga_pgroup(true);
            return tamm::ProcGroup::create_world_coll();
        }
    };
} // namespace NWQS
