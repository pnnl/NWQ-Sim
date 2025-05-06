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

    class TN_CUDA : public QuantumState
    {
    public:
        TN_CUDA(IdxType _n_qubits)
        : QuantumState(SimType::TN),
          n_qubits(_n_qubits),
          max_bond_dim(10)
        {
            int tmp_argc = 0;
            char** tmp_argv = nullptr;
            MPI_Init(&tmp_argc, &tmp_argv);
            GA_Initialize();
        
            tamm::ProcGroup::self_ga_pgroup(true);
            auto pg = tamm::ProcGroup::create_world_coll();
        
            ec = tamm::ExecutionContext{pg,
                                        tamm::DistributionKind::nw,
                                        tamm::MemoryManagerKind::ga};
        
            bond_tis.resize(n_qubits + 1);
            for(IdxType i = 0; i <= n_qubits; ++i) {
                IdxType dim = (i == 0 || i == n_qubits ? 1 : max_bond_dim);
                tamm::IndexSpace is({0, (tamm::Index)dim - 1});
                bond_tis[i] = tamm::TiledIndexSpace(is, 1);
            }
        
            phys_tis.resize(n_qubits);
            for(IdxType i = 0; i < n_qubits; ++i) {
                tamm::IndexSpace is({0, 1});
                phys_tis[i] = tamm::TiledIndexSpace(is, 1);
            }
        
            mps_tensors.reserve(n_qubits);
            for(IdxType i = 0; i < n_qubits; ++i) {
                tamm::Tensor<Cplx> T({bond_tis[i], phys_tis[i], bond_tis[i+1]});
                T.allocate(&ec);
                T.loop_nest().iterate([&](auto const& idxs) {
                    Cplx v = (idxs[0] == 0 && idxs[1] == 0 && idxs[2] == 0)
                             ? Cplx(1.0,0.0)
                             : Cplx(0.0,0.0);
                    T.put(idxs, gsl::span<Cplx>(&v,1));
                });
                mps_tensors.push_back(std::move(T));
            }
        }

        // Virtual destructor inherits from QuantumState
        ~TN_CUDA() noexcept override 
        {
            SAFE_FREE_HOST(results);
            GA_Terminate();
            MPI_Finalize();
        }

        void reset_state() override {
            for(IdxType i = 0; i < n_qubits; ++i) {
                auto& T = mps_tensors[i];
                T.loop_nest().iterate([&](auto const& idxs) {
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
            assert(circuit->num_qubits() == n_qubits);
            // fuse the circuit into a list of tensor-apply gates
            auto gates = fuse_circuit_sv(circuit);

            for(const auto& g : gates)
            {
                if(g.op_name == OP::C1)
                {
                    std::array<Cplx,4> U;
                    for(int idx=0; idx<4; ++idx)
                    {
                        U[idx] = Cplx(g.gm_real[idx], g.gm_imag[idx]);
                    }

                    apply_one_qubit(U, g.qubit);
                }

                else if(g.op_name == OP::C2)
                {
                    std::array<Cplx,16> U4;
                    for(int idx=0; idx<16; ++idx)
                        U4[idx] = Cplx(g.gm_real[idx], g.gm_imag[idx]);

                    apply_two_qubit(U4, g.qubit, g.ctrl);
                }
                else
                {
                    throw std::logic_error("TN_CUDA: unrecognized gate type");
                }
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

        IdxType* measure_all(IdxType repetition) override {
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
            memset(results, 0, sizeof(IdxType) * repetition);
        
            std::mt19937_64                        rng{std::random_device{}()};
            std::uniform_real_distribution<double> dist(0.0, 1.0);
        
            for(IdxType rep = 0; rep < repetition; ++rep) {
                IdxType bitstring = 0;
                auto work = mps_tensors;
        
                for(IdxType site = 0; site < n_qubits; ++site) {
                    std::vector<tamm::LabeledTensor<Cplx>> factors;
                    factors.reserve(n_qubits * 2);
        
                    for(IdxType k = 0; k < n_qubits; ++k) {
                        auto Ak  = work[k](
                            "l" + std::to_string(k),
                            "p" + std::to_string(k),
                            "l" + std::to_string(k+1)
                        );
                        auto AkH = tamm::conj(Ak)(
                            "l" + std::to_string(k),
                            "p" + std::to_string(k) + "'",
                            "l" + std::to_string(k+1)
                        );
                        factors.push_back(Ak);
                        factors.push_back(AkH);
                    }
        
                    tamm::Tensor<Cplx> R({phys_tis[site], phys_tis[site]});
                    R.allocate(&ec);
        
                    auto expr = factors[0];
                    for(size_t i = 1; i < factors.size(); ++i) {
                        expr = expr * factors[i];
                    }
        
                    {
                        tamm::Scheduler sch{ec};
                        sch(
                          R("p'","p") = expr,
                          "compute_rdm_site",
                          tamm::ExecutionHW::GPU
                        );
                        sch.execute(tamm::ExecutionHW::GPU);
                    }
        
                    Cplx v0, v1;
                    R.get({0,0}, gsl::span<Cplx>(&v0,1));
                    R.get({1,1}, gsl::span<Cplx>(&v1,1));
                    double p0 = std::real(v0), p1 = std::real(v1);
                    R.deallocate();
        
                    bool outcome1 = dist(rng) >= p0;
                    if(outcome1) bitstring |= (IdxType(1) << site);
        
                    std::array<Cplx,4> P = {
                      Cplx(outcome1 ? 0.0 : 1.0/std::sqrt(p0), 0.0),
                      Cplx(0.0, 0.0),
                      Cplx(0.0, 0.0),
                      Cplx(outcome1 ? 1.0/std::sqrt(p1) : 0.0, 0.0)
                    };
        
                    tamm::Tensor<Cplx> G({phys_tis[site],phys_tis[site]});
                    G.allocate(&ec);
                    G.loop_nest().iterate([&](auto const& idxs){
                        Cplx val = P[idxs[0]*2 + idxs[1]];
                        G.put(idxs, gsl::span<Cplx>(&val,1));
                    });
        
                    auto& T = work[site];
                    tamm::Tensor<Cplx> Tnew({bond_tis[site],phys_tis[site],bond_tis[site+1]});
                    Tnew.allocate(&ec);
                    {
                        tamm::Scheduler sch2{ec};
                        sch2(
                          Tnew("l","p'","r") = G("p'","p") * T("l","p","r"),
                          "proj_site",
                          tamm::ExecutionHW::GPU
                        );
                        sch2.execute(tamm::ExecutionHW::GPU);
                    }
        
                    T.deallocate();
                    work[site] = std::move(Tnew);
                    G.deallocate();
        
                    results[rep] = bitstring;
                }
            }
        
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

        tamm::ExecutionContext ec;
        std::vector<IdxType> bond_dims;
        std::vector<IdxType> phys_dims;
        std::vector<tamm::TiledIndexSpace> bond_tis;
        std::vector<tamm::TiledIndexSpace> phys_tis;
        std::vector<tamm::Tensor<Cplx>> mps_tensors;
        IdxType* result = nullptr;



        void apply_one_qubit(const std::array<Cplx,4>& U, IdxType site) 
        {
            tamm::Tensor<Cplx> G({phys_tis[site], phys_tis[site]});
            G.allocate(&ec);
            G.loop_nest().iterate([&](auto const& idxs) 
            {
                Cplx v = U[idxs[0]*2 + idxs[1]];
                G.put(idxs, gsl::span<Cplx>(&v,1));
            });
        
            auto& T = mps_tensors[site];
            tamm::Tensor<Cplx> Tnew({bond_tis[site], phys_tis[site], bond_tis[site+1]});
            Tnew.allocate(&ec);
        
            tamm::Scheduler sch{ec};
            sch(Tnew("l","p'","r") = G("p'","p") * T("l","p","r"),
                "apply_one_qubit", tamm::ExecutionHW::GPU);
            sch.execute(tamm::ExecutionHW::GPU);
        
            T.deallocate();
            mps_tensors[site] = std::move(Tnew);
        }

        void apply_two_qubit(const std::array<Cplx,16>& U4, IdxType q0, IdxType q1) 
        {
            IdxType Dl = bond_dims[q0];
            IdxType Dr = bond_dims[q1+1];
            IdxType d  = phys_dims[q0];
        
            tamm::Tensor<Cplx> M({bond_tis[q0], phys_tis[q0], phys_tis[q1], bond_tis[q1+1]});
            M.allocate(&ec);
            {
                tamm::Scheduler sch{ec};
                sch(M("l","p0","p1","r") = mps_tensors[q0]("l","p0","b")
                                         * mps_tensors[q1]("b","p1","r"),
                    "merge_two", tamm::ExecutionHW::GPU);
                sch.execute(tamm::ExecutionHW::GPU);
            }
        
            tamm::Tensor<Cplx> G4({phys_tis[q0], phys_tis[q1], phys_tis[q0], phys_tis[q1]});
            G4.allocate(&ec);
            G4.loop_nest().iterate([&](auto const& idxs)
            {
                int row = idxs[0]*d + idxs[1];
                int col = idxs[2]*d + idxs[3];
                Cplx v = U4[row*4 + col];
                G4.put(idxs, gsl::span<Cplx>(&v,1));
            });
        
            tamm::Tensor<Cplx> M2({bond_tis[q0], phys_tis[q0], phys_tis[q1], bond_tis[q1+1]});
            M2.allocate(&ec);
            {
                tamm::Scheduler sch{ec};
                sch(M2("l","p0p","p1p","r") = G4("p0p","p1p","p0","p1") * M("l","p0","p1","r"),
                    "apply_two", tamm::ExecutionHW::GPU);
                sch.execute(tamm::ExecutionHW::GPU);
            }
        
            M.deallocate();
            G4.deallocate();
        
            Eigen::Index rows = Dl * d;
            Eigen::Index cols = d  * Dr;
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> mat(rows, cols);
            M2.loop_nest().iterate([&](auto const& idxs)
            {
                auto l  = idxs[0], p0 = idxs[1], p1 = idxs[2], r = idxs[3];
                Cplx v;
                M2.get(idxs, gsl::span<Cplx>(&v,1));
                mat(l*d + p0, p1*Dr + r) = v;
            });
        
            Eigen::JacobiSVD<decltype(mat)> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
            IdxType chi = std::min<IdxType>(max_bond_dim,
                                             static_cast<IdxType>(svd.singularValues().size()));
        
            auto U  = svd.matrixU().leftCols(chi);
            auto S  = svd.singularValues().head(chi).asDiagonal();
            auto Vh = svd.matrixV().leftCols(chi).adjoint();
            auto US = U * S;
        
            tamm::Tensor<Cplx> Ti_new({bond_tis[q0], phys_tis[q0], bond_tis[q0+1]});
            Ti_new.allocate(&ec);
            Ti_new.loop_nest().iterate([&](auto const& idxs){
                auto l  = idxs[0], p0 = idxs[1], b = idxs[2];
                Cplx v = (b < chi ? US(l*d + p0, b) : Cplx(0.0,0.0));
                Ti_new.put(idxs, gsl::span<Cplx>(&v,1));
            });
        
            tamm::Tensor<Cplx> Tj_new({bond_tis[q0+1], phys_tis[q1], bond_tis[q1+1]});
            Tj_new.allocate(&ec);
            Tj_new.loop_nest().iterate([&](auto const& idxs){
                auto b  = idxs[0], p1 = idxs[1], r = idxs[2];
                Cplx v = (b < chi ? Vh(b, p1*Dr + r) : Cplx(0.0,0.0));
                Tj_new.put(idxs, gsl::span<Cplx>(&v,1));
            });
        
            mps_tensors[q0].deallocate();
            mps_tensors[q1].deallocate();
            mps_tensors[q0] = std::move(Ti_new);
            mps_tensors[q1] = std::move(Tj_new);
        
            M2.deallocate();
        }

} // namespace NWQS
