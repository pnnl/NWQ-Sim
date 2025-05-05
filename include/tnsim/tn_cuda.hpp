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

namespace NWQSim
{
    class TN_CUDA : public QuantumState
    {
    public:
        TN_CUDA(IdxType _n_qubits)
        : QuantumState(SimType::TN),
          n_qubits(_n_qubits)
        {
            int64_t bond_dim = 10;
            MPI_Init(&argc,&argv);
            GA_Initialize();

            tamm::ProcGroup::self_ga_pgroup(true);
            auto pg = tamm::ProcGroup::create_world_coll();
            
            // Use the LOCAL memory manager, so that TAMM will use RMM+GPU
            tamm::ExecutionContext ec{pg,
                            tamm::DistributionKind::nw,
                            tamm::MemoryManagerKind::ga};

            # Setting up MPS
            bont_tis.resize(n_qubits+1);
            for(IdxType i = 0l i <=n_qubits ? 1 : max_bond_dims)
            {
                IdxType dim = (i==0 || i==n_qubits ? 1 : max_bond_dim);
                tamm::IndexSpace is ({0, (tamm::Index)dim - 1});
                bond_tis[i] = tamm::TiledIndexSpace(is, 1);
            }

            phys_tis.resize(n_qubits);
            for(IdxType i=0; i < n_qubits; ++i) 
            {
                tamm::IndexSpace is ({0,1});
                phys_tis[i] = tamm::TiledIndexSpace(is, 1);
            }

            # Initialize to all zeros

            mps_tensors.reserve(n_qubits);
            for(IdxType i=0; i < n_qubits; ++i)
            {
                tamm::Tensor<ValType> T({bond_tis[i], phys_tis[i], bond_tis[i+1]});
                T.allocate(ec.get());
                T.loop_nest().iterate([&](auto const& idxs)
                {
                    ValType v = (idxs[0]==0 && idxs[1]==0 && idxs[2]==0) ? ValType(1.0) : ValType(0.0);
                    T.put(idxs, gsl::span<ValType>($v,1));
                });
                mps_tensors.push_back(std::move(T));
            }
        }

        // Virtual destructor inherits from QuantumState
        ~TN_CUDA() override 
        {
            SAFE_FREE_HOST(results);
            GA_Terminate();
            MPI_Finalize();
        }

        void reset_state() override
        {
            for(IdxType i = 0; i < n_qubits; ++i) {
                auto& T = mps_tensors[i];
                T.loop_nest().iterate([&](auto const& idxs) {
                    // idxs = { left_bond, phys, right_bond }
                    ValType v = (idxs[0] == 0 && idxs[1] == 0 && idxs[2] == 0)
                               ? ValType(1.0)
                               : ValType(0.0);
                    T.put(idxs, gsl::span<ValType>(&v, 1));
              });        
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
                    std::array<ValType,4> U;
                    for(int idx=0; idx<4; ++idx)
                    {
                        U[idx] = ValType(g.gm_real[idx], g.gm_imag[idx]);
                    }

                    apply_one_qubit(U, g.qubit);
                }

                else if(g.op_name == OP::C2)
                {
                    std::array<ValType,16> U4;
                    for(int idx=0; idx<16; ++idx)
                        U4[idx] = ValType(g.gm_real[idx], g.gm_imag[idx]);

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

        IdxType* measure_all(IdxType repetition) override
        {
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
            memset(results, 0, sizeof(IdxType) * repetition);
            std::mt19937_64 rng{ std::random_device{}() };
            std::uniform_real_distribution<double> dist(0.0,1.0);
            for(IdxType rep=0;rep<repetition;++rep){
                IdxType bitstring=0;
                std::vector<tamm::Tensor<ValType>> work=mps_tensors;
                for(IdxType site=0;site<n_qubits;++site){
                    tamm::Tensor<ValType> R({phys_tis[site],phys_tis[site]});
                    R.allocate(ec.get());
                    {
                        tamm::Scheduler sch{*ec};
                        auto expr=work[0]("l0","p0","l1")*tamm::conj(work[0])("l0","p0'","l1");
                        for(IdxType k=1;k<n_qubits;++k){
                            expr=expr*work[k]("l"+std::to_string(k),"p"+std::to_string(k),"l"+std::to_string(k+1))
                                      *tamm::conj(work[k])("l"+std::to_string(k),"p"+std::to_string(k)+"'","l"+std::to_string(k+1));
                        }
                        sch(R("p'","p")=expr,"compute_rdm_site",tamm::ExecutionHW::GPU);
                        sch.execute(tamm::ExecutionHW::GPU);
                    }
                    double p0,p1;
                    ValType v0,v1;
                    R.get({0,0},gsl::span<ValType>(&v0,1));
                    R.get({1,1},gsl::span<ValType>(&v1,1));
                    p0=std::real(v0);
                    p1=std::real(v1);
                    R.deallocate();
                    double u=dist(rng);
                    int b=(u<p0?0:1);
                    if(b==1)bitstring|=(IdxType(1)<<site);
                    std::array<ValType,4> P={
                        ValType(b==0?1.0/std::sqrt(p0):0.0,0.0),
                        ValType(0.0,0.0),
                        ValType(0.0,0.0),
                        ValType(b==1?1.0/std::sqrt(p1):0.0,0.0)
                    };
                    {
                        tamm::Tensor<ValType> G({phys_tis[site],phys_tis[site]});
                        G.allocate(ec.get());
                        G.loop_nest().iterate([&](auto const& idxs){
                            ValType vv=P[idxs[0]*2+idxs[1]];
                            G.put(idxs,gsl::span<ValType>(&vv,1));
                        });
                        auto &T=work[site];
                        tamm::Tensor<ValType> Tnew({bond_tis[site],phys_tis[site],bond_tis[site+1]});
                        Tnew.allocate(ec.get());
                        tamm::Scheduler sch2{*ec};
                        sch2(Tnew("l","p'","r")=G("p'","p")*T("l","p","r"),"proj_site",tamm::ExecutionHW::GPU);
                        sch2.execute(tamm::ExecutionHW::GPU);
                        T.deallocate();
                        work[site]=std::move(Tnew);
                        G.deallocate();
                    }
                }
                results[rep]=bitstring;
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

        std::unique_ptr<tamm::ExecutionContext> ec;
        std::vector<IdxType> bond_dims;
        std::vector<IdxType> phys_dims;
        std::vector<tamm::TiledIndexSpace> bond_tis;
        std::vector<tamm::TiledIndexSpace> phys_tis;
        std::vector<tamm::Tensor<ValType>> mps_tensors;
        IdxType* result = nullptr;


        void apply_one_qubit(const std::array<ValType,4>& U, IdxType site)
        {
            tamm::Tensor<ValType> G({phys_tis[site], phys_tis[site]});
            G.alloate(ec.get());
            G.loop_nest().iterate([&](auto const& idxs)
            {
                ValType v = U[idxs[0]*2 + idxs[1]];
                G.put(idxs, gsl::span<ValType>(&v,1));
            });

            auto& T = mps_tensors[site];
            tamm::Tensor<ValType> Tnew({ bond_tis[site], phys_tis[site], bond_tis[site+1] });
            Tnew.allocate(ec.get());

            tamm::Scheduler sch{*ec};
            sch( Tnew("l", "p'","r") = G("p'","p") * T("l","p","r"),
                    "apply_one_qubit", tamm::ExecutionHW::GPU);
            sch.execute(tamm::ExecutionHW::GPU);

            T.deallocate();
            mps_tensors[site] = std::move(Tnew);
        }

        void apply_two_qubit(const std::array<ValType,16> U4, IdxType q0, IdxType q1)
        {
            IdxType Dl = bond_dims[i];
            IdxType Dr = bond_dims[j+1];
            IdxType d  = phys_dims[i];  // = phys_dims[j] = 2
    
            // 1) Merge T_i and T_j into M[l,p0,p1,r]
            tamm::Tensor<ValType> M({ bond_tis[i],
                                      phys_tis[i],
                                      phys_tis[j],
                                      bond_tis[j+1] });
            M.allocate(ec.get());
            {
              tamm::Scheduler sch{*ec};
              sch( M("l","p0","p1","r") =
                   mps_tensors[i]("l","p0","b") *
                   mps_tensors[j]("b","p1","r"),
                   "merge_two", tamm::ExecutionHW::GPU );
              sch.execute(tamm::ExecutionHW::GPU);
            }

                    tamm::Tensor<ValType> G4({ phys_tis[i],
                                   phys_tis[j],
                                   phys_tis[i],
                                   phys_tis[j] });
        G4.allocate(ec.get());
        G4.loop_nest().iterate([&](auto const& idxs){
          // idxs = { p0', p1', p0, p1 }
          int row = idxs[0]*d + idxs[1];
          int col = idxs[2]*d + idxs[3];
          ValType v = U4[row*4 + col];
          G4.put(idxs, gsl::span<ValType>(&v,1));
        });

        // 3) Apply G4: produce M2[l,p0',p1',r]
        tamm::Tensor<ValType> M2({ bond_tis[i],
                                   phys_tis[i],
                                   phys_tis[j],
                                   bond_tis[j+1] });
        M2.allocate(ec.get());
        {
          tamm::Scheduler sch{*ec};
          sch( M2("l","p0p","p1p","r") =
               G4("p0p","p1p","p0","p1") *
               M ("l","p0","p1","r"),
               "apply_two", tamm::ExecutionHW::GPU );
          sch.execute(tamm::ExecutionHW::GPU);
        }

        M .deallocate();
        G4.deallocate();

        Eigen::Index rows = Dl * d;
        Eigen::Index cols = d  * Dr;
        Eigen::Matrix<ValType, Eigen::Dynamic, Eigen::Dynamic> mat(rows, cols);

        M2.loop_nest().iterate([&](auto const& idxs){
          auto l  = idxs[0], p0 = idxs[1], p1 = idxs[2], r = idxs[3];
          ValType v;
          M2.get(idxs, gsl::span<ValType>(&v,1));
          mat(l*d + p0, p1*Dr + r) = v;
        });

        Eigen::JacobiSVD<decltype(mat)> svd(
          mat, Eigen::ComputeThinU | Eigen::ComputeThinV
        );
        IdxType chi = std::min<IdxType>(
          max_bond_dim,
          static_cast<IdxType>(svd.singularValues().size())
        );

        auto U  = svd.matrixU().leftCols(chi);
        auto S  = svd.singularValues().head(chi).asDiagonal();
        auto Vh = svd.matrixV().leftCols(chi).adjoint();

        Eigen::Matrix<ValType, Eigen::Dynamic, Eigen::Dynamic> US = U * S;

        tamm::Tensor<ValType> Ti_new({ bond_tis[i],
                                       phys_tis[i],
                                       bond_tis[i+1] });
        Ti_new.allocate(ec.get());
        Ti_new.loop_nest().iterate([&](auto const& idxs){
          auto l  = idxs[0], p0 = idxs[1], b = idxs[2];
          ValType v = (b < chi
                       ? US(l*d + p0, b)
                       : ValType(0.0));
          Ti_new.put(idxs, gsl::span<ValType>(&v,1));
        });

        tamm::Tensor<ValType> Tj_new({ bond_tis[i+1],
                                       phys_tis[j],
                                       bond_tis[j+1] });
        Tj_new.allocate(ec.get());
        Tj_new.loop_nest().iterate([&](auto const& idxs){
          auto b  = idxs[0], p1 = idxs[1], r = idxs[2];
          ValType v = (b < chi
                       ? Vh(b, p1*Dr + r)
                       : ValType(0.0));
          Tj_new.put(idxs, gsl::span<ValType>(&v,1));
        });

        mps_tensors[i].deallocate();
        mps_tensors[j].deallocate();
        mps_tensors[i] = std::move(Ti_new);
        mps_tensors[j] = std::move(Tj_new);

        M2.deallocate();
    };

} // namespace NWQS
