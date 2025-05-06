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
          max_bond_dim(10),
          pg(init_pg()),
          ec(pg, tamm::DistributionKind::dense, tamm::MemoryManagerKind::ga)
        {
            bond_tis.resize(n_qubits + 1);
            bond_dims.resize(n_qubits + 1);
            for(IdxType i = 0; i <= n_qubits; ++i) {
                IdxType dim = (i == 0 || i == n_qubits ? 1 : max_bond_dim);
                bond_dims[i] = dim; 
                tamm::IndexSpace is(tamm::Range{0, static_cast<tamm::Index>(dim)});
                bond_tis[i] = tamm::TiledIndexSpace(is, 1);
            }
        
            phys_tis.resize(n_qubits);
            phys_dims.resize(n_qubits);
            for(IdxType i = 0; i < n_qubits; ++i) {
                phys_dims[i] = 2;
                tamm::IndexSpace is(tamm::Range{0, static_cast<tamm::Index>(2)});
                phys_tis[i] = tamm::TiledIndexSpace(is, 1);
            }
        
            mps_tensors.reserve(n_qubits);
            for(IdxType i = 0; i < n_qubits; ++i) {
                std::cerr << "[TN_CUDA] allocate T["<<i<<"] dims = ("
                << bond_dims[i] << "," 
                << phys_dims[i] << "," 
                << bond_dims[i+1] << ")\n";
                tamm::Tensor<Cplx> T{bond_tis[i], phys_tis[i], bond_tis[i+1]};
                printf("Created first tensor!");
                T.set_dense();
                T.allocate(&ec);
                printf("Allocated first tensor!");
                T.loop_nest().iterate([&](auto const& idxs) {
                    Cplx v = (idxs[0] == 0 && idxs[1] == 0 && idxs[2] == 0)
                             ? Cplx(1.0,0.0)
                             : Cplx(0.0,0.0);
                    T.put(idxs, gsl::span<Cplx>(&v,1));
                });
                mps_tensors.push_back(std::move(T));
            }

            printf("Constructor finished!");
        }

        // Virtual destructor inherits from QuantumState
        ~TN_CUDA() noexcept override 
        {
            SAFE_FREE_HOST(results);
            //GA_Terminate();
            //MPI_Finalize();
            //tamm::finalize();
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
            printf("In Sim!");
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
                else if (g.op_name == OP::MA)
                {
                    measure_all(g.qubit); 
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
            std::cerr << "[DEBUG] Enter measure_all(repetition=" << repetition << ")\n";
        
            // 1) Prepare results buffer
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
            std::memset(results, 0, sizeof(IdxType) * repetition);
            std::cerr << "[DEBUG] results buffer allocated at " << results << "\n";
        
            // 2) RNG setup
            std::mt19937_64                        rng{std::random_device{}()};
            std::uniform_real_distribution<double> dist(0.0, 1.0);
        
            // 3) Repeat for each shot
            for(IdxType rep = 0; rep < 1; ++rep) {
              std::cerr << "[DEBUG] === rep " << rep << " ===\n";
              IdxType bitstring = 0;
              auto work = mps_tensors;  // copy so we can update in place
        
              // 4) Loop over sites
              for(IdxType site = 0; site < n_qubits; ++site) {
                std::cerr << "[DEBUG] --- site " << site << " ---\n";
        
                // 4a) Build the list of A_k and its conjugate
                std::string l_lbl   = "l"  + std::to_string(site);
                std::string p_lbl   = "p"  + std::to_string(site);
                std::string pp_lbl  = p_lbl + "'";
                std::string r_lbl   = "l"  + std::to_string(site+1);
        
                auto Ak  = work[site](l_lbl, p_lbl,   r_lbl);
                auto AkH = tamm::conj(Ak)(l_lbl, pp_lbl, r_lbl);
        
                std::cerr << "[DEBUG] Built factors Ak@" << Ak.tensor().base_ptr()
                          << " and AkH@" << AkH.tensor().base_ptr() << "\n";
        
                // 4b) Allocate R and a fresh Rtmp
                tamm::Tensor<Cplx> R   ({phys_tis[site], phys_tis[site]});
                tamm::Tensor<Cplx> Rtmp({phys_tis[site], phys_tis[site]});
                R   .set_dense(); R   .allocate(&ec);
                Rtmp.set_dense(); Rtmp.allocate(&ec);
                std::cerr << "[DEBUG] Allocated R@" << R.base_ptr()
                          << " and Rtmp@" << Rtmp.base_ptr() << "\n";
        
                // 4c) Build the full expr = Ak * AkH
                auto expr = Ak * AkH;
        
                // 4d) Contract into Rtmp
                {
                  std::cerr << "[DEBUG] Scheduling contraction into Rtmp\n";
                  tamm::Scheduler sch{ec};
                  sch(
                    Rtmp(pp_lbl, p_lbl) = expr,
                    "compute_rdm_site_into_Rtmp",
                    tamm::ExecutionHW::GPU
                  );
                  sch.execute(tamm::ExecutionHW::GPU);
                }
                std::cerr << "[DEBUG] Completed contraction into Rtmp\n";
        
                // 4e) Copy Rtmp → R
                {
                  std::cerr << "[DEBUG] Copying Rtmp into R\n";
                  tamm::Scheduler sch{ec};
                  sch(
                    R(pp_lbl, p_lbl) = Rtmp(pp_lbl, p_lbl),
                    "move_Rtmp_into_R",
                    tamm::ExecutionHW::GPU
                  );
                  sch.execute(tamm::ExecutionHW::GPU);
                }
                Rtmp.deallocate();
                std::cerr << "[DEBUG] Deallocated Rtmp\n";
        
                // 4f) Extract diagonal, sample
                Cplx v0, v1;
                R.get({0,0}, gsl::span<Cplx>(&v0,1));
                R.get({1,1}, gsl::span<Cplx>(&v1,1));
                double p0 = std::real(v0), p1 = std::real(v1);
                std::cerr << "[DEBUG] p0=" << p0 << ", p1=" << p1 << "\n";
                R.deallocate();
        
                bool outcome1 = (dist(rng) >= p0);
                std::cerr << "[DEBUG] sampled outcome=" << outcome1 << "\n";
                if(outcome1) bitstring |= (IdxType(1) << site);
        
                // 4g) Build & apply projector G
                std::array<Cplx,4> P = {
                  Cplx(outcome1 ? 0.0 : 1.0/std::sqrt(p0), 0.0),
                  Cplx(0.0,0.0),
                  Cplx(0.0,0.0),
                  Cplx(outcome1 ? 1.0/std::sqrt(p1) : 0.0, 0.0)
                };
                tamm::Tensor<Cplx> G({phys_tis[site], phys_tis[site]});
                G.set_dense(); G.allocate(&ec);
                G.loop_nest().iterate([&](auto const& idxs){
                  G.put(idxs, gsl::span<Cplx>(&P[idxs[0]*2 + idxs[1]],1));
                });
                std::cerr << "[DEBUG] Projector G built\n";
        
                {
                  auto& T = work[site];
                  tamm::Tensor<Cplx> Tnew({bond_tis[site],
                                           phys_tis[site],
                                           bond_tis[site+1]});
                  Tnew.set_dense(); Tnew.allocate(&ec);
                  tamm::Scheduler sch2{ec};
                  sch2(
                    Tnew("l","p'","r") = G("p'","p") * T("l","p","r"),
                    "proj_site",
                    tamm::ExecutionHW::CPU
                  );
                  sch2.execute(tamm::ExecutionHW::CPU);
                  T.deallocate();
                  work[site] = std::move(Tnew);
                }
                G.deallocate();
                std::cerr << "[DEBUG] Collapsed MPS tensor " << site << "\n";
        
                results[rep] = bitstring;
                std::cerr << "[DEBUG] partial bitstring=0x"
                          << std::hex << bitstring << std::dec << "\n";
              } // per-site
        
              std::cerr << "[DEBUG] rep " << rep
                        << " complete, result=0x"
                        << std::hex << results[rep] << std::dec << "\n";
            } // per-rep
        
            std::cerr << "[DEBUG] Exit measure_all\n";
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



        void apply_one_qubit(const std::array<Cplx,4>& U, IdxType site) 
        {
            printf("In One qubit gaate!");
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
            // 1) Entry log
            std::cerr << "[DEBUG] Enter apply_two_qubit(q0=" << q0
                      << ", q1=" << q1 << ")" << std::endl;
        
            // 2) Local dims from stored vectors
            IdxType Dl = bond_dims[q0];
            IdxType Dr = bond_dims[q1+1];
            IdxType d  = phys_dims[q0];
            std::cerr << "[DEBUG] Dl (bond_dims[" << q0 << "]) = " << Dl
                      << ", Dr (bond_dims[" << (q1+1) << "]) = " << Dr
                      << ", d (phys_dims[" << q0 << "]) = " << d << std::endl;
        
            // 3) Sanity‐check all IndexSpace extents
            auto check_tis = [&](auto const& tis, const char* name) {
                auto extent = tis.index_space().num_indices();
                std::cerr << "[DEBUG] " << name << " extent = " << extent << std::endl;
            };
            check_tis(bond_tis[q0],      "bond_tis[q0]");
            check_tis(phys_tis[q0],      "phys_tis[q0]");
            check_tis(phys_tis[q1],      "phys_tis[q1]");
            check_tis(bond_tis[q1+1],    "bond_tis[q1+1]");
        
            // 4) Merge tensors M = T[q0] × T[q1]
            {
                std::cerr << "[DEBUG] Constructing M("
                          << Dl << "×" << d << "×" << d << "×" << Dr << ")"
                          << std::endl;
                tamm::Tensor<Cplx> M({bond_tis[q0], phys_tis[q0],
                                      phys_tis[q1], bond_tis[q1+1]});
                std::cerr << "[DEBUG] Calling M.set_dense()" << std::endl;
                M.set_dense();
                std::cerr << "[DEBUG] Allocating M" << std::endl;
                M.allocate(&ec);
                std::cerr << "[DEBUG] Allocation of M succeeded" << std::endl;
        
                std::cerr << "[DEBUG] Submitting merge-two scheduler" << std::endl;
                {
                    tamm::Scheduler sch{ec};
                    sch(M("l","p0","p1","r") =
                           mps_tensors[q0]("l","p0","b") *
                           mps_tensors[q1]("b","p1","r"),
                        "merge_two", tamm::ExecutionHW::GPU);
                    sch.execute(tamm::ExecutionHW::GPU);
                }
                std::cerr << "[DEBUG] merge-two scheduler complete" << std::endl;
        
                // 5) Build two‐qubit gate tensor G4
                std::cerr << "[DEBUG] Constructing G4(2×2×2×2)" << std::endl;
                tamm::Tensor<Cplx> G4({phys_tis[q0], phys_tis[q1],
                                       phys_tis[q0], phys_tis[q1]});
                G4.set_dense();
                std::cerr << "[DEBUG] Allocating G4" << std::endl;
                G4.allocate(&ec);
                std::cerr << "[DEBUG] Filling G4 from U4[]" << std::endl;
                G4.loop_nest().iterate([&](auto const& idxs){
                    int row = idxs[0]*d + idxs[1];
                    int col = idxs[2]*d + idxs[3];
                    Cplx v = U4[row*4 + col];
                    G4.put(idxs, gsl::span<Cplx>(&v,1));
                });
                std::cerr << "[DEBUG] G4 initialization complete" << std::endl;
        
                // 6) Apply G4 to M → M2
                std::cerr << "[DEBUG] Constructing M2(same shape as M)" << std::endl;
                tamm::Tensor<Cplx> M2({bond_tis[q0], phys_tis[q0],
                                       phys_tis[q1], bond_tis[q1+1]});
                M2.set_dense();
                std::cerr << "[DEBUG] Allocating M2" << std::endl;
                M2.allocate(&ec);
                std::cerr << "[DEBUG] Submitting apply-two scheduler" << std::endl;
                {
                    tamm::Scheduler sch{ec};
                    sch(M2("l","p0p","p1p","r") =
                           G4("p0p","p1p","p0","p1") * M("l","p0","p1","r"),
                        "apply_two", tamm::ExecutionHW::GPU);
                    sch.execute(tamm::ExecutionHW::GPU);
                }
                std::cerr << "[DEBUG] apply-two scheduler complete" << std::endl;
        
                // 7) Tear down temporaries M and G4
                std::cerr << "[DEBUG] Deallocating M and G4" << std::endl;
                M.deallocate();
                G4.deallocate();
        
                // 8) Move data from M2 into Eigen for SVD
                std::cerr << "[DEBUG] Copying M2 into Eigen matrix of size ("
                          << (Dl*d) << "×" << (d*Dr) << ")" << std::endl;
                Eigen::Index rows = Dl * d;
                Eigen::Index cols = d  * Dr;
                Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> mat(rows, cols);
                M2.loop_nest().iterate([&](auto const& idxs){
                    IdxType l  = idxs[0], p0 = idxs[1], p1 = idxs[2], r = idxs[3];
                    Cplx v;
                    M2.get(idxs, gsl::span<Cplx>(&v,1));
                    mat(l*d + p0, p1*Dr + r) = v;
                });
                std::cerr << "[DEBUG] Eigen matrix population complete" << std::endl;
        
                // 9) SVD and truncation
                std::cerr << "[DEBUG] Computing SVD" << std::endl;
                Eigen::JacobiSVD<decltype(mat)> svd(mat,
                    Eigen::ComputeThinU | Eigen::ComputeThinV);
                IdxType chi = std::min<IdxType>(max_bond_dim,
                                                 static_cast<IdxType>(svd.singularValues().size()));
                std::cerr << "[DEBUG] Truncation chi = " << chi << std::endl;
                auto U  = svd.matrixU().leftCols(chi);
                auto S  = svd.singularValues().head(chi).asDiagonal();
                auto Vh = svd.matrixV().leftCols(chi).adjoint();
                auto US = U * S;
        
                // 10) Build new MPS tensors Ti_new and Tj_new
                std::cerr << "[DEBUG] Constructing Ti_new("
                          << Dl << "×" << d << "×" << chi << ")" << std::endl;
                tamm::Tensor<Cplx> Ti_new({bond_tis[q0], phys_tis[q0],
                                           tamm::TiledIndexSpace(tamm::range(0,chi),1)});
                Ti_new.set_dense();
                Ti_new.allocate(&ec);
                Ti_new.loop_nest().iterate([&](auto const& idxs){
                    IdxType l  = idxs[0], p0 = idxs[1], b = idxs[2];
                    Cplx v = (b < chi ? US(l*d + p0, b) : Cplx{0,0});
                    Ti_new.put(idxs, gsl::span<Cplx>(&v,1));
                });
                std::cerr << "[DEBUG] Ti_new complete" << std::endl;
        
                std::cerr << "[DEBUG] Constructing Tj_new("
                          << chi << "×" << d << "×" << Dr << ")" << std::endl;
                tamm::Tensor<Cplx> Tj_new({
                    tamm::TiledIndexSpace(tamm::range(0,chi),1),
                    phys_tis[q1],
                    bond_tis[q1+1]
                });
                Tj_new.set_dense();
                Tj_new.allocate(&ec);
                Tj_new.loop_nest().iterate([&](auto const& idxs){
                    IdxType b  = idxs[0], p1 = idxs[1], r = idxs[2];
                    Cplx v = (b < chi ? Vh(b, p1*Dr + r) : Cplx{0,0});
                    Tj_new.put(idxs, gsl::span<Cplx>(&v,1));
                });
                std::cerr << "[DEBUG] Tj_new complete" << std::endl;
        
                // 11) Replace old tensors
                std::cerr << "[DEBUG] Deallocating old mps_tensors[" << q0
                          << "] and [" << q1 << "]" << std::endl;
                mps_tensors[q0].deallocate();
                mps_tensors[q1].deallocate();
                mps_tensors[q0] = std::move(Ti_new);
                mps_tensors[q1] = std::move(Tj_new);
        
                // 12) Final deallocation of M2
                std::cerr << "[DEBUG] Deallocating M2 and exit apply_two_qubit" << std::endl;
                M2.deallocate();
            }
        }

        static tamm::ProcGroup init_pg() {
            int argc = 0; char** argv = nullptr;
            //MPI_Init(&argc,&argv);
            //GA_Initialize();
            tamm::initialize(argc, argv);
            //tamm::ProcGroup::self_ga_pgroup(true);
            return tamm::ProcGroup::create_world_coll();
        }
    };
} // namespace NWQS
