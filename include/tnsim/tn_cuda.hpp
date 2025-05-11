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
                tamm::IndexSpace is{tamm::range(dim)};
                bond_tis[i] = tamm::TiledIndexSpace(is, 1);
            }
        
            phys_tis.resize(n_qubits);
            phys_dims.resize(n_qubits);
            for(IdxType i = 0; i < n_qubits; ++i) {
                phys_dims[i] = 2;
                tamm::IndexSpace is{tamm::range(2)};
                phys_tis[i] = tamm::TiledIndexSpace(is, 1);
            }
        
            mps_tensors.reserve(n_qubits);
            for(IdxType i = 0; i < n_qubits; ++i) {
                tamm::Tensor<Cplx> T{bond_tis[i], phys_tis[i], bond_tis[i+1]};
                T.set_dense();
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

         //Virtual destructor inherits from QuantumState
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

                    apply_two_qubit(U4, g.ctrl, g.qubit);
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
            // 1) Prepare results buffer
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
            std::memset(results, 0, sizeof(IdxType) * repetition);
        
            // 2) RNG setup
            std::mt19937_64                        rng{std::random_device{}()};
            std::uniform_real_distribution<double> dist(0.0, 1.0);
        
            // 3) One‐time left‐canonical sweep on the live MPS
            left_canonicalize(mps_tensors);
        
            // Temp host buffers
            std::vector<Cplx> env, tmp0, tmp1;
        
            // 4) Repeat for each shot
            for(IdxType rep = 0; rep < repetition; ++rep) {
                IdxType packed = 0;
        
                // 4a) Deep‐copy the canonical MPS into work
                std::vector<tamm::Tensor<Cplx>> work;
                work.reserve(n_qubits);
                for(IdxType site = 0; site < n_qubits; ++site) {
                    auto& T0 = mps_tensors[site];
                    tamm::Tensor<Cplx> Tcopy({bond_tis[site],
                                              phys_tis[site],
                                              bond_tis[site+1]});
                    Tcopy.set_dense();
                    Tcopy.allocate(&ec);
                    tamm::Scheduler sch_copy{ec};
                    sch_copy(Tcopy("l","p","r") = T0("l","p","r"),
                             "mps_deep_copy", tamm::ExecutionHW::GPU);
                    sch_copy.execute(tamm::ExecutionHW::GPU);
                    work.push_back(std::move(Tcopy));
                }
        
                // 4b) Initialize the left‐environment
                env.assign(1, Cplx(1.0, 0.0));
        
                // 5) Sample qubits in right→left order
                for(IdxType j = 1; j <= n_qubits; ++j) {
                    IdxType site = n_qubits - j;               // rightmost first
                    IdxType Dr   = bond_dims[site+1];
        
                    tmp0.assign(Dr, Cplx(0.0,0.0));
                    tmp1.assign(Dr, Cplx(0.0,0.0));
        
                    // 5a) Contract env into the two physical slices
                    auto& T = work[site];
                    T.loop_nest().iterate([&](auto const& idxs){
                        IdxType l = idxs[0], s = idxs[1], r = idxs[2];
                        Cplx    v; T.get(idxs, gsl::span<Cplx>(&v,1));
                        if(s == 0) tmp0[r] += env[l] * v;
                        else       tmp1[r] += env[l] * v;
                    });
        
                    // 5b) Compute probabilities
                    double w0 = 0.0, w1 = 0.0;
                    for(IdxType r = 0; r < Dr; ++r) {
                        w0 += std::norm(tmp0[r]);
                        w1 += std::norm(tmp1[r]);
                    }
                    double total = w0 + w1;
                    double p0    = w0 / total;
        
                    // 5c) Sample outcome
                    bool outcome1 = (dist(rng) >= p0);
        
                    // 5d) Pack bit into its physical position (first sample → bit 0, etc.)
                    if(outcome1) packed |= (IdxType(1) << (j - 1));
        
                    // 5e) Update env
                    double invnorm = 1.0 / std::sqrt(outcome1 ? w1 : w0);
                    env.clear();
                    env.reserve(Dr);
                    auto& chosen = outcome1 ? tmp1 : tmp0;
                    for(IdxType r = 0; r < Dr; ++r) {
                        env.push_back(chosen[r] * invnorm);
                    }
                }
        
                // 6) Store result
                results[rep] = packed;
        
                // 7) Cleanup work
                for(auto& Tw : work) {
                    Tw.deallocate();
                }
                work.clear();
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

        void apply_two_qubit(const std::array<Cplx,16>& U4, IdxType q0, IdxType q1) {
            // 0) Initial diagnostics
            printf("[DEBUG] Qubit sites: %zu, %zu\n", (size_t)q0, (size_t)q1);
            std::array<Cplx,16> U4_eff = U4;  // copy to mutable buffer
        
            // Dump U4_eff
            {
                IdxType d = phys_dims[q0];
                printf("[DEBUG] U4_eff (row-major %zux%zu):\n", (size_t)d, (size_t)d);
                for(size_t i = 0; i < d*d; ++i) {
                    Cplx v = U4_eff[i];
                    printf("  U4_eff[%2zu] = (%.6f, %.6f)\n",
                           i, (double)std::real(v), (double)std::imag(v));
                }
            }
        
            // Bond dims
            printf("[DEBUG] bond_dims =");
            for(auto bd : bond_dims) printf(" %zu", (size_t)bd);
            printf("\n");
        
            // 1) Local dimensions
            IdxType Dl = bond_dims[q0];
            IdxType Dr = bond_dims[q1+1];
            IdxType d  = phys_dims[q0];
            printf("[DEBUG] Dl=%zu, d=%zu, Dr=%zu\n", (size_t)Dl, (size_t)d, (size_t)Dr);
        
            // 2) Merge MPS cores into M(l,p0,p1,r)
            tamm::Tensor<Cplx> M({ bond_tis[q0], phys_tis[q0], phys_tis[q1], bond_tis[q1+1] });
            M.set_dense(); M.allocate(&ec);
            {
                tamm::Scheduler sch{ec};
                sch(M("l","p0","p1","r") =
                      mps_tensors[q0]("l","p0","b") *
                      mps_tensors[q1]("b","p1","r"),
                    "merge_two", tamm::ExecutionHW::GPU);
                sch.execute(tamm::ExecutionHW::GPU);
            }
            printf("[DEBUG] Merged cores q%zu and q%zu into M\n", (size_t)q0, (size_t)q1);
        
            // 3) Build two-qubit gate tensor G4(p0',p1',p0,p1)
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
                int row = p0 * d + p1;
                int col = p0p * d + p1p;
                G4.put(idxs, gsl::span<Cplx>(&U4_eff[row*(d*d) + col],1));
            });
            printf("[DEBUG] Constructed G4\n");
        
            // 4) Apply G4 → M2
            tamm::Tensor<Cplx> M2({ bond_tis[q0], phys_tis[q0], phys_tis[q1], bond_tis[q1+1] });
            M2.set_dense(); M2.allocate(&ec);
            {
                tamm::Scheduler sch{ec};
                sch(M2("l","p0p","p1p","r") =
                      G4("p0p","p1p","p0","p1") * M("l","p0","p1","r"),
                    "apply_two", tamm::ExecutionHW::GPU);
                sch.execute(tamm::ExecutionHW::GPU);
            }
            M.deallocate(); G4.deallocate();
            printf("[DEBUG] Applied G4 to M → M2\n");
        
            // 5) Flatten M2 into Eigen matrix (Dl·d)×(d·Dr)
            Eigen::Index rows = Dl * d;
            Eigen::Index cols = d  * Dr;
            printf("[DEBUG] Preparing Eigen matrix (%lld × %lld)\n",
                   (long long)rows, (long long)cols);
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> mat(rows, cols);
            M2.loop_nest().iterate([&](auto const& idxs){
                IdxType l  = idxs[0];
                IdxType p0 = idxs[1];
                IdxType p1 = idxs[2];
                IdxType r  = idxs[3];
                Cplx v; M2.get(idxs, gsl::span<Cplx>(&v,1));
                mat(l*d + p0, p1*Dr + r) = v;
            });
            printf("[DEBUG] Populated Eigen matrix\n");
        
            // 6) SVD + truncate
            Eigen::JacobiSVD<decltype(mat)> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
            auto svals = svd.singularValues();
            IdxType chi = std::min<IdxType>(max_bond_dim, (IdxType)svals.size());
            printf("[DEBUG] SVD singular values =");
            for(int i = 0; i < std::min<int>(5, (int)svals.size()); ++i)
                printf(" %.6f", (double)svals(i));
            double ssum = svals.sum();
            printf(" ... sum=%.6f, chi=%zu\n", ssum, (size_t)chi);
        
            auto Umat  = svd.matrixU().leftCols(chi);
            auto Sdiag = svals.head(chi).asDiagonal();
            auto Vh    = svd.matrixV().leftCols(chi).adjoint();
        
            // 7) Update bond metadata
            bond_dims[q0+1] = chi;
            {
                tamm::IndexSpace is_new{tamm::range(chi)};
                bond_tis[q0+1]   = tamm::TiledIndexSpace(is_new,1);
            }
            printf("[DEBUG] Updated bond_dims[%zu] = %zu\n", (size_t)(q0+1), (size_t)chi);
        
            // 8) Build left-canonical Ti_new
            tamm::Tensor<Cplx> Ti_new({
                bond_tis[q0], phys_tis[q0],
                tamm::TiledIndexSpace(tamm::range(chi),1)
            });
            Ti_new.set_dense(); Ti_new.allocate(&ec);
            Ti_new.loop_nest().iterate([&](auto const& idxs){
                IdxType l = idxs[0], p0 = idxs[1], b = idxs[2];
                Cplx val = Umat(l*d + p0, b);
                Ti_new.put(idxs, gsl::span<Cplx>(&val,1));
            });
        
            // 9) Build center Tj_new = Sdiag·Vh
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
            printf("[DEBUG] Constructed Ti_new and Tj_new\n");
        
            // 10) Replace old tensors
            mps_tensors[q0].deallocate();
            mps_tensors[q1].deallocate();
            mps_tensors[q0] = std::move(Ti_new);
            mps_tensors[q1] = std::move(Tj_new);
            M2.deallocate();
            printf("[DEBUG] Replaced mps_tensors at sites %zu and %zu\n", (size_t)q0, (size_t)q1);
        
            // 11) Post-gate norm checks (optional)
            {
                tamm::Tensor<Cplx> R({ bond_tis[q0], phys_tis[q0], bond_tis[q0+1] });
                R.set_dense(); R.allocate(&ec);
                tamm::Scheduler sch{ec};
                sch(R("l","p","r") =
                      mps_tensors[q0]("l","p","r") *
                      tamm::conj(mps_tensors[q0])("l","p","r"),
                    "norm_q0", tamm::ExecutionHW::GPU);
                sch.execute(tamm::ExecutionHW::GPU);
                double norm0 = 0;
                R.loop_nest().iterate([&](auto const& idxs){
                    Cplx v; R.get(idxs, gsl::span<Cplx>(&v,1));
                    norm0 += std::real(v);
                });
                R.deallocate();
                printf("[DEBUG] post-gate norm(q%zu) = %f\n", (size_t)q0, norm0);
            }
            {
                tamm::Tensor<Cplx> R({ bond_tis[q1], phys_tis[q1], bond_tis[q1+1] });
                R.set_dense(); R.allocate(&ec);
                tamm::Scheduler sch{ec};
                sch(R("l","p","r") =
                      mps_tensors[q1]("l","p","r") *
                      tamm::conj(mps_tensors[q1])("l","p","r"),
                    "norm_q1", tamm::ExecutionHW::GPU);
                sch.execute(tamm::ExecutionHW::GPU);
                double norm1 = 0;
                R.loop_nest().iterate([&](auto const& idxs){
                    Cplx v; R.get(idxs, gsl::span<Cplx>(&v,1));
                    norm1 += std::real(v);
                });
                R.deallocate();
                printf("[DEBUG] post-gate norm(q%zu) = %f\n", (size_t)q1, norm1);
            }
        
            printf("[DEBUG] apply_two_qubit completed for sites %zu, %zu\n\n", (size_t)q0, (size_t)q1);
        }

        void left_canonicalize(std::vector<tamm::Tensor<Cplx>>& MPS) {
          using Eigen::Index;
        
          for(IdxType i = 0; i < n_qubits - 1; ++i) {
            //–– 1) Read dimensions
            IdxType Dl     = bond_dims[i];
            IdxType Dr_old = bond_dims[i+1];
            IdxType d      = phys_dims[i];
            printf("[LC] sweeping bond %zu (sites %zu↔%zu): Dl=%zu, d=%zu, Dr_old=%zu\n",
                   (size_t)i, (size_t)i, (size_t)(i+1),
                   (size_t)Dl, (size_t)d, (size_t)Dr_old);
        
            //–– 2) Build flattened Mmat ∈ ℂ^{(Dl·d)×Dr_old}
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Mmat(Dl*d, Dr_old);
            MPS[i].loop_nest().iterate([&](auto const& idxs){
              IdxType l = idxs[0], p = idxs[1], r = idxs[2];
              Cplx    v; MPS[i].get(idxs, gsl::span<Cplx>(&v,1));
              Mmat(l*d + p, r) = v;
            });
            printf("[LC]   built Mmat of size %zux%zu\n",
                   (size_t)(Dl*d), (size_t)Dr_old);
        
            //–– Debug: print every entry of Mmat
            printf("[LC]   dumping Mmat entries:\n");
            for(Index row = 0; row < Mmat.rows(); ++row) {
              for(Index col = 0; col < Mmat.cols(); ++col) {
                Cplx v = Mmat(row,col);
                printf("    Mmat(%lld,%lld) = (%.6f, %.6f)\n",
                       (long long)row, (long long)col,
                       (double)std::real(v), (double)std::imag(v));
              }
            }
        
            //–– 3) Compute thin SVD: Mmat = Umat · Sdiag · Vh
            Eigen::JacobiSVD<decltype(Mmat)> svd(
              Mmat, Eigen::ComputeThinU | Eigen::ComputeThinV
            );
            auto svals = svd.singularValues();
            Index full_rank = (Index)svals.size();
            Index chi = std::min<Index>((Index)max_bond_dim, full_rank);
            printf("[LC]   SVD full_rank=%lld, truncating to χ=%zu\n",
                   (long long)full_rank, (size_t)chi);
        
            auto Umat  = svd.matrixU().leftCols(chi);           // (Dl·d)×χ
            auto Sdiag = svals.head(chi).asDiagonal();          // χ×χ
            auto Vh    = svd.matrixV().leftCols(chi).adjoint(); // χ×Dr_old
            printf("[LC]   extracted Umat(%zux%zu), Sdiag(%zux%zu), Vh(%zux%zu)\n",
                   (size_t)(Dl*d), (size_t)chi,
                   (size_t)chi,     (size_t)chi,
                   (size_t)chi,     (size_t)Dr_old);
        
            //–– 4) Update bond dimension & create new left tensor
            bond_dims[i+1] = chi;
            {
              tamm::IndexSpace is_new{tamm::range(chi)};
              bond_tis[i+1]   = tamm::TiledIndexSpace(is_new, 1);
            }
            tamm::Tensor<Cplx> Tleft({bond_tis[i], phys_tis[i], bond_tis[i+1]});
            Tleft.set_dense(); Tleft.allocate(&ec);
            Tleft.loop_nest().iterate([&](auto const& idxs){
              IdxType l = idxs[0], p = idxs[1], b = idxs[2];
              Cplx    v = Umat(l*d + p, b);
              Tleft.put(idxs, gsl::span<Cplx>(&v,1));
            });
            printf("[LC]   built Tleft(site=%zu) with new bond dim=%zu\n",
                   (size_t)i, (size_t)chi);
        
            //–– Debug: dump Tleft entries
            printf("[LC]   dumping Tleft entries:\n");
            Tleft.loop_nest().iterate([&](auto const& idxs){
              Cplx v; Tleft.get(idxs, gsl::span<Cplx>(&v,1));
              printf("    Tleft(%zu,%zu,%zu) = (%.6f, %.6f)\n",
                     (size_t)idxs[0], (size_t)idxs[1], (size_t)idxs[2],
                     (double)std::real(v), (double)std::imag(v));
            });
        
            //–– 5) Absorb S⋅Vh into the next tensor
            auto M2 = Sdiag * Vh;  // χ×Dr_old
            printf("[LC]   formed M2 matrix of size %zux%zu\n",
                   (size_t)chi, (size_t)Dr_old);
        
            auto& Tnext_old = MPS[i+1];
            IdxType d_next   = phys_dims[i+1];
            IdxType Dr_right = bond_dims[i+2];
        
            // 5a) Flatten old Tnext to Dr_old × (d_next·Dr_right)
            Eigen::Matrix<Cplx, Eigen::Dynamic, Eigen::Dynamic> Tnext_mat(
              Dr_old, d_next * Dr_right
            );
            Tnext_old.loop_nest().iterate([&](auto const& idxs){
              IdxType b = idxs[0], p = idxs[1], r = idxs[2];
              Cplx    v; Tnext_old.get(idxs, gsl::span<Cplx>(&v,1));
              Tnext_mat(b, p*Dr_right + r) = v;
            });
            printf("[LC]   flattened Tnext_old(site=%zu) to %zux%zu\n",
                   (size_t)(i+1),
                   (size_t)Dr_old,
                   (size_t)(d_next * Dr_right)
            );
        
            //–– Debug: dump Tnext_mat entries
            printf("[LC]   dumping Tnext_mat entries:\n");
            for(Index br = 0; br < Tnext_mat.rows(); ++br) {
              for(Index cr = 0; cr < Tnext_mat.cols(); ++cr) {
                Cplx v = Tnext_mat(br,cr);
                printf("    Tnext_mat(%lld,%lld) = (%.6f, %.6f)\n",
                       (long long)br, (long long)cr,
                       (double)std::real(v), (double)std::imag(v));
              }
            }
        
            // 5b) Multiply: new_flat = M2 (χ×Dr_old) × Tnext_mat (Dr_old×(d_next·Dr_right))
            auto Tnext_mat2 = M2 * Tnext_mat;  // χ×(d_next·Dr_right)
            printf("[LC]   applied M2 → new matrix %zux%zu\n",
                   (size_t)chi,
                   (size_t)(d_next * Dr_right)
            );
        
            //–– 5c) Reshape back into Tnext(site=i+1)
            tamm::Tensor<Cplx> Tnext({bond_tis[i+1], phys_tis[i+1], bond_tis[i+2]});
            Tnext.set_dense(); Tnext.allocate(&ec);
            Tnext.loop_nest().iterate([&](auto const& idxs){
              IdxType b = idxs[0], p = idxs[1], r = idxs[2];
              Cplx    v = Tnext_mat2(b, p*Dr_right + r);
              Tnext.put(idxs, gsl::span<Cplx>(&v,1));
            });
            printf("[LC]   rebuilt Tnext(site=%zu)\n", (size_t)(i+1));
        
            //–– Debug: dump Tnext entries
            printf("[LC]   dumping Tnext entries:\n");
            Tnext.loop_nest().iterate([&](auto const& idxs){
              Cplx v; Tnext.get(idxs, gsl::span<Cplx>(&v,1));
              printf("    Tnext(%zu,%zu,%zu) = (%.6f, %.6f)\n",
                     (size_t)idxs[0], (size_t)idxs[1], (size_t)idxs[2],
                     (double)std::real(v), (double)std::imag(v));
            });
        
            //–– 6) Swap in the new cores
            MPS[i]     .deallocate();
            Tnext_old  .deallocate();
            MPS[i]     = std::move(Tleft);
            MPS[i+1]   = std::move(Tnext);
            printf("[LC]   swapped in new cores for sites %zu and %zu\n\n",
                   (size_t)i, (size_t)(i+1));
          }
        
          printf("[LC] left_canonicalize complete; sites 0…%zu are left-iso, centre at site %zu\n",
                 (size_t)(n_qubits-2), (size_t)(n_qubits-1));
        }


        static tamm::ProcGroup init_pg() {
            int argc = 0; char** argv = nullptr;
            //MPI_Init(&argc,&argv);
            //GA_Initialize();
            //tamm::initialize(argc, argv);
            //tamm::ProcGroup::self_ga_pgroup(true);
            return tamm::ProcGroup::create_world_coll();
        }
    };
} // namespace NWQS
