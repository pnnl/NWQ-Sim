#ifndef VQE_CUDA_STATE
#define VQE_CUDA_STATE
#include "svsim/sv_cuda_mpi.cuh"
#include "vqe_state.hpp"
#include "observable/pauli_operator.hpp"
#include "utils.hpp"
#include "circuit/ansatz.hpp"
#include "circuit/measurement.hpp"
#include "gradient/sa_gradient.hpp"
#include "observable/hamiltonian.hpp"
#include "nlopt.hpp"
#include <memory>
#include <cmath>

namespace NWQSim
{
  namespace VQE {
    class SV_CUDA_VQE: public VQEState, public SV_CUDA_MPI {
      public:
        SV_CUDA_VQE(std::shared_ptr<Ansatz> a, 
                   const Hamiltonian& h, 
                   nlopt::algorithm optimizer_algorithm,
                   Callback _callback,
                   IdxType seed = 0,
                   OptimizerSettings opt_settings = OptimizerSettings()): 
                                      SV_CUDA(a->num_qubits()),
                                      VQEState(a, h, optimizer_algorithm, _callback, seed, opt_settings) {};
      virtual void call_simulator(std::shared_ptr<Ansatz> ansatz) override {        
        reset_state();
        sim(ansatz);
      };
    __global__ void simulation_kernel_cuda_mpi(SV_CUDA_MPI *sv_gpu, IdxType n_gates) override
    {
        IdxType cur_index = 0;
        IdxType lg2_m_gpu = sv_gpu->lg2_m_gpu;
        grid_group grid = this_grid();
        bool already_sync = false;

        for (IdxType t = 0; t < n_gates; t++)
        {
            OP op_name = (sv_gpu->gates_gpu)[t].op_name;
            IdxType qubit = (sv_gpu->gates_gpu)[t].qubit;

            IdxType ctrl = (sv_gpu->gates_gpu)[t].ctrl;
            ValType *gm_real = (sv_gpu->gates_gpu)[t].gm_real;
            ValType *gm_imag = (sv_gpu->gates_gpu)[t].gm_imag;

            IdxType repetition = (sv_gpu->gates_gpu)[t].qubit;

            // only need sync when operating on remote qubits
            if ((ctrl >= lg2_m_gpu) || (qubit >= lg2_m_gpu))
            {
                if (!already_sync) // do not need repeated sync
                {
                    if (threadIdx.x == 0 && blockIdx.x == 0)
                        nvshmem_barrier_all();
                    grid.sync();
                }
            }
            already_sync = false;

            // only need sync when operating on remote qubits
            if (op_name == OP::C1)
            {
                sv_gpu->C1V2_GATE(gm_real, gm_imag, qubit);
                // sv_gpu->C1V1_GATE(gm_real, gm_imag, qubit);
            }
            else if (op_name == OP::C2)
            {
                if ((ctrl >= lg2_m_gpu) && (qubit >= lg2_m_gpu))
                {
                    sv_gpu->SWAP_GATE(0, ctrl);
                    BARR_NVSHMEM;
                    sv_gpu->C2V1_GATE(gm_real, gm_imag, 0, qubit);
                    BARR_NVSHMEM;
                    sv_gpu->SWAP_GATE(0, ctrl);
                }
                else
                {
                    sv_gpu->C2V1_GATE(gm_real, gm_imag, ctrl, qubit);
                }
            }
            else if (op_name == OP::RESET)
            {
                sv_gpu->RESET_GATE(qubit);
            }
            else if (op_name == OP::M)
            {
                sv_gpu->M_GATE(qubit, cur_index);
                cur_index++;
            }
            else if (op_name == OP::MA)
            {
                sv_gpu->MA_GATE(repetition, cur_index);
                cur_index += repetition;
            }

            // only need sync when operating on remote qubits
            if ((ctrl >= lg2_m_gpu) || (qubit >= lg2_m_gpu))
            {
                if (threadIdx.x == 0 && blockIdx.x == 0)
                    nvshmem_barrier_all();
                already_sync = true;
            }
            grid.sync();
        }
    }
      virtual ValType getPauliExpectation(const PauliOperator& op) override {
          IdxType qubit = 0;
          IdxType xmask = 0;
          IdxType zmask = 0;
          IdxType y_phase = 0;
          IdxType max_x = 0;
          for (auto pauli_op: *op.getOps()) {
            switch (pauli_op)
            {
              case PauliOp::X:
                xmask = xmask | (1 << qubit);
                max_x = qubit;
                break;
              case PauliOp::Y:
                xmask = xmask | (1 << qubit);
                zmask = zmask | (1 << qubit);
                max_x = qubit;
                y_phase += 1;
                break;
              case PauliOp::Z:
                zmask = zmask | (1ll << qubit);
                break;
              default:
                break;
            }
            qubit++;
          }
          ValType expectation = 0.0;
          if (xmask == 0) {
            for (IdxType i = 0; i < dim; i++) {
              ValType local_exp = sv_real[i] * sv_real[i] - sv_imag[i] * sv_imag[i];
              if (count_ones(zmask & i) & 1) {
                local_exp *= -1;
              }
              expectation += local_exp;
            }
            return expectation;
          }
          ValType sign = (y_phase / 2) % 2 ? -1: 1;
          ValType phase = y_phase % 2;
          size_t mask_u = ~((1lu << (max_x + 1)) - 1);
          size_t mask_l = (1lu << (max_x)) - 1;
          for (IdxType i = 0; i < half_dim; i++) {
            IdxType idx0 = ((i << 1) & mask_u) | (i & mask_l);
            IdxType idx1 = xmask ^ idx0;
            ValType v0, v1;
            if (phase) {
              v0 = -sv_imag[idx1] * sv_real[idx0] + sv_imag[idx0] * sv_real[idx1];
              v1 = -sv_imag[idx0] * sv_real[idx1] + sv_imag[idx1] * sv_real[idx0];
            } else {
              v0 = sv_real[idx1] * sv_real[idx0] + sv_imag[idx1] * sv_imag[idx0];
              v1 = sv_real[idx1] * sv_real[idx0] + sv_imag[idx0] * sv_imag[idx1];
            }
            v0 *= sign;
            v1 *= sign;
            ValType thisval = v0;
            if ((count_ones(idx0 & zmask) & 1) != 0) {
              thisval *= -1;
            }
            if ((count_ones(idx1 & zmask) & 1) != 0) {
              thisval -= v1;
            } else {
              thisval += v1;
            }
            expectation += thisval;
          }
          return expectation;
        }       

    };
  };
} // namespace NWQSim

#endif