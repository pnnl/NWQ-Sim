#include "backend/statevector_gpu.hpp"

#ifdef VQE_ENABLE_HIP

#include "backend/gate_builder.hpp"
#include "backend/gate_fusion.hpp"
#include "backend/statevector_cpu.hpp"

#include <hip/hip_cooperative_groups.h>
#include <hip/hip_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

constexpr int THREADS_CTA_HIP = 256;

namespace NWQSim
{
  using IdxType = long long int;
  constexpr double S2I = 0.70710678118654752440;
}

namespace vqe::backend
{

  namespace cg = cooperative_groups;

  namespace
  {

    inline void check_hip(hipError_t err, const char *what)
    {
      if (err != hipSuccess)
      {
        throw std::runtime_error(std::string(what) + ": " + hipGetErrorString(err));
      }
    }

    struct device_state
    {
      NWQSim::IdxType dimension = 0;
      NWQSim::IdxType half_dimension = 0;
      std::uint32_t num_qubits = 0;
      double *real = nullptr;
      double *imag = nullptr;
    };

    struct expectation_parameters
    {
      std::size_t total_work = 0;
      std::uint64_t x_mask = 0;
      std::uint64_t z_mask = 0;
      double sign = 0.0;
      int odd_y = 0;
      std::uint64_t mask_lower = 0;
      std::uint64_t mask_upper = 0;
    };

    __device__ __constant__ double kGateHReal[4] = {NWQSim::S2I, NWQSim::S2I, NWQSim::S2I, -NWQSim::S2I};
    __device__ __constant__ double kGateHImag[4] = {0.0, 0.0, 0.0, 0.0};

    __device__ __constant__ double kGateYForwardReal[4] = {0.5, 0.5, 0.5, -0.5};
    __device__ __constant__ double kGateYForwardImag[4] = {-0.5, 0.5, -0.5, -0.5};
    __device__ __constant__ double kGateYBackwardReal[4] = {0.5, 0.5, 0.5, -0.5};
    __device__ __constant__ double kGateYBackwardImag[4] = {0.5, 0.5, -0.5, 0.5};

    __device__ __constant__ double kGateCNOTReal[16] = {
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 1.0, 0.0};
    __device__ __constant__ double kGateCNOTImag[16] = {0.0};

    __device__ inline void apply_single_device(const device_state &st,
                                               cg::grid_group &grid,
                                               std::uint32_t qubit,
                                               const double *gm_real,
                                               const double *gm_imag)
    {
      const NWQSim::IdxType q = static_cast<NWQSim::IdxType>(qubit);
      const NWQSim::IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
      const NWQSim::IdxType half_dim = st.half_dimension;
      const NWQSim::IdxType dim = st.dimension;
      if (half_dim == 0)
      {
        if (tid == 0 && q == 0 && dim == 1)
        {
          double elr = st.real[0];
          double eli = st.imag[0];
          double in_real[2] = {elr, 0.0};
          double in_imag[2] = {eli, 0.0};
          double out_real[2] = {0.0, 0.0};
          double out_imag[2] = {0.0, 0.0};
          for (int row = 0; row < 2; ++row)
          {
            double sum_r = 0.0;
            double sum_i = 0.0;
            for (int col = 0; col < 2; ++col)
            {
              const int idx = row * 2 + col;
              const double mr = gm_real[idx];
              const double mi = gm_imag[idx];
              const double vr = in_real[col];
              const double vi = in_imag[col];
              sum_r += mr * vr - mi * vi;
              sum_i += mr * vi + mi * vr;
            }
            out_real[row] = sum_r;
            out_imag[row] = sum_i;
          }
          st.real[0] = out_real[0];
          st.imag[0] = out_imag[0];
        }
        grid.sync();
        return;
      }
      const NWQSim::IdxType stride_limit = half_dim;
      for (NWQSim::IdxType i = tid; i < stride_limit; i += blockDim.x * gridDim.x)
      {
        NWQSim::IdxType outer = (i >> q);
        NWQSim::IdxType inner = (i & ((((NWQSim::IdxType)1) << q) - 1));
        NWQSim::IdxType offset = (outer << (q + 1));
        NWQSim::IdxType pos0 = offset + inner;
        NWQSim::IdxType pos1 = offset + inner + (((NWQSim::IdxType)1) << q);

        const double el0_real = st.real[pos0];
        const double el0_imag = st.imag[pos0];
        const double el1_real = st.real[pos1];
        const double el1_imag = st.imag[pos1];

        const double m00r = gm_real[0];
        const double m01r = gm_real[1];
        const double m10r = gm_real[2];
        const double m11r = gm_real[3];
        const double m00i = gm_imag[0];
        const double m01i = gm_imag[1];
        const double m10i = gm_imag[2];
        const double m11i = gm_imag[3];

        const double out0_r = m00r * el0_real - m00i * el0_imag + m01r * el1_real - m01i * el1_imag;
        const double out0_i = m00r * el0_imag + m00i * el0_real + m01r * el1_imag + m01i * el1_real;
        const double out1_r = m10r * el0_real - m10i * el0_imag + m11r * el1_real - m11i * el1_imag;
        const double out1_i = m10r * el0_imag + m10i * el0_real + m11r * el1_imag + m11i * el1_real;

        st.real[pos0] = out0_r;
        st.imag[pos0] = out0_i;
        st.real[pos1] = out1_r;
        st.imag[pos1] = out1_i;
      }
      grid.sync();
    }

    __device__ inline void apply_two_device(const device_state &st,
                                            cg::grid_group &grid,
                                            std::uint32_t control,
                                            std::uint32_t target,
                                            const double *gm_real,
                                            const double *gm_imag)
    {
      const NWQSim::IdxType q0 = static_cast<NWQSim::IdxType>(control);
      const NWQSim::IdxType q1 = static_cast<NWQSim::IdxType>(target);
      const NWQSim::IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
      const NWQSim::IdxType per_work = st.dimension >> 2;
      if (per_work == 0)
      {
        grid.sync();
        return;
      }
      const NWQSim::IdxType max_q = q0 > q1 ? q0 : q1;
      const NWQSim::IdxType min_q = q0 < q1 ? q0 : q1;
      const NWQSim::IdxType q0dim = ((NWQSim::IdxType)1) << max_q;
      const NWQSim::IdxType q1dim = ((NWQSim::IdxType)1) << min_q;
      const NWQSim::IdxType mider_factor = (q0dim + q1dim + q1dim - 1) >> (min_q + 1);
      const NWQSim::IdxType inner_factor = q1dim;
      const NWQSim::IdxType qubit0_dim = ((NWQSim::IdxType)1) << q0;
      const NWQSim::IdxType qubit1_dim = ((NWQSim::IdxType)1) << q1;

      for (NWQSim::IdxType i = tid; i < per_work; i += blockDim.x * gridDim.x)
      {
        NWQSim::IdxType outer = ((i / inner_factor) / mider_factor) * (q0dim + q0dim);
        NWQSim::IdxType mider = ((i / inner_factor) % mider_factor) * (q1dim + q1dim);
        NWQSim::IdxType inner = i % inner_factor;
        NWQSim::IdxType pos0 = outer + mider + inner;
        NWQSim::IdxType pos1 = pos0 + qubit1_dim;
        NWQSim::IdxType pos2 = pos0 + qubit0_dim;
        NWQSim::IdxType pos3 = pos0 + q0dim + q1dim;

        const double el0_real = st.real[pos0];
        const double el0_imag = st.imag[pos0];
        const double el1_real = st.real[pos1];
        const double el1_imag = st.imag[pos1];
        const double el2_real = st.real[pos2];
        const double el2_imag = st.imag[pos2];
        const double el3_real = st.real[pos3];
        const double el3_imag = st.imag[pos3];

        double vec_real[4] = {el0_real, el1_real, el2_real, el3_real};
        double vec_imag[4] = {el0_imag, el1_imag, el2_imag, el3_imag};
        double out_real[4] = {0.0, 0.0, 0.0, 0.0};
        double out_imag[4] = {0.0, 0.0, 0.0, 0.0};

#pragma unroll
        for (int row = 0; row < 4; ++row)
        {
          double sum_r = 0.0;
          double sum_i = 0.0;
#pragma unroll
          for (int col = 0; col < 4; ++col)
          {
            const int idx = row * 4 + col;
            const double mr = gm_real[idx];
            const double mi = gm_imag[idx];
            const double vr = vec_real[col];
            const double vi = vec_imag[col];
            sum_r += mr * vr - mi * vi;
            sum_i += mr * vi + mi * vr;
          }
          out_real[row] = sum_r;
          out_imag[row] = sum_i;
        }

        st.real[pos0] = out_real[0];
        st.imag[pos0] = out_imag[0];
        st.real[pos1] = out_real[1];
        st.imag[pos1] = out_imag[1];
        st.real[pos2] = out_real[2];
        st.imag[pos2] = out_imag[2];
        st.real[pos3] = out_real[3];
        st.imag[pos3] = out_imag[3];
      }
      grid.sync();
    }

    __device__ inline void apply_global_phase(device_state &st,
                                              cg::grid_group &grid,
                                              double angle)
    {
      const double half = angle * 0.5;
      const double c = cos(half);
      const double s = sin(half);
      const NWQSim::IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
      const NWQSim::IdxType dim = st.dimension;
      for (NWQSim::IdxType idx = tid; idx < dim; idx += blockDim.x * gridDim.x)
      {
        double r = st.real[idx];
        double i = st.imag[idx];
        st.real[idx] = c * r + s * i;
        st.imag[idx] = -s * r + c * i;
      }
      grid.sync();
    }

    __device__ inline void apply_rz_device(device_state &st,
                                           cg::grid_group &grid,
                                           std::uint32_t qubit,
                                           double angle)
    {
      const double half = angle * 0.5;
      const double c = cos(half);
      const double s = sin(half);
      double gm_real[4] = {c, 0.0, 0.0, c};
      double gm_imag[4] = {-s, 0.0, 0.0, s};
      apply_single_device(st, grid, qubit, gm_real, gm_imag);
    }

    __device__ inline void apply_pauli_device(device_state &st,
                                              cg::grid_group &grid,
                                              const sim_gate &gate)
    {
      const std::uint32_t count = gate.pauli_qubit_count;
      if (count == 0)
      {
        apply_global_phase(st, grid, gate.angle);
        return;
      }

      for (std::uint32_t idx = 0; idx < count; ++idx)
      {
        const std::uint32_t qubit = gate.pauli_qubits[idx];
        const std::uint8_t basis = gate.pauli_basis[idx];
        if (basis == 1)
        {
          apply_single_device(st, grid, qubit, kGateHReal, kGateHImag);
        }
        else if (basis == 2)
        {
          apply_single_device(st, grid, qubit, kGateYForwardReal, kGateYForwardImag);
        }
      }

      for (std::uint32_t idx = count; idx > 1; --idx)
      {
        const std::uint32_t ctrl = gate.pauli_qubits[idx - 1];
        const std::uint32_t tgt = gate.pauli_qubits[idx - 2];
        apply_two_device(st, grid, ctrl, tgt, kGateCNOTReal, kGateCNOTImag);
      }

      apply_rz_device(st, grid, gate.pauli_qubits[0], gate.angle);

      for (std::uint32_t idx = 1; idx < count; ++idx)
      {
        const std::uint32_t ctrl = gate.pauli_qubits[idx];
        const std::uint32_t tgt = gate.pauli_qubits[idx - 1];
        apply_two_device(st, grid, ctrl, tgt, kGateCNOTReal, kGateCNOTImag);
      }

      for (std::uint32_t idx = 0; idx < count; ++idx)
      {
        const std::uint32_t qubit = gate.pauli_qubits[idx];
        const std::uint8_t basis = gate.pauli_basis[idx];
        if (basis == 1)
        {
          apply_single_device(st, grid, qubit, kGateHReal, kGateHImag);
        }
        else if (basis == 2)
        {
          apply_single_device(st, grid, qubit, kGateYBackwardReal, kGateYBackwardImag);
        }
      }
    }

    __global__ void apply_gates_kernel(device_state state,
                                       const sim_gate *gates,
                                       std::size_t gate_count)
    {
      cg::grid_group grid = cg::this_grid();
      for (std::size_t idx = 0; idx < gate_count; ++idx)
      {
        const sim_gate &gate = gates[idx];
        if (gate.op == sim_gate::kind::single)
        {
          apply_single_device(state, grid, gate.target, gate.real.data(), gate.imag.data());
        }
        else if (gate.op == sim_gate::kind::two)
        {
          apply_two_device(state, grid, gate.control, gate.target, gate.real.data(), gate.imag.data());
        }
        else if (gate.op == sim_gate::kind::pauli)
        {
          apply_pauli_device(state, grid, gate);
        }
        grid.sync();
      }
    }

    __global__ void expectation_terms_kernel(device_state state,
                                             const expectation_parameters *params,
                                             const double2 *coefficients,
                                             double2 *results,
                                             std::size_t term_count)
    {
      const std::size_t term_idx = static_cast<std::size_t>(blockIdx.x);
      if (term_idx >= term_count)
      {
        return;
      }

      const expectation_parameters param = params[term_idx];
      const double2 coeff = coefficients[term_idx];
      const std::size_t total_work = param.total_work;
      const std::uint64_t x_mask = param.x_mask;
      const std::uint64_t z_mask = param.z_mask;
      const double sign = param.sign;
      const bool odd_y = param.odd_y != 0;
      const std::uint64_t mask_lower = param.mask_lower;
      const std::uint64_t mask_upper = param.mask_upper;

      double accum = 0.0;
      const double *real = state.real;
      const double *imag = state.imag;

      if (x_mask == 0ULL)
      {
        for (std::size_t idx = threadIdx.x; idx < total_work; idx += blockDim.x)
        {
          const double r = real[idx];
          const double im = imag[idx];
          double contrib = r * r + im * im;
          if ((__popcll(static_cast<unsigned long long>(idx & z_mask)) & 1) != 0)
          {
            contrib = -contrib;
          }
          accum += contrib;
        }
      }
      else
      {
        for (std::size_t i = threadIdx.x; i < total_work; i += blockDim.x)
        {
          const std::uint64_t base = static_cast<std::uint64_t>(i);
          const std::uint64_t idx0 = ((base << 1ULL) & mask_upper) | (base & mask_lower);
          const std::uint64_t idx1 = idx0 ^ x_mask;

          const double a0_r = real[static_cast<std::size_t>(idx0)];
          const double a0_i = imag[static_cast<std::size_t>(idx0)];
          const double a1_r = real[static_cast<std::size_t>(idx1)];
          const double a1_i = imag[static_cast<std::size_t>(idx1)];

          double v0;
          double v1;

          if (odd_y)
          {
            v0 = -a1_i * a0_r + a0_i * a1_r;
            v1 = -a0_i * a1_r + a1_i * a0_r;
          }
          else
          {
            v0 = a1_r * a0_r + a1_i * a0_i;
            v1 = a1_r * a0_r + a0_i * a1_i;
          }

          v0 *= sign;
          v1 *= sign;

          double contrib = v0;
          if ((__popcll(static_cast<unsigned long long>(idx0 & z_mask)) & 1) != 0)
          {
            contrib = -contrib;
          }

          if ((__popcll(static_cast<unsigned long long>(idx1 & z_mask)) & 1) != 0)
          {
            contrib -= v1;
          }
          else
          {
            contrib += v1;
          }

          accum += contrib;
        }
      }

      extern __shared__ double shared[];
      shared[threadIdx.x] = accum;
      __syncthreads();

      for (unsigned int offset = blockDim.x >> 1; offset > 0; offset >>= 1)
      {
        if (threadIdx.x < offset)
        {
          shared[threadIdx.x] += shared[threadIdx.x + offset];
        }
        __syncthreads();
      }

      if (threadIdx.x == 0)
      {
        const double expectation_value = shared[0];
        results[term_idx] = make_double2(coeff.x * expectation_value, coeff.y * expectation_value);
      }
    }

  } // namespace

  struct statevector_gpu::impl
  {
    explicit impl(std::size_t qubits)
    {
      num_qubits = qubits;
      dimension = qubits == 0 ? std::size_t{1} : (std::size_t{1} << qubits);
      half_dimension = qubits < 1 ? std::size_t{0} : (dimension >> 1);

      check_hip(hipSetDevice(0), "hipSetDevice");
      check_hip(hipGetDeviceProperties(&device_properties, 0), "hipGetDeviceProperties init");
      check_hip(hipMalloc(&d_real, dimension * sizeof(double)), "hipMalloc real");
      check_hip(hipMalloc(&d_imag, dimension * sizeof(double)), "hipMalloc imag");
      check_hip(hipMemset(d_real, 0, dimension * sizeof(double)), "hipMemset real");
      check_hip(hipMemset(d_imag, 0, dimension * sizeof(double)), "hipMemset imag");
      double one = 1.0;
      check_hip(hipMemcpy(d_real, &one, sizeof(double), hipMemcpyHostToDevice), "hipMemcpy real init");
      check_hip(hipStreamCreateWithFlags(&compute_stream, hipStreamNonBlocking), "hipStreamCreate compute");
      check_hip(hipStreamCreateWithFlags(&expectation_stream, hipStreamNonBlocking), "hipStreamCreate expectation");

      compute_block = dim3(THREADS_CTA_HIP, 1, 1);
      int blocks_per_sm = 0;
      check_hip(hipOccupancyMaxActiveBlocksPerMultiprocessor(&blocks_per_sm, apply_gates_kernel, THREADS_CTA_HIP, 0), "hipOccupancyMaxActiveBlocksPerMultiprocessor init");
      const int max_blocks = std::max(1, blocks_per_sm * device_properties.multiProcessorCount);
      const std::size_t work_items = std::max<std::size_t>(1, dimension);
      const int required_blocks = static_cast<int>((work_items + THREADS_CTA_HIP - 1) / THREADS_CTA_HIP);
      const int total_blocks = std::max(1, std::min(max_blocks, required_blocks));
      compute_grid = dim3(static_cast<unsigned int>(total_blocks), 1, 1);
    }

    ~impl()
    {
      if (h_expectation_params_pinned != nullptr)
      {
        hipHostFree(h_expectation_params_pinned);
      }
      if (h_expectation_coeffs_pinned != nullptr)
      {
        hipHostFree(h_expectation_coeffs_pinned);
      }
      if (h_expectation_results_pinned != nullptr)
      {
        hipHostFree(h_expectation_results_pinned);
      }
      if (d_expectation_params != nullptr)
      {
        hipFree(d_expectation_params);
      }
      if (d_expectation_coeffs != nullptr)
      {
        hipFree(d_expectation_coeffs);
      }
      if (d_expectation_results != nullptr)
      {
        hipFree(d_expectation_results);
      }
      if (d_gates != nullptr)
      {
        hipFree(d_gates);
      }
      if (d_real != nullptr)
      {
        hipFree(d_real);
      }
      if (d_imag != nullptr)
      {
        hipFree(d_imag);
      }
      if (compute_stream != nullptr)
      {
        hipStreamDestroy(compute_stream);
      }
      if (expectation_stream != nullptr)
      {
        hipStreamDestroy(expectation_stream);
      }
    }

    void reset()
    {
      check_hip(hipMemset(d_real, 0, dimension * sizeof(double)), "hipMemset real");
      check_hip(hipMemset(d_imag, 0, dimension * sizeof(double)), "hipMemset imag");
      double one = 1.0;
      check_hip(hipMemcpy(d_real, &one, sizeof(double), hipMemcpyHostToDevice), "hipMemcpy reset");
      amplitudes_valid = false;
      cached_expectation_valid = false;
    }

    void ensure_gate_capacity(std::size_t count)
    {
      if (count <= gate_capacity)
      {
        return;
      }
      if (d_gates != nullptr)
      {
        hipFree(d_gates);
      }
      gate_capacity = std::max<std::size_t>(count, 1);
      check_hip(hipMalloc(&d_gates, gate_capacity * sizeof(sim_gate)), "hipMalloc gates");
    }

    void apply(const circuit &circ)
    {
      build_simulation_gates(circ, gate_buffer);
      fuse_simulation_gates(gate_buffer, num_qubits, fuse_buffer, fuse_tmp1, fuse_tmp2, fuse_tmp3, fuse_chunk, fused_buffer);

      const std::size_t gate_count = fused_buffer.size();
      if (gate_count == 0)
      {
        amplitudes_valid = false;
        return;
      }

      ensure_gate_capacity(gate_count);
      check_hip(hipMemcpy(d_gates, fused_buffer.data(), gate_count * sizeof(sim_gate), hipMemcpyHostToDevice), "hipMemcpy gates");

      device_state state{};
      state.dimension = static_cast<NWQSim::IdxType>(dimension);
      state.half_dimension = static_cast<NWQSim::IdxType>(half_dimension);
      state.num_qubits = static_cast<std::uint32_t>(num_qubits);
      state.real = d_real;
      state.imag = d_imag;

      void *args[] = {&state, &d_gates, const_cast<std::size_t *>(&gate_count)};
      check_hip(hipLaunchCooperativeKernel(reinterpret_cast<void *>(apply_gates_kernel), compute_grid, compute_block, args, 0, compute_stream), "hipLaunchCooperativeKernel apply gates");
      check_hip(hipStreamSynchronize(compute_stream), "hipStreamSynchronize apply");

      amplitudes_valid = false;
    }

    const std::vector<std::complex<double>> &amplitudes() const
    {
      if (!amplitudes_valid)
      {
        std::vector<double> host_real(dimension);
        std::vector<double> host_imag(dimension);
        check_hip(hipMemcpy(host_real.data(), d_real, dimension * sizeof(double), hipMemcpyDeviceToHost), "hipMemcpy amplitudes real");
        check_hip(hipMemcpy(host_imag.data(), d_imag, dimension * sizeof(double), hipMemcpyDeviceToHost), "hipMemcpy amplitudes imag");
        amplitude_cache.resize(dimension);
        for (std::size_t idx = 0; idx < dimension; ++idx)
        {
          amplitude_cache[idx] = {host_real[idx], host_imag[idx]};
        }
        amplitudes_valid = true;
      }
      return amplitude_cache;
    }

    void load_state(const std::vector<std::complex<double>> &amplitudes_vec)
    {
      if (amplitudes_vec.size() != dimension)
      {
        throw std::invalid_argument("amplitude vector size mismatch");
      }
      std::vector<double> host_real(dimension);
      std::vector<double> host_imag(dimension);
      for (std::size_t idx = 0; idx < dimension; ++idx)
      {
        host_real[idx] = amplitudes_vec[idx].real();
        host_imag[idx] = amplitudes_vec[idx].imag();
      }
      check_hip(hipMemcpy(d_real, host_real.data(), dimension * sizeof(double), hipMemcpyHostToDevice), "hipMemcpy load real");
      check_hip(hipMemcpy(d_imag, host_imag.data(), dimension * sizeof(double), hipMemcpyHostToDevice), "hipMemcpy load imag");
      amplitudes_valid = false;
      cached_expectation_valid = false;
    }

    void copy_state_from(const impl &other)
    {
      if (other.dimension != dimension)
      {
        throw std::invalid_argument("state dimension mismatch in copy_state_from");
      }
      check_hip(hipMemcpy(d_real, other.d_real, dimension * sizeof(double), hipMemcpyDeviceToDevice), "hipMemcpy copy real");
      check_hip(hipMemcpy(d_imag, other.d_imag, dimension * sizeof(double), hipMemcpyDeviceToDevice), "hipMemcpy copy imag");
      amplitudes_valid = false;
      cached_expectation_valid = false;
    }

    std::complex<double> expectation(const pauli_term &term) const
    {
      expectation_host_complex.resize(1);
      compute_expectation_terms_gpu(&term, 1, expectation_host_complex.data());
      return expectation_host_complex[0];
    }

    std::complex<double> expectation(const std::vector<pauli_term> &terms) const
    {
      if (terms.empty())
      {
        return {0.0, 0.0};
      }

      expectation_host_complex.resize(terms.size());
      compute_expectation_terms_gpu(terms.data(), terms.size(), expectation_host_complex.data());

      std::complex<double> total{0.0, 0.0};
      for (const auto &value : expectation_host_complex)
      {
        total += value;
      }
      return total;
    }

    void ensure_expectation_term_capacity(std::size_t count) const
    {
      if (count <= expectation_term_capacity)
      {
        return;
      }
      if (count > expectation_term_capacity)
      {
        if (d_expectation_params != nullptr)
        {
          hipFree(d_expectation_params);
        }
        if (d_expectation_coeffs != nullptr)
        {
          hipFree(d_expectation_coeffs);
        }
        if (d_expectation_results != nullptr)
        {
          hipFree(d_expectation_results);
        }
        expectation_term_capacity = std::max<std::size_t>(count, 1);
        check_hip(hipMalloc(&d_expectation_params, expectation_term_capacity * sizeof(expectation_parameters)), "hipMalloc expectation params");
        check_hip(hipMalloc(&d_expectation_coeffs, expectation_term_capacity * sizeof(double2)), "hipMalloc expectation coeffs");
        check_hip(hipMalloc(&d_expectation_results, expectation_term_capacity * sizeof(double2)), "hipMalloc expectation results");
        ensure_expectation_host_capacity(expectation_term_capacity);
        cached_expectation_valid = false;
      }
    }

    void ensure_expectation_host_capacity(std::size_t count) const
    {
      if (count <= expectation_host_pinned_capacity)
      {
        return;
      }
      if (count > expectation_host_pinned_capacity)
      {
        if (h_expectation_params_pinned != nullptr)
        {
          hipHostFree(h_expectation_params_pinned);
        }
        if (h_expectation_coeffs_pinned != nullptr)
        {
          hipHostFree(h_expectation_coeffs_pinned);
        }
        if (h_expectation_results_pinned != nullptr)
        {
          hipHostFree(h_expectation_results_pinned);
        }
        expectation_host_pinned_capacity = std::max<std::size_t>(count, 1);
        check_hip(hipHostMalloc(reinterpret_cast<void **>(&h_expectation_params_pinned), expectation_host_pinned_capacity * sizeof(expectation_parameters)), "hipHostMalloc expectation params");
        check_hip(hipHostMalloc(reinterpret_cast<void **>(&h_expectation_coeffs_pinned), expectation_host_pinned_capacity * sizeof(double2)), "hipHostMalloc expectation coeffs");
        check_hip(hipHostMalloc(reinterpret_cast<void **>(&h_expectation_results_pinned), expectation_host_pinned_capacity * sizeof(double2)), "hipHostMalloc expectation results");
        cached_expectation_valid = false;
      }
    }

    void compute_expectation_terms_gpu(const pauli_term *terms,
                                       std::size_t term_count,
                                       std::complex<double> *out) const
    {
      if (term_count == 0)
      {
        return;
      }

      ensure_expectation_term_capacity(term_count);
      ensure_expectation_host_capacity(term_count);

      expectation_parameters *host_params = h_expectation_params_pinned;
      double2 *host_coeffs = h_expectation_coeffs_pinned;
      if (!cached_expectation_valid || cached_expectation_terms != terms || cached_expectation_count != term_count)
      {
        for (std::size_t idx = 0; idx < term_count; ++idx)
        {
          const auto &term = terms[idx];
          expectation_parameters param{};
          param.x_mask = term.x_mask;
          param.z_mask = term.z_mask;
          if (term.x_mask == 0ULL)
          {
            param.sign = 1.0;
            param.odd_y = 0;
            param.mask_lower = 0ULL;
            param.mask_upper = 0ULL;
            param.total_work = dimension;
          }
          else
          {
            const std::uint64_t y_mask = term.x_mask & term.z_mask;
            const unsigned int y_count = static_cast<unsigned int>(__builtin_popcountll(static_cast<unsigned long long>(y_mask)));
            param.sign = ((y_count / 2U) & 1U) ? -1.0 : 1.0;
            param.odd_y = static_cast<int>(y_count & 1U);
            const unsigned int leading = static_cast<unsigned int>(__builtin_clzll(term.x_mask));
            const unsigned int max_x = 63U - leading;
            param.mask_lower = (max_x == 0U) ? 0ULL : ((std::uint64_t{1} << max_x) - 1ULL);
            if (max_x + 1U >= 64U)
            {
              param.mask_upper = 0ULL;
            }
            else
            {
              param.mask_upper = ~((std::uint64_t{1} << (max_x + 1U)) - 1ULL);
            }
            param.total_work = half_dimension;
          }

          host_params[idx] = param;
          host_coeffs[idx] = make_double2(term.coefficient.real(), term.coefficient.imag());
        }

        check_hip(hipMemcpyAsync(d_expectation_params, host_params, term_count * sizeof(expectation_parameters), hipMemcpyHostToDevice, expectation_stream), "hipMemcpyAsync expectation params");
        check_hip(hipMemcpyAsync(d_expectation_coeffs, host_coeffs, term_count * sizeof(double2), hipMemcpyHostToDevice, expectation_stream), "hipMemcpyAsync expectation coeffs");

        cached_expectation_terms = terms;
        cached_expectation_count = term_count;
        cached_expectation_valid = true;
      }

      device_state state{};
      state.dimension = static_cast<NWQSim::IdxType>(dimension);
      state.half_dimension = static_cast<NWQSim::IdxType>(half_dimension);
      state.num_qubits = static_cast<std::uint32_t>(num_qubits);
      state.real = d_real;
      state.imag = d_imag;

      const int threads = THREADS_CTA_HIP;
      const std::size_t shared = static_cast<std::size_t>(threads) * sizeof(double);
      hipLaunchKernelGGL(expectation_terms_kernel, dim3(static_cast<unsigned int>(term_count)), dim3(threads), shared, expectation_stream, state, d_expectation_params, d_expectation_coeffs, d_expectation_results, term_count);
      check_hip(hipGetLastError(), "expectation_terms_kernel launch");

      check_hip(hipMemcpyAsync(h_expectation_results_pinned, d_expectation_results, term_count * sizeof(double2), hipMemcpyDeviceToHost, expectation_stream), "hipMemcpyAsync expectation results");
      check_hip(hipStreamSynchronize(expectation_stream), "hipStreamSynchronize expectation");

      const double2 *result_data = h_expectation_results_pinned;
      for (std::size_t idx = 0; idx < term_count; ++idx)
      {
        out[idx] = std::complex<double>(result_data[idx].x, result_data[idx].y);
      }
    }

    std::size_t num_qubits = 0;
    std::size_t dimension = 0;
    std::size_t half_dimension = 0;
    double *d_real = nullptr;
    double *d_imag = nullptr;
    sim_gate *d_gates = nullptr;
    std::size_t gate_capacity = 0;

    mutable std::vector<std::complex<double>> amplitude_cache;
    mutable bool amplitudes_valid = false;

    mutable expectation_parameters *d_expectation_params = nullptr;
    mutable double2 *d_expectation_coeffs = nullptr;
    mutable double2 *d_expectation_results = nullptr;
    mutable std::size_t expectation_term_capacity = 0;
    mutable expectation_parameters *h_expectation_params_pinned = nullptr;
    mutable double2 *h_expectation_coeffs_pinned = nullptr;
    mutable double2 *h_expectation_results_pinned = nullptr;
    mutable std::size_t expectation_host_pinned_capacity = 0;
    hipStream_t compute_stream = nullptr;
    hipStream_t expectation_stream = nullptr;
    mutable std::vector<std::complex<double>> expectation_host_complex;
    hipDeviceProp_t device_properties{};
    dim3 compute_grid{};
    dim3 compute_block{};

    std::vector<sim_gate> gate_buffer;
    std::vector<sim_gate> fuse_tmp1;
    std::vector<sim_gate> fuse_tmp2;
    std::vector<sim_gate> fuse_tmp3;
    std::vector<sim_gate> fuse_buffer;
    std::vector<sim_gate> fuse_chunk;
    std::vector<sim_gate> fused_buffer;

    mutable const pauli_term *cached_expectation_terms = nullptr;
    mutable std::size_t cached_expectation_count = 0;
    mutable bool cached_expectation_valid = false;
  };

  statevector_gpu::statevector_gpu(std::size_t num_qubits)
      : impl_(std::make_unique<impl>(num_qubits)) {}

  statevector_gpu::~statevector_gpu() = default;

  void statevector_gpu::reset()
  {
    impl_->reset();
  }

  void statevector_gpu::apply(const circuit &circ)
  {
    impl_->apply(circ);
  }

  std::complex<double> statevector_gpu::expectation(const pauli_term &term) const
  {
    return impl_->expectation(term);
  }

  std::complex<double> statevector_gpu::expectation(const std::vector<pauli_term> &terms) const
  {
    return impl_->expectation(terms);
  }

  const std::vector<std::complex<double>> &statevector_gpu::amplitudes() const
  {
    return impl_->amplitudes();
  }

  void statevector_gpu::load_state(const std::vector<std::complex<double>> &amplitudes_vec)
  {
    impl_->load_state(amplitudes_vec);
  }

  void statevector_gpu::copy_state_from(const statevector_gpu &other)
  {
    impl_->copy_state_from(*other.impl_);
  }

  std::size_t statevector_gpu::num_qubits() const
  {
    return impl_->num_qubits;
  }

} // namespace vqe::backend

#endif // VQE_ENABLE_HIP
