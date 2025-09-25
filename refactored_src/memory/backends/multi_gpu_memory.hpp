#pragma once

#if MEMORY_HAS_NVSHMEM && MEMORY_HAS_CUDA

#include <cstring>
#include <cmath>
#include <random>
#include <vector>

#include "gpu_memory.hpp"

namespace NWQSim
{
    namespace Memory
    {
        namespace MultiGPU
        {
            inline DeviceGate *upload_gates(MemoryHandle &handle, const std::vector<DeviceGate> &gates)
            {
                return GPU::upload_gates(handle, gates);
            }
        }

        inline void initialize_multi_gpu(MemoryHandle &handle, SimType method, IdxType n_qubits, IdxType nodes)
        {
            static bool nvshmem_initialized = false;
            if (!nvshmem_initialized)
            {
                check_nvshmem(nvshmem_init(), "nvshmem_init");
                nvshmem_initialized = true;
            }

            int total_pes = nvshmem_n_pes();
            if (nodes > 0 && static_cast<int>(nodes) != total_pes)
            {
                throw std::runtime_error("Requested node count does not match NVSHMEM PE count");
            }

            int rank = nvshmem_my_pe();

            GateKernels::StateView &view = handle.view;

            IdxType log_dim = (method == SimType::SV) ? n_qubits : (2 * n_qubits);
            IdxType dim = (log_dim > 0) ? (IdxType(1) << log_dim) : 1;
            if (dim % total_pes != 0)
            {
                throw std::runtime_error("Global state dimension must be divisible by number of NVSHMEM PEs");
            }

            int tmp = total_pes;
            IdxType gpu_scale = 0;
            while (tmp > 1)
            {
                if (tmp % 2 != 0)
                    throw std::runtime_error("Number of GPUs must be a power of two for MultiGPU backend");
                tmp >>= 1;
                ++gpu_scale;
            }

            IdxType lg2_m_gpu = log_dim - gpu_scale;
            IdxType local_elems = dim >> gpu_scale;
            size_t slice_bytes = static_cast<size_t>(local_elems) * sizeof(ValType);
            IdxType half_dim = (log_dim > 0) ? (dim >> 1) : 0;

            check_cuda(cudaSetDevice(0), "cudaSetDevice (multi GPU)");

            ValType *host_real = nullptr;
            ValType *host_imag = nullptr;
            check_cuda(cudaMallocHost(reinterpret_cast<void **>(&host_real), slice_bytes), "cudaMallocHost multi real");
            register_allocation(handle, host_real, MemoryHandle::AllocationKind::CudaHost, slice_bytes);
            check_cuda(cudaMallocHost(reinterpret_cast<void **>(&host_imag), slice_bytes), "cudaMallocHost multi imag");
            register_allocation(handle, host_imag, MemoryHandle::AllocationKind::CudaHost, slice_bytes);
            std::memset(host_real, 0, slice_bytes);
            std::memset(host_imag, 0, slice_bytes);
            if (rank == 0 && local_elems > 0)
                host_real[0] = 1.0;

            ValType *real = static_cast<ValType *>(nvshmem_malloc(slice_bytes));
            register_allocation(handle, real, MemoryHandle::AllocationKind::Nvshmem, slice_bytes);
            ValType *imag = static_cast<ValType *>(nvshmem_malloc(slice_bytes));
            register_allocation(handle, imag, MemoryHandle::AllocationKind::Nvshmem, slice_bytes);
            ValType *buffer_real = static_cast<ValType *>(nvshmem_malloc(slice_bytes));
            register_allocation(handle, buffer_real, MemoryHandle::AllocationKind::Nvshmem, slice_bytes);
            ValType *buffer_imag = static_cast<ValType *>(nvshmem_malloc(slice_bytes));
            register_allocation(handle, buffer_imag, MemoryHandle::AllocationKind::Nvshmem, slice_bytes);

            if (!real || !imag || !buffer_real || !buffer_imag)
            {
                throw std::runtime_error("nvshmem_malloc failed for MultiGPU initialization");
            }

            check_cuda(cudaMemcpy(real, host_real, slice_bytes, cudaMemcpyHostToDevice), "cudaMemcpy multi real");
            check_cuda(cudaMemcpy(imag, host_imag, slice_bytes, cudaMemcpyHostToDevice), "cudaMemcpy multi imag");
            check_cuda(cudaMemset(buffer_real, 0, slice_bytes), "cudaMemset buffer real");
            check_cuda(cudaMemset(buffer_imag, 0, slice_bytes), "cudaMemset buffer imag");

            view.data_real = real;
            view.data_imag = imag;
            view.buffer_real = buffer_real;
            view.buffer_imag = buffer_imag;
            view.dim = dim;
            view.half_dim = half_dim;
            view.gpu_scale = gpu_scale;
            view.lg2_m_gpu = lg2_m_gpu;
            view.m_gpu = static_cast<IdxType>(IdxType(1) << lg2_m_gpu);
            view.rank = rank;

            handle.extras.host_real = host_real;
            handle.extras.host_imag = host_imag;
            handle.extras.host_bytes = slice_bytes;
            handle.nodes = total_pes;
        }

        inline void reset_multi_gpu(MemoryHandle &handle)
        {
            GateKernels::StateView &view = handle.view;
            size_t slice_bytes = handle.extras.host_bytes;

            std::memset(handle.extras.host_real, 0, slice_bytes);
            std::memset(handle.extras.host_imag, 0, slice_bytes);
            if (view.rank == 0 && slice_bytes > 0)
                handle.extras.host_real[0] = 1.0;

            check_cuda(cudaMemcpy(view.data_real, handle.extras.host_real, slice_bytes, cudaMemcpyHostToDevice), "cudaMemcpy reset real");
            check_cuda(cudaMemcpy(view.data_imag, handle.extras.host_imag, slice_bytes, cudaMemcpyHostToDevice), "cudaMemcpy reset imag");
            check_cuda(cudaMemset(view.buffer_real, 0, slice_bytes), "cudaMemset buffer real");
            check_cuda(cudaMemset(view.buffer_imag, 0, slice_bytes), "cudaMemset buffer imag");
        }

        inline void load_multi_gpu(MemoryHandle &handle, const ValType *real, const ValType *imag)
        {
            GateKernels::StateView &view = handle.view;
            size_t slice_bytes = handle.extras.host_bytes;
            std::memcpy(handle.extras.host_real, real, slice_bytes);
            std::memcpy(handle.extras.host_imag, imag, slice_bytes);
            check_cuda(cudaMemcpy(view.data_real, handle.extras.host_real, slice_bytes, cudaMemcpyHostToDevice), "cudaMemcpy load real");
            check_cuda(cudaMemcpy(view.data_imag, handle.extras.host_imag, slice_bytes, cudaMemcpyHostToDevice), "cudaMemcpy load imag");
        }

        inline void dump_multi_gpu(const MemoryHandle &handle, ValType *real_out, ValType *imag_out)
        {
            const GateKernels::StateView &view = handle.view;
            size_t slice_bytes = handle.extras.host_bytes;
            check_cuda(cudaMemcpy(handle.extras.host_real, view.data_real, slice_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy dump real");
            check_cuda(cudaMemcpy(handle.extras.host_imag, view.data_imag, slice_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy dump imag");
            std::memcpy(real_out, handle.extras.host_real, slice_bytes);
            std::memcpy(imag_out, handle.extras.host_imag, slice_bytes);
        }

        namespace MultiGPU
        {
            inline IdxType prepare_measurements(MemoryHandle &handle,
                                                const std::vector<DeviceGate> &gates,
                                                uint64_t seed = 0)
            {
                IdxType slots = 0;
                for (const auto &g : gates)
                {
                    if (g.type == DeviceGate::Type::Measure)
                        ++slots;
                    else if (g.type == DeviceGate::Type::MeasureAll)
                        slots += g.repetition;
                }

                GPU::free_measurement_buffers(handle);

                if (slots == 0)
                    return 0;

                auto &gpu = handle.extras.gpu;

                check_cuda(cudaMallocHost(reinterpret_cast<void **>(&gpu.results_host), slots * sizeof(IdxType)), "cudaMallocHost results");
                std::memset(gpu.results_host, 0, slots * sizeof(IdxType));

                gpu.results_gpu = static_cast<IdxType *>(nvshmem_malloc(slots * sizeof(IdxType)));
                if (!gpu.results_gpu)
                    throw std::runtime_error("nvshmem_malloc failed for measurement results");
                check_cuda(cudaMemset(gpu.results_gpu, 0, slots * sizeof(IdxType)), "cudaMemset results");

                check_cuda(cudaMallocHost(reinterpret_cast<void **>(&gpu.randoms_host), slots * sizeof(ValType)), "cudaMallocHost randoms");
                std::mt19937 rng(seed ? seed : std::random_device{}());
                std::uniform_real_distribution<ValType> dist(0.0, 1.0);
                for (IdxType i = 0; i < slots; ++i)
                    gpu.randoms_host[i] = dist(rng);

                check_cuda(cudaMalloc(reinterpret_cast<void **>(&gpu.randoms_gpu), slots * sizeof(ValType)), "cudaMalloc randoms");
                check_cuda(cudaMemcpy(gpu.randoms_gpu, gpu.randoms_host, slots * sizeof(ValType), cudaMemcpyHostToDevice), "cudaMemcpy randoms");

                gpu.measurement_slots = slots;
                gpu.results_device_is_nvshmem = true;

                return slots;
            }

            inline IdxType *device_measure_results(const MemoryHandle &handle)
            {
                return handle.extras.gpu.results_gpu;
            }

            inline ValType *device_measure_randoms(const MemoryHandle &handle)
            {
                return handle.extras.gpu.randoms_gpu;
            }

            inline IdxType *host_measure_results(MemoryHandle &handle)
            {
                return handle.extras.gpu.results_host;
            }

            inline ValType *host_measure_randoms(MemoryHandle &handle)
            {
                return handle.extras.gpu.randoms_host;
            }

            inline IdxType measurement_slots(const MemoryHandle &handle)
            {
                return handle.extras.gpu.measurement_slots;
            }

            inline void copy_measure_results_to_host(MemoryHandle &handle)
            {
                GPU::copy_measure_results_to_host(handle);
            }
        } // namespace MultiGPU
    } // namespace Memory
} // namespace NWQSim

#endif // MEMORY_HAS_NVSHMEM && MEMORY_HAS_CUDA
