#pragma once

#include <cstring>

#if MEMORY_HAS_CUDA

#include <random>
#include <vector>

namespace NWQSim
{
    namespace Memory
    {
        namespace GPU
        {
            inline void free_measurement_buffers(MemoryHandle &handle)
            {
                auto &gpu = handle.extras.gpu;
                if (gpu.results_gpu)
                {
                    if (gpu.results_device_is_nvshmem)
                        nvshmem_free(gpu.results_gpu);
                    else
                        cudaFree(gpu.results_gpu);
                }
                gpu.results_gpu = nullptr;
                gpu.results_device_is_nvshmem = false;

                if (gpu.randoms_gpu)
                    cudaFree(gpu.randoms_gpu);
                gpu.randoms_gpu = nullptr;

                if (gpu.results_host)
                    cudaFreeHost(gpu.results_host);
                gpu.results_host = nullptr;

                if (gpu.randoms_host)
                    cudaFreeHost(gpu.randoms_host);
                gpu.randoms_host = nullptr;

                gpu.measurement_slots = 0;
            }
        } // namespace GPU

        inline void initialize_gpu(MemoryHandle &handle, SimType method, IdxType n_qubits)
        {
            GateKernels::StateView &view = handle.view;

            check_cuda(cudaSetDevice(0), "cudaSetDevice");

            IdxType log_dim = (method == SimType::SV) ? n_qubits : (2 * n_qubits);
            IdxType dim = (log_dim > 0) ? (IdxType(1) << log_dim) : 1;
            IdxType half_dim = (log_dim > 0) ? (dim >> 1) : 0;

            const size_t state_bytes = static_cast<size_t>(dim) * sizeof(ValType);

            ValType *host_real = nullptr;
            ValType *host_imag = nullptr;
            check_cuda(cudaMallocHost(reinterpret_cast<void **>(&host_real), state_bytes), "cudaMallocHost real");
            register_allocation(handle, host_real, MemoryHandle::AllocationKind::CudaHost, state_bytes);
            check_cuda(cudaMallocHost(reinterpret_cast<void **>(&host_imag), state_bytes), "cudaMallocHost imag");
            register_allocation(handle, host_imag, MemoryHandle::AllocationKind::CudaHost, state_bytes);
            std::memset(host_real, 0, state_bytes);
            std::memset(host_imag, 0, state_bytes);
            if (dim > 0)
                host_real[0] = 1.0;

            ValType *dev_real = nullptr;
            ValType *dev_imag = nullptr;
            check_cuda(cudaMalloc(reinterpret_cast<void **>(&dev_real), state_bytes), "cudaMalloc real");
            register_allocation(handle, dev_real, MemoryHandle::AllocationKind::CudaDevice, state_bytes);
            check_cuda(cudaMalloc(reinterpret_cast<void **>(&dev_imag), state_bytes), "cudaMalloc imag");
            register_allocation(handle, dev_imag, MemoryHandle::AllocationKind::CudaDevice, state_bytes);

            size_t aux_bytes = state_bytes + sizeof(ValType);
            ValType *buffer_real = nullptr;
            ValType *buffer_imag = nullptr;
            check_cuda(cudaMalloc(reinterpret_cast<void **>(&buffer_real), aux_bytes), "cudaMalloc buffer real");
            register_allocation(handle, buffer_real, MemoryHandle::AllocationKind::CudaDevice, aux_bytes);
            check_cuda(cudaMalloc(reinterpret_cast<void **>(&buffer_imag), aux_bytes), "cudaMalloc buffer imag");
            register_allocation(handle, buffer_imag, MemoryHandle::AllocationKind::CudaDevice, aux_bytes);

            check_cuda(cudaMemcpy(dev_real, host_real, state_bytes, cudaMemcpyHostToDevice), "cudaMemcpy real");
            check_cuda(cudaMemcpy(dev_imag, host_imag, state_bytes, cudaMemcpyHostToDevice), "cudaMemcpy imag");
            check_cuda(cudaMemset(buffer_real, 0, aux_bytes), "cudaMemset buffer real");
            check_cuda(cudaMemset(buffer_imag, 0, aux_bytes), "cudaMemset buffer imag");

            view.data_real = dev_real;
            view.data_imag = dev_imag;
            view.buffer_real = buffer_real;
            view.buffer_imag = buffer_imag;
            view.dim = dim;
            view.half_dim = half_dim;
            view.gpu_scale = 0;
            view.lg2_m_gpu = log_dim;
            view.m_gpu = dim;
            view.rank = 0;

            handle.extras.host_real = host_real;
            handle.extras.host_imag = host_imag;
            handle.extras.host_bytes = state_bytes;
        }

        inline void reset_gpu(MemoryHandle &handle)
        {
            GateKernels::StateView &view = handle.view;
            size_t state_bytes = handle.extras.host_bytes;
            if (state_bytes == 0)
                return;

            std::memset(handle.extras.host_real, 0, state_bytes);
            std::memset(handle.extras.host_imag, 0, state_bytes);
            if (view.dim > 0)
                handle.extras.host_real[0] = 1.0;

            size_t aux_bytes = state_bytes + sizeof(ValType);
            check_cuda(cudaMemcpy(view.data_real, handle.extras.host_real, state_bytes, cudaMemcpyHostToDevice), "cudaMemcpy real");
            check_cuda(cudaMemcpy(view.data_imag, handle.extras.host_imag, state_bytes, cudaMemcpyHostToDevice), "cudaMemcpy imag");
            check_cuda(cudaMemset(view.buffer_real, 0, aux_bytes), "cudaMemset buffer real");
            check_cuda(cudaMemset(view.buffer_imag, 0, aux_bytes), "cudaMemset buffer imag");
        }

        inline void load_gpu(MemoryHandle &handle, const ValType *real, const ValType *imag)
        {
            GateKernels::StateView &view = handle.view;
            size_t state_bytes = handle.extras.host_bytes;
            std::memcpy(handle.extras.host_real, real, state_bytes);
            std::memcpy(handle.extras.host_imag, imag, state_bytes);
            check_cuda(cudaMemcpy(view.data_real, handle.extras.host_real, state_bytes, cudaMemcpyHostToDevice), "cudaMemcpy real");
            check_cuda(cudaMemcpy(view.data_imag, handle.extras.host_imag, state_bytes, cudaMemcpyHostToDevice), "cudaMemcpy imag");
        }

        inline void dump_gpu(const MemoryHandle &handle, ValType *real_out, ValType *imag_out)
        {
            const GateKernels::StateView &view = handle.view;
            size_t state_bytes = handle.extras.host_bytes;
            check_cuda(cudaMemcpy(handle.extras.host_real, view.data_real, state_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy real D2H");
            check_cuda(cudaMemcpy(handle.extras.host_imag, view.data_imag, state_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy imag D2H");
            std::memcpy(real_out, handle.extras.host_real, state_bytes);
            std::memcpy(imag_out, handle.extras.host_imag, state_bytes);
        }

        namespace GPU
        {
            inline DeviceGate *upload_gates(MemoryHandle &handle, const std::vector<DeviceGate> &gates)
            {
                auto &gpu = handle.extras.gpu;
                size_t bytes = gates.size() * sizeof(DeviceGate);

                if (bytes == 0)
                {
                    if (gpu.gates_gpu)
                        cudaFree(gpu.gates_gpu);
                    gpu.gates_gpu = nullptr;
                    gpu.gates_bytes = 0;
                    gpu.gates_count = 0;
                    return nullptr;
                }

                if (!gpu.gates_gpu || gpu.gates_bytes != bytes)
                {
                    if (gpu.gates_gpu)
                        cudaFree(gpu.gates_gpu);
                    check_cuda(cudaMalloc(reinterpret_cast<void **>(&gpu.gates_gpu), bytes), "cudaMalloc gates");
                    gpu.gates_bytes = bytes;
                }

                check_cuda(cudaMemcpy(gpu.gates_gpu, gates.data(), bytes, cudaMemcpyHostToDevice), "cudaMemcpy gates");
                gpu.gates_count = static_cast<IdxType>(gates.size());
                return gpu.gates_gpu;
            }

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

                if (slots == 0)
                {
                    free_measurement_buffers(handle);
                    return 0;
                }

                free_measurement_buffers(handle);

                auto &gpu = handle.extras.gpu;

                check_cuda(cudaMallocHost(reinterpret_cast<void **>(&gpu.results_host), slots * sizeof(IdxType)), "cudaMallocHost results");
                std::memset(gpu.results_host, 0, slots * sizeof(IdxType));

                check_cuda(cudaMalloc(reinterpret_cast<void **>(&gpu.results_gpu), slots * sizeof(IdxType)), "cudaMalloc results");
                check_cuda(cudaMemset(gpu.results_gpu, 0, slots * sizeof(IdxType)), "cudaMemset results");

                check_cuda(cudaMallocHost(reinterpret_cast<void **>(&gpu.randoms_host), slots * sizeof(ValType)), "cudaMallocHost randoms");
                std::mt19937 rng(seed ? seed : std::random_device{}());
                std::uniform_real_distribution<ValType> dist(0.0, 1.0);
                for (IdxType i = 0; i < slots; ++i)
                    gpu.randoms_host[i] = dist(rng);

                check_cuda(cudaMalloc(reinterpret_cast<void **>(&gpu.randoms_gpu), slots * sizeof(ValType)), "cudaMalloc randoms");
                check_cuda(cudaMemcpy(gpu.randoms_gpu, gpu.randoms_host, slots * sizeof(ValType), cudaMemcpyHostToDevice), "cudaMemcpy randoms");

                gpu.measurement_slots = slots;
                gpu.results_device_is_nvshmem = false;
                return slots;
            }

            inline DeviceGate *device_gates(const MemoryHandle &handle)
            {
                return handle.extras.gpu.gates_gpu;
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
                auto &gpu = handle.extras.gpu;
                if (!gpu.results_gpu || !gpu.results_host || gpu.measurement_slots == 0)
                    return;
                size_t bytes = static_cast<size_t>(gpu.measurement_slots) * sizeof(IdxType);
                check_cuda(cudaMemcpy(gpu.results_host, gpu.results_gpu, bytes, cudaMemcpyDeviceToHost), "cudaMemcpy results D2H");
            }
        } // namespace GPU
    } // namespace Memory
} // namespace NWQSim

#endif // MEMORY_HAS_CUDA
