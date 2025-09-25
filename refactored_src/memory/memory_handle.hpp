#pragma once

#include <array>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "state_view.hpp"
#include "../nwqsim_core.hpp"

#include <cmath>

#if __has_include(<cuda_runtime_api.h>)
#define MEMORY_HAS_CUDA 1
#include <cuda_runtime_api.h>
#else
#define MEMORY_HAS_CUDA 0
#endif

#if __has_include(<hip/hip_runtime_api.h>)
#define MEMORY_HAS_HIP 1
#include <hip/hip_runtime_api.h>
#else
#define MEMORY_HAS_HIP 0
#endif

#if __has_include(<nvshmem.h>) && __has_include(<nvshmemx.h>)
#define MEMORY_HAS_NVSHMEM 1
#include <nvshmem.h>
#include <nvshmemx.h>
#else
#define MEMORY_HAS_NVSHMEM 0
#endif

namespace NWQSim
{
    struct DeviceGate
    {
        enum class Type
        {
            Matrix,
            Measure,
            MeasureAll
        };

        Type type = Type::Matrix;
        IdxType qubit0 = 0;
        IdxType qubit1 = -1;
        IdxType qubit2 = -1;
        IdxType qubit3 = -1;
        IdxType repetition = 1;
        std::array<ValType, 16> gm_real{};
        std::array<ValType, 16> gm_imag{};
        IdxType matrix_dim = 2; // 2 -> single-qubit (4 entries), 4 -> two-qubit (16 entries)
    };
    namespace Memory
    {
        inline void check_cuda(
#if MEMORY_HAS_CUDA
            cudaError_t result, const char *what)
        {
            if (result != cudaSuccess)
            {
                throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(result));
            }
        }
#else
            int, const char *)
        {
        }
#endif

        inline void check_hip(
#if MEMORY_HAS_HIP
            hipError_t result, const char *what)
        {
            if (result != hipSuccess)
            {
                throw std::runtime_error(std::string(what) + ": " + hipGetErrorString(result));
            }
        }
#else
            int, const char *)
        {
        }
#endif

        inline void check_nvshmem(int result, const char *what)
        {
#if MEMORY_HAS_NVSHMEM
            if (result != 0)
            {
                throw std::runtime_error(std::string(what) + ": NVSHMEM error code " + std::to_string(result));
            }
#else
            (void)result;
            (void)what;
#endif
        }

        enum class Backend
        {
            CPU,
            OMP,
            GPU,
            MultiGPU
        };

        struct MemoryConfig
        {
            Backend backend = Backend::CPU;
            SimType method = SimType::SV;
            IdxType n_qubits = 0;
            std::vector<DeviceGate> gates;
            IdxType nodes = 1;
            uint64_t seed = 0;
        };

        class MemoryHandle
        {
        public:
            GateKernels::StateView view{};
            Backend backend = Backend::CPU;
            SimType method = SimType::SV;
            IdxType n_qubits = 0;
            IdxType nodes = 1;

            struct Extras
            {
                ValType *host_real = nullptr;
                ValType *host_imag = nullptr;
                size_t host_bytes = 0;
#if MEMORY_HAS_CUDA
            struct GPUExtras
            {
                DeviceGate *gates_gpu = nullptr;
                size_t gates_bytes = 0;
                IdxType gates_count = 0;
                IdxType *results_host = nullptr;
                IdxType *results_gpu = nullptr;
                ValType *randoms_host = nullptr;
                ValType *randoms_gpu = nullptr;
                IdxType measurement_slots = 0;
                bool results_device_is_nvshmem = false;
            } gpu;
#endif
        } extras;

            MemoryHandle() = default;

            MemoryHandle(MemoryHandle &&other) noexcept
            {
                *this = std::move(other);
            }

            MemoryHandle &operator=(MemoryHandle &&other) noexcept
            {
                if (this != &other)
                {
                    release();
                    view = other.view;
                    backend = other.backend;
                    method = other.method;
                    n_qubits = other.n_qubits;
                    nodes = other.nodes;
                    extras = other.extras;
                    allocations_ = std::move(other.allocations_);

                    other.view = {};
                    other.backend = Backend::CPU;
                    other.method = SimType::SV;
                    other.n_qubits = 0;
                    other.nodes = 1;
                    other.extras = {};
                }
                return *this;
            }

            ~MemoryHandle()
            {
                release();
            }

            MemoryHandle(const MemoryHandle &) = delete;
            MemoryHandle &operator=(const MemoryHandle &) = delete;

        private:
            enum class AllocationKind
            {
                Host,
#if MEMORY_HAS_CUDA
                CudaHost,
                CudaDevice,
#endif
#if MEMORY_HAS_HIP
                HipHost,
                HipDevice,
#endif
#if MEMORY_HAS_NVSHMEM
                Nvshmem,
#endif
            };

            struct Allocation
            {
                void *ptr = nullptr;
                AllocationKind kind = AllocationKind::Host;
                size_t bytes = 0;
            };

            std::vector<Allocation> allocations_;

            void release()
            {
                for (auto it = allocations_.rbegin(); it != allocations_.rend(); ++it)
                {
                    if (!it->ptr)
                        continue;
                    switch (it->kind)
                    {
                    case AllocationKind::Host:
                        std::free(it->ptr);
                        break;
#if MEMORY_HAS_CUDA
                    case AllocationKind::CudaHost:
                        cudaFreeHost(it->ptr);
                        break;
                    case AllocationKind::CudaDevice:
                        cudaFree(it->ptr);
                        break;
#endif
#if MEMORY_HAS_HIP
                    case AllocationKind::HipHost:
                        hipHostFree(it->ptr);
                        break;
                    case AllocationKind::HipDevice:
                        hipFree(it->ptr);
                        break;
#endif
#if MEMORY_HAS_NVSHMEM
                    case AllocationKind::Nvshmem:
                        nvshmem_free(it->ptr);
                        break;
#endif
                    }
                }
                allocations_.clear();
#if MEMORY_HAS_CUDA
                if (extras.gpu.gates_gpu)
                    cudaFree(extras.gpu.gates_gpu);
                extras.gpu.gates_gpu = nullptr;
                extras.gpu.gates_bytes = 0;

                if (extras.gpu.results_gpu)
                {
                    if (extras.gpu.results_device_is_nvshmem)
                        nvshmem_free(extras.gpu.results_gpu);
                    else
                        cudaFree(extras.gpu.results_gpu);
                }
                extras.gpu.results_gpu = nullptr;
                extras.gpu.results_device_is_nvshmem = false;

                if (extras.gpu.randoms_gpu)
                    cudaFree(extras.gpu.randoms_gpu);
                extras.gpu.randoms_gpu = nullptr;

                if (extras.gpu.results_host)
                    cudaFreeHost(extras.gpu.results_host);
                extras.gpu.results_host = nullptr;

                if (extras.gpu.randoms_host)
                    cudaFreeHost(extras.gpu.randoms_host);
                extras.gpu.randoms_host = nullptr;

                extras.gpu.measurement_slots = 0;
#endif
                view = {};
                extras = {};
            }

            void register_allocation(void *ptr, AllocationKind kind, size_t bytes)
            {
                if (!ptr)
                    throw std::bad_alloc();
                allocations_.push_back({ptr, kind, bytes});
            }

            friend void register_allocation(MemoryHandle &, void *, AllocationKind, size_t);
            friend MemoryHandle create_memory(const MemoryConfig &config);
            friend void reset_state(MemoryHandle &handle);
            friend void load_state(MemoryHandle &handle, const ValType *real, const ValType *imag);
            friend void dump_state(const MemoryHandle &handle, ValType *real_out, ValType *imag_out);

            friend void initialize_cpu(MemoryHandle &, SimType, IdxType);
            friend void reset_cpu(MemoryHandle &);
            friend void load_cpu(MemoryHandle &, const ValType *, const ValType *);
            friend void dump_cpu(const MemoryHandle &, ValType *, ValType *);

#if MEMORY_HAS_CUDA
            friend void initialize_gpu(MemoryHandle &, SimType, IdxType);
            friend void reset_gpu(MemoryHandle &);
            friend void load_gpu(MemoryHandle &, const ValType *, const ValType *);
            friend void dump_gpu(const MemoryHandle &, ValType *, ValType *);
#endif
#if MEMORY_HAS_NVSHMEM && MEMORY_HAS_CUDA
            friend void initialize_multi_gpu(MemoryHandle &, SimType, IdxType, IdxType);
            friend void reset_multi_gpu(MemoryHandle &);
            friend void load_multi_gpu(MemoryHandle &, const ValType *, const ValType *);
            friend void dump_multi_gpu(const MemoryHandle &, ValType *, ValType *);
#endif
        };

        inline void register_allocation(MemoryHandle &handle, void *ptr, MemoryHandle::AllocationKind kind, size_t bytes)
        {
            handle.register_allocation(ptr, kind, bytes);
        }

        MemoryHandle create_memory(const MemoryConfig &config);
        void reset_state(MemoryHandle &handle);
        void load_state(MemoryHandle &handle, const ValType *real, const ValType *imag);
        void dump_state(const MemoryHandle &handle, ValType *real_out, ValType *imag_out);

    } // namespace Memory
} // namespace NWQSim

#include "backends/cpu_memory.hpp"
#if MEMORY_HAS_CUDA
#include "backends/gpu_memory.hpp"
#endif
#if MEMORY_HAS_NVSHMEM && MEMORY_HAS_CUDA
#include "backends/multi_gpu_memory.hpp"
#endif

namespace NWQSim
{
    namespace Memory
    {
        MemoryHandle create_memory(const MemoryConfig &config)
        {
            MemoryHandle handle;
            handle.backend = config.backend;
            handle.method = config.method;
            handle.n_qubits = config.n_qubits;
            handle.nodes = config.nodes;

            switch (config.backend)
            {
            case Backend::CPU:
            case Backend::OMP:
                initialize_cpu(handle, config.method, config.n_qubits);
                break;
            case Backend::GPU:
#if MEMORY_HAS_CUDA
                initialize_gpu(handle, config.method, config.n_qubits);
                if (!config.gates.empty())
                    GPU::upload_gates(handle, config.gates);
                GPU::prepare_measurements(handle, config.gates, config.seed);
                break;
#else
                throw std::runtime_error("GPU backend requested but CUDA runtime is unavailable during build");
#endif
            case Backend::MultiGPU:
#if MEMORY_HAS_NVSHMEM && MEMORY_HAS_CUDA
                initialize_multi_gpu(handle, config.method, config.n_qubits, config.nodes);
                if (!config.gates.empty())
                    MultiGPU::upload_gates(handle, config.gates);
                MultiGPU::prepare_measurements(handle, config.gates, config.seed);
                break;
#else
                throw std::runtime_error("MultiGPU backend requested but NVSHMEM/CUDA support is unavailable during build");
#endif
            default:
                throw std::runtime_error("Unsupported backend for memory creation");
            }

            return handle;
        }

        void reset_state(MemoryHandle &handle)
        {
            switch (handle.backend)
            {
            case Backend::CPU:
            case Backend::OMP:
                reset_cpu(handle);
                break;
            case Backend::GPU:
#if MEMORY_HAS_CUDA
                reset_gpu(handle);
                break;
#else
                throw std::runtime_error("GPU reset requested without CUDA support");
#endif
            case Backend::MultiGPU:
#if MEMORY_HAS_NVSHMEM && MEMORY_HAS_CUDA
                reset_multi_gpu(handle);
                break;
#else
                throw std::runtime_error("MultiGPU reset requested without NVSHMEM/CUDA support");
#endif
            }
        }

        void load_state(MemoryHandle &handle, const ValType *real, const ValType *imag)
        {
            switch (handle.backend)
            {
            case Backend::CPU:
            case Backend::OMP:
                load_cpu(handle, real, imag);
                break;
            case Backend::GPU:
#if MEMORY_HAS_CUDA
                load_gpu(handle, real, imag);
                break;
#else
                throw std::runtime_error("GPU load requested without CUDA support");
#endif
            case Backend::MultiGPU:
#if MEMORY_HAS_NVSHMEM && MEMORY_HAS_CUDA
                load_multi_gpu(handle, real, imag);
                break;
#else
                throw std::runtime_error("MultiGPU load requested without NVSHMEM/CUDA support");
#endif
            }
        }

        void dump_state(const MemoryHandle &handle, ValType *real_out, ValType *imag_out)
        {
            switch (handle.backend)
            {
            case Backend::CPU:
            case Backend::OMP:
                dump_cpu(handle, real_out, imag_out);
                break;
            case Backend::GPU:
#if MEMORY_HAS_CUDA
                dump_gpu(handle, real_out, imag_out);
                break;
#else
                throw std::runtime_error("GPU dump requested without CUDA support");
#endif
            case Backend::MultiGPU:
#if MEMORY_HAS_NVSHMEM && MEMORY_HAS_CUDA
                dump_multi_gpu(handle, real_out, imag_out);
                break;
#else
                throw std::runtime_error("MultiGPU dump requested without NVSHMEM/CUDA support");
#endif
            }
        }

#if MEMORY_HAS_CUDA
        inline DeviceGate *upload_gates(MemoryHandle &handle, const std::vector<DeviceGate> &gates)
        {
            switch (handle.backend)
            {
            case Backend::GPU:
                return GPU::upload_gates(handle, gates);
#if MEMORY_HAS_NVSHMEM && MEMORY_HAS_CUDA
            case Backend::MultiGPU:
                return MultiGPU::upload_gates(handle, gates);
#endif
            default:
                throw std::runtime_error("Gate upload is only supported on GPU-based backends");
            }
        }

        inline IdxType prepare_measurements(MemoryHandle &handle,
                                             const std::vector<DeviceGate> &gates,
                                             uint64_t seed = 0)
        {
            switch (handle.backend)
            {
            case Backend::GPU:
                return GPU::prepare_measurements(handle, gates, seed);
#if MEMORY_HAS_NVSHMEM && MEMORY_HAS_CUDA
            case Backend::MultiGPU:
                return MultiGPU::prepare_measurements(handle, gates, seed);
#endif
            default:
                throw std::runtime_error("Measurement preparation is only available for GPU-based backends");
            }
        }

        inline void copy_measure_results_to_host(MemoryHandle &handle)
        {
            switch (handle.backend)
            {
            case Backend::GPU:
                GPU::copy_measure_results_to_host(handle);
                return;
#if MEMORY_HAS_NVSHMEM && MEMORY_HAS_CUDA
            case Backend::MultiGPU:
                MultiGPU::copy_measure_results_to_host(handle);
                return;
#endif
            default:
                throw std::runtime_error("Measurement result copy is only available for GPU-based backends");
            }
        }

        inline IdxType measurement_slots(const MemoryHandle &handle)
        {
            return handle.extras.gpu.measurement_slots;
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

        inline DeviceGate *device_gates(const MemoryHandle &handle)
        {
            return handle.extras.gpu.gates_gpu;
        }

        inline IdxType device_gate_count(const MemoryHandle &handle)
        {
            return handle.extras.gpu.gates_count;
        }
#endif // MEMORY_HAS_CUDA
    } // namespace Memory
} // namespace NWQSim
