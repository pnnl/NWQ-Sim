#pragma once

#include <cstring>

namespace NWQSim
{
    namespace Memory
    {
        inline void initialize_cpu(MemoryHandle &handle, SimType method, IdxType n_qubits)
        {
            GateKernels::StateView &view = handle.view;

            IdxType log_dim = (method == SimType::SV) ? n_qubits : (2 * n_qubits);
            IdxType dim = (log_dim > 0) ? (IdxType(1) << log_dim) : 1;
            IdxType half_dim = (log_dim > 0) ? (dim >> 1) : 0;

            const size_t state_bytes = static_cast<size_t>(dim) * sizeof(ValType);

            ValType *real = static_cast<ValType *>(std::malloc(state_bytes));
            register_allocation(handle, real, MemoryHandle::AllocationKind::Host, state_bytes);
            std::memset(real, 0, state_bytes);

            ValType *imag = static_cast<ValType *>(std::malloc(state_bytes));
            register_allocation(handle, imag, MemoryHandle::AllocationKind::Host, state_bytes);
            std::memset(imag, 0, state_bytes);

            ValType *aux = nullptr;
            if (method == SimType::SV)
            {
                size_t aux_bytes = state_bytes + sizeof(ValType);
                aux = static_cast<ValType *>(std::malloc(aux_bytes));
                register_allocation(handle, aux, MemoryHandle::AllocationKind::Host, aux_bytes);
                std::memset(aux, 0, aux_bytes);
            }
            else
            {
                IdxType diag_entries = (IdxType(1) << n_qubits) + 1;
                size_t diag_bytes = static_cast<size_t>(diag_entries) * sizeof(ValType);
                aux = static_cast<ValType *>(std::malloc(diag_bytes));
                register_allocation(handle, aux, MemoryHandle::AllocationKind::Host, diag_bytes);
                std::memset(aux, 0, diag_bytes);
            }

            real[0] = 1.0;

            view.data_real = real;
            view.data_imag = imag;
            view.buffer_real = aux;
            view.buffer_imag = nullptr;
            view.dim = dim;
            view.half_dim = half_dim;
            view.gpu_scale = 0;
            view.lg2_m_gpu = log_dim;
            view.m_gpu = dim;
            view.rank = 0;
        }

        inline void reset_cpu(MemoryHandle &handle)
        {
            GateKernels::StateView &view = handle.view;
            if (!view.data_real || !view.data_imag)
                return;

            std::memset(view.data_real, 0, static_cast<size_t>(view.dim) * sizeof(ValType));
            std::memset(view.data_imag, 0, static_cast<size_t>(view.dim) * sizeof(ValType));
            if (view.buffer_real)
            {
                size_t aux_bytes = (handle.method == SimType::SV)
                                       ? (static_cast<size_t>(view.dim) * sizeof(ValType) + sizeof(ValType))
                                       : (static_cast<size_t>((IdxType(1) << handle.n_qubits) + 1) * sizeof(ValType));
                std::memset(view.buffer_real, 0, aux_bytes);
            }

            if (view.dim > 0)
                view.data_real[0] = 1.0;
        }

        inline void load_cpu(MemoryHandle &handle, const ValType *real, const ValType *imag)
        {
            GateKernels::StateView &view = handle.view;
            size_t bytes = static_cast<size_t>(view.dim) * sizeof(ValType);
            std::memcpy(view.data_real, real, bytes);
            std::memcpy(view.data_imag, imag, bytes);
        }

        inline void dump_cpu(const MemoryHandle &handle, ValType *real_out, ValType *imag_out)
        {
            const GateKernels::StateView &view = handle.view;
            size_t bytes = static_cast<size_t>(view.dim) * sizeof(ValType);
            std::memcpy(real_out, view.data_real, bytes);
            std::memcpy(imag_out, view.data_imag, bytes);
        }
    } // namespace Memory
} // namespace NWQSim

