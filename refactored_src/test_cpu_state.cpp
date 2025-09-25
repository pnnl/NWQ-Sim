#include <cassert>
#include <cmath>
#include <iostream>

#include "gate_kernels/cpu_gates.hpp"
#include "memory/memory_handle.hpp"

int main()
{
    using namespace NWQSim;
    using namespace NWQSim::GateKernels;

    Memory::MemoryConfig config{};
    config.backend = Memory::Backend::CPU;
    config.method = SimType::SV;
    config.n_qubits = 1;
    auto handle = Memory::create_memory(config);
    StateView &state = handle.view;

    const ValType x_real[4] = {0.0, 1.0, 1.0, 0.0};
    const ValType x_imag[4] = {0.0, 0.0, 0.0, 0.0};

    CPU::apply_c1_gate(x_real, x_imag, 0, state);

    constexpr double eps = 1e-12;
    assert(std::fabs(state.data_real[0]) < eps);
    assert(std::fabs(state.data_imag[0]) < eps);
    assert(std::fabs(state.data_real[1] - 1.0) < eps);
    assert(std::fabs(state.data_imag[1]) < eps);

    std::cout << "CPU gate test passed\n";
    return 0;
}
