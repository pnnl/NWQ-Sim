#include <cassert>
#include <cmath>
#include <iostream>

#include "simulator/sv_sim.hpp"

int main()
{
    using namespace NWQSim;

    constexpr ValType eps = 1e-12;

    // Simple single-qubit superposition
    Circuit single_qubit(1);
    single_qubit.h(0);

    SVSim single_sim(1, Memory::Backend::CPU, 1234);
    auto single_result = single_sim.simulate(single_qubit);
    const auto &single_state = single_sim.state();

    assert(single_result.single_qubit_measurements.empty());
    assert(single_result.measure_all_results.empty());
    assert(std::fabs(single_state.real[0] - S2I) < eps);
    assert(std::fabs(single_state.imag[0]) < eps);
    assert(std::fabs(single_state.real[1] - S2I) < eps);
    assert(std::fabs(single_state.imag[1]) < eps);

    // Two-qubit Bell state
    Circuit bell(2);
    bell.h(0);
    bell.cx(0, 1);

    SVSim bell_sim(2, Memory::Backend::CPU, 42);
    auto bell_result = bell_sim.simulate(bell);
    const auto &bell_state = bell_sim.state();

    assert(bell_result.single_qubit_measurements.empty());
    assert(std::fabs(bell_state.real[0] - S2I) < eps);
    assert(std::fabs(bell_state.imag[0]) < eps);
    assert(std::fabs(bell_state.real[3] - S2I) < eps);
    assert(std::fabs(bell_state.imag[3]) < eps);
    assert(std::fabs(bell_state.real[1]) < eps);
    assert(std::fabs(bell_state.imag[1]) < eps);
    assert(std::fabs(bell_state.real[2]) < eps);
    assert(std::fabs(bell_state.imag[2]) < eps);

    // Measurement on |0> should deterministically return 0
    Circuit measure_zero(1);
    measure_zero.measure(0);
    SVSim measure_sim(1, Memory::Backend::CPU, 99);
    auto measure_result = measure_sim.simulate(measure_zero);
    assert(measure_result.single_qubit_measurements.size() == 1);
    assert(measure_result.single_qubit_measurements.front() == 0);

    // Fusion of sequential single-qubit gates should condense to one device gate
    Circuit double_h(1);
    double_h.h(0);
    double_h.h(0);
    SVSim fusion_single_sim(1, Memory::Backend::CPU, 17);
    auto fusion_single_result = fusion_single_sim.simulate(double_h);
    const auto &fusion_single_state = fusion_single_sim.state();
    assert(fusion_single_result.metadata.total_device_gates == 1);
    assert(std::fabs(fusion_single_state.real[0] - 1.0) < eps);
    assert(std::fabs(fusion_single_state.imag[0]) < eps);
    assert(std::fabs(fusion_single_state.real[1]) < eps);
    assert(std::fabs(fusion_single_state.imag[1]) < eps);

    SVSim fusion_single_disabled(1, Memory::Backend::CPU, 33);
    fusion_single_disabled.set_fusion_enabled(false);
    auto fusion_single_disabled_result = fusion_single_disabled.simulate(double_h);
    assert(fusion_single_disabled_result.metadata.total_device_gates == 2);

    // Fusion of repeated two-qubit gates on the same pair
    Circuit double_cx(2);
    double_cx.cx(0, 1);
    double_cx.cx(0, 1);
    SVSim fusion_two_sim(2, Memory::Backend::CPU, 21);
    auto fusion_two_result = fusion_two_sim.simulate(double_cx);
    const auto &fusion_two_state = fusion_two_sim.state();
    assert(fusion_two_result.metadata.total_device_gates == 1);
    assert(std::fabs(fusion_two_state.real[0] - 1.0) < eps);
    assert(std::fabs(fusion_two_state.imag[0]) < eps);

    SVSim fusion_two_disabled(2, Memory::Backend::CPU, 45);
    fusion_two_disabled.set_fusion_enabled(false);
    auto fusion_two_disabled_result = fusion_two_disabled.simulate(double_cx);
    assert(fusion_two_disabled_result.metadata.total_device_gates == 2);

    Circuit param(1);
    auto theta_idx = param.declare_parameter(PI);
    param.rx_param(0, theta_idx);
    SVSim param_sim(1, Memory::Backend::CPU, 12);
    (void)param_sim.simulate(param);
    const auto &param_state = param_sim.state();
    assert(std::fabs(param_state.real[0]) < eps);
    assert(std::fabs(param_state.real[1]) < eps);
    assert(std::fabs(param_state.imag[1] + 1.0) < eps);

    param.set_parameter(theta_idx, PI / 2);
    (void)param_sim.simulate(param);
    const auto &updated_state = param_sim.state();
    assert(std::fabs(updated_state.real[0] - std::cos(PI / 4)) < eps);

    param.set_parameters({0.0});
    (void)param_sim.simulate(param);
    const auto &reset_state = param_sim.state();
    assert(std::fabs(reset_state.real[0] - 1.0) < eps);

    std::cout << "Refactored SV simulator tests passed\n";
    return 0;
}
