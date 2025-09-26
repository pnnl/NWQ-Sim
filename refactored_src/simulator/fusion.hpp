#pragma once

#include "memory/memory_handle.hpp"

#include <vector>

namespace NWQSim
{
    std::vector<DeviceGate> fuse_device_gates(const std::vector<DeviceGate> &input,
                                              IdxType n_qubits);
}
