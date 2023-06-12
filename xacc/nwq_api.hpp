#pragma once
#include <vector>
#include <memory>
#include <string>

namespace NWQSim
{

    class NWQBackend
    {
    public:
        virtual void init(int n_qubits) = 0;
        virtual void add_gate(std::string gate_name, const std::vector<int> &qubits,
                              const std::vector<double> &params = {}) = 0;
        virtual std::vector<int64_t> measure(int shots) = 0;
        virtual void finalize() = 0;

        virtual double calcExpectationValueZ(const std::vector<size_t> &in_bits) = 0;

        
    };

    extern std::shared_ptr<NWQBackend> get_backend();
} // namespace