#include "../nwq_api.hpp"
#include "../../single/dmsim/src/util.h"
#include "../../single/dmsim/src/dmsim_cpu_sin.hpp"

#include <map>
#include <string>
#include <functional>

namespace NWQSim
{
    class NWQBackendDMCpuSin : public NWQBackend
    {
    public:
        virtual void init(int n_qubits) override
        {
            m_sim = std::make_shared<::NWQSim::Simulation>(n_qubits);
        }
        virtual void add_gate(std::string gate_name, const std::vector<int> &qubits,
                              const std::vector<double> &params) override
        {
            std::map<std::string, std::function<void()>> op_map = {
                {"X", [&]()
                 { m_sim->X(qubits[0]); }},
                {"RZ", [&]()
                 { m_sim->RZ(params[0], qubits[0]); }},
                {"CX", [&]()
                 { m_sim->CX(qubits[0], qubits[1]); }},
                {"RESET", [&]()
                 { m_sim->RESET(qubits[0]); }}};

            // Assuming 'op' is of type std::string
            auto it = op_map.find(gate_name);
            if (it != op_map.end())
            {
                it->second();
            }
            else
            {
                throw std::logic_error("Undefined gate is called!");
            }
        }

        virtual std::vector<int64_t> measure(int shots) override
        {
            // Run the simulation
            auto *res = m_sim->measure_all(shots);

            std::vector<int64_t> result;
            result.reserve(shots);
            for (int i = 0; i < shots; ++i)
            {
                result.emplace_back(res[i]);
            }

            return result;
        }

        virtual void finalize() override
        {
            m_sim->reset_sim();
        }

    private:
        std::shared_ptr<::NWQSim::Simulation> m_sim;
    };

    std::shared_ptr<NWQBackend> get_backend()
    {
        return std::make_shared<NWQBackendDMCpuSin>();
    }

}

/// hpc/home/lian599/svsim_qsharp_noise