#include "../nwq_api.hpp"
#include "../../single/svsim/src/util.h"
#include "../../single/svsim/src/svsim_nvgpu_sin.cuh"

#include <map>
#include <string>
#include <functional>

namespace NWQSim
{
    class NWQBackendSVNVGPU : public NWQBackend
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
                {"U", [&]()
                 { m_sim->U(params[0], params[1], params[2], qubits[0]); }},
                {"U1", [&]()
                 { m_sim->U1(params[0], qubits[0]); }},
                {"U2", [&]()
                 { m_sim->U2(params[0], params[1], qubits[0]); }},
                {"U3", [&]()
                 { m_sim->U3(params[0], params[1], params[2], qubits[0]); }},
                {"X", [&]()
                 { m_sim->X(qubits[0]); }},
                {"Y", [&]()
                 { m_sim->Y(qubits[0]); }},
                {"Z", [&]()
                 { m_sim->Z(qubits[0]); }},
                {"H", [&]()
                 { m_sim->H(qubits[0]); }},
                {"S", [&]()
                 { m_sim->S(qubits[0]); }},
                {"SDG", [&]()
                 { m_sim->SDG(qubits[0]); }},
                {"T", [&]()
                 { m_sim->T(qubits[0]); }},
                {"TDG", [&]()
                 { m_sim->TDG(qubits[0]); }},
                {"RX", [&]()
                 { m_sim->RX(params[0], qubits[0]); }},
                {"RY", [&]()
                 { m_sim->RY(params[0], qubits[0]); }},
                {"RZ", [&]()
                 { m_sim->RZ(params[0], qubits[0]); }},
                {"CX", [&]()
                 { m_sim->CX(qubits[0], qubits[1]); }},
                {"CY", [&]()
                 { m_sim->CY(qubits[0], qubits[1]); }},
                {"CZ", [&]()
                 { m_sim->CZ(qubits[0], qubits[1]); }},
                {"CH", [&]()
                 { m_sim->CH(qubits[0], qubits[1]); }},
                {"CCX", [&]()
                 { m_sim->CCX(qubits[0], qubits[1], qubits[2]); }},
                {"CRX", [&]()
                 { m_sim->CRX(params[0], qubits[0], qubits[1]); }},
                {"CRY", [&]()
                 { m_sim->CRY(params[0], qubits[0], qubits[1]); }},
                {"CRZ", [&]()
                 { m_sim->CRZ(params[0], qubits[0], qubits[1]); }},
                {"CU", [&]()
                 { m_sim->CU(params[0], params[1], params[2], params[3], qubits[0], qubits[1]); }},
                {"CU1", [&]()
                 { m_sim->CU(0, 0, params[0], 0, qubits[0], qubits[1]); }},
                {"CU3", [&]()
                 { m_sim->CU(params[0], params[1], params[2], 0, qubits[0], qubits[1]); }},
                {"RESET", [&]()
                 { m_sim->RESET(qubits[0]); }},
                {"SWAP", [&]()
                 { m_sim->SWAP(qubits[0], qubits[1]); }},
                {"RI", [&]()
                 { m_sim->RI(params[0], qubits[0]); }},
                {"P", [&]()
                 { m_sim->P(params[0], qubits[0]); }},
                {"CS", [&]()
                 { m_sim->CS(qubits[0], qubits[1]); }},
                {"CSDG", [&]()
                 { m_sim->CSDG(qubits[0], qubits[1]); }},
                {"CT", [&]()
                 { m_sim->CT(qubits[0], qubits[1]); }},
                {"CTDG", [&]()
                 { m_sim->CTDG(qubits[0], qubits[1]); }},
                {"CSX", [&]()
                 { m_sim->CSX(qubits[0], qubits[1]); }},
                {"CP", [&]()
                 { m_sim->CP(params[0], qubits[0], qubits[1]); }},
                {"CSWAP", [&]()
                 { m_sim->CSWAP(qubits[0], qubits[1], qubits[2]); }},
                {"ID", [&]()
                 { m_sim->ID(qubits[0]); }},
                {"RXX", [&]()
                 { m_sim->RXX(params[0], qubits[0], qubits[1]); }},
                {"RYY", [&]()
                 { m_sim->RYY(params[0], qubits[0], qubits[1]); }},
                {"RZZ", [&]()
                 { m_sim->RZZ(params[0], qubits[0], qubits[1]); }},
            };

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

        virtual double calcExpectationValueZ(const std::vector<size_t> &in_bits) override
        {

            const auto hasEvenParity = [](unsigned long long x, const std::vector<size_t> &in_qubitIndices) -> bool
            {
                size_t count = 0;
                for (const auto &bitIdx : in_qubitIndices)
                {
                    if (x & (1ULL << bitIdx))
                    {
                        count++;
                    }
                }
                return (count % 2) == 0;
            };

            m_sim->sim();

            cudaSafeCall(cudaMemcpy(m_sim->sv_real_cpu, m_sim->sv_real, m_sim->sv_size, cudaMemcpyDeviceToHost));

            cudaSafeCall(cudaMemcpy(m_sim->sv_imag_cpu, m_sim->sv_imag, m_sim->sv_size, cudaMemcpyDeviceToHost));

            double result = 0.0;

            for (unsigned long long i = 0; i < m_sim->dim; ++i)
            {
                result += (hasEvenParity(i, in_bits) ? 1.0 : -1.0) *
                          (m_sim->sv_real_cpu[i] * m_sim->sv_real_cpu[i] + m_sim->sv_imag_cpu[i] * m_sim->sv_imag_cpu[i]);
            }

            return result;
        }

    private:
        std::shared_ptr<::NWQSim::Simulation> m_sim;
    };

    std::shared_ptr<NWQBackend> get_backend()
    {
        return std::make_shared<NWQBackendSVNVGPU>();
    }

}

/// hpc/home/lian599/svsim_qsharp_noise