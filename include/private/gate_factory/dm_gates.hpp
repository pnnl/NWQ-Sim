#pragma once

#include <vector>
#include <complex>

#include "device_noise.hpp"
#include "../sim_gate.hpp"
#include "../config.hpp"

namespace NWQSim
{

    std::vector<DMGate> getDMGates(const std::vector<Gate> &gates, const IdxType n_qubits)
    {
        std::vector<DMGate> sim_dm_gates;

        for (const auto &g : gates)
        {

            if (g.op_name == OP::RESET)
            {
                sim_dm_gates.push_back(DMGate(OP::RESET, g.qubit, g.ctrl));
            }
            else if (g.op_name == OP::M)
            {
                if (Config::ENABLE_NOISE)
                {
                    std::complex<double> noisy_operator[4][4] = {};
                    getMeasureSP(noisy_operator[0], g.qubit);

                    DMGate noisy_dm_gate(OP::C2, g.qubit, g.ctrl);
                    noisy_dm_gate.set_gm(noisy_operator[0], 4);

                    sim_dm_gates.push_back(noisy_dm_gate);
                }

                sim_dm_gates.push_back(DMGate(OP::M, g.qubit, g.ctrl));
            }
            else if (g.op_name == OP::MA)
            {
                if (Config::ENABLE_NOISE)
                {
                    for (IdxType i = 0; i < n_qubits; i++)
                    {
                        std::complex<double> noisy_operator[4][4] = {};
                        getMeasureSP(noisy_operator[0], i);

                        DMGate noisy_dm_gate(OP::C2, i, g.ctrl);
                        noisy_dm_gate.set_gm(noisy_operator[0], 4);

                        sim_dm_gates.push_back(noisy_dm_gate);
                    }
                }
                sim_dm_gates.push_back(DMGate(OP::MA, g.repetition, g.ctrl));
            }
            else
            {
                sim_dm_gates.push_back(generateDMGate(g.op_name, g.qubit, g.ctrl, g.theta));
            }
        }
        return sim_dm_gates;
    }
}
