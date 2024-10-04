// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Muqing Zheng and Ang Li
// Pacific Northwest National Laboratory(PNNL), U.S.
// GitHub repo: http://www.github.com/pnnl/SV-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
#pragma once

#include <fstream>   // std::ifstream
#include <vector>    // std::vector, to tell if a gate is 1 or 2 qubits
#include <algorithm> // std::find
#include <complex>

#include "noise_util.hpp"
#include "noise_model.hpp"
#include "../sim_gate.hpp"

#include "../../gate.hpp"
#include "../../config.hpp"

namespace NWQSim
{
    DMGate generateDMGate(NoiseModel noise_model, OP gate_op, int q1, int q2, double theta);
    void getMeasureSP(NoiseModel noise_model, std::complex<double> *gate_sp, int qubit_index, bool relaxation_noise = false);

    std::vector<DMGate> getDMGates(const std::vector<Gate> &gates, const IdxType n_qubits)
    {
        std::vector<DMGate> sim_dm_gates;

        NoiseModel noise_model(Config::device_noise_file, Config::device_layout_file, Config::device_layout_str);

        if (n_qubits > noise_model.get_num_qubits())
        {
            std::string msg = "Error: Circuit uses " + std::to_string(n_qubits) + " qubits, more than " + std::to_string(noise_model.get_num_qubits()) + "qubits in the device!!\n";
            throw std::logic_error(msg.c_str());
        }

        for (const auto &g : gates)
        {
            if (g.op_name == MOD_NOISE)
            {
                noise_model.modify_noise(g.mod_op, g.mod_noise, g.mod_value, g.mod_qubits);
            }
            else if (g.op_name == OP::RESET)
            {
                sim_dm_gates.push_back(DMGate(OP::RESET, g.qubit, g.ctrl));
            }
            else if (g.op_name == OP::M)
            {
                if (Config::ENABLE_NOISE)
                {
                    std::complex<double> noisy_operator[4][4] = {};
                    getMeasureSP(noise_model, noisy_operator[0], g.qubit);

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
                        getMeasureSP(noise_model, noisy_operator[0], i);

                        DMGate noisy_dm_gate(OP::C2, i, g.ctrl);
                        noisy_dm_gate.set_gm(noisy_operator[0], 4);

                        sim_dm_gates.push_back(noisy_dm_gate);
                    }
                }
                sim_dm_gates.push_back(DMGate(OP::MA, g.repetition, g.ctrl));
            }
            else
            {
                sim_dm_gates.push_back(generateDMGate(noise_model, g.op_name, g.qubit, g.ctrl, g.theta));
            }
        }
        return sim_dm_gates;
    }

    // Create gate error superoperator
    // Ref https://github.com/Qiskit/qiskit-aer/blob/f778f4b53b4223c86f5e41b50794d176208c76c6/qiskit/providers/aer/noise/device/models.py#L171
    /**
     * @brief Computer a noise superoperator for a given basis gate on the given qubit(s). The noise parameters are read from a nlohmann::json object
     *        NOTE: if a gate wants to be recognized as 2-qubit gate, please add the name to std::vector<std::string> gates2q at line 12.
     *
     * @param gate_name gate name in string. It must be the name appear in the backend_config.For IBMQ backend, it usually be one of {"id", "sx", "rz", "x", "reset", "cx"}
     * @param noise_sp where the noise superoperator is saved
     * @param q1 the (first) qubit that the gate applied on. For 1-qubit gate, just fill this parameter and ignore the next parameter.
     * @param q2 the second qubit that a two-qubit gate applied on. It is only used for a two-qubit gate.
     */
    DMGate generateDMGate(NoiseModel noise_model, OP gate_op, int q1, int q2, double theta)
    {
        std::string gate_name(OP_NAMES[gate_op]);
        std::transform(gate_name.begin(), gate_name.end(), gate_name.begin(), [](unsigned char c)
                       { return std::tolower(c); });

        if (q1 < 0)
        {
            throw std::invalid_argument("Check your first qubit. Qubit index must be non-negative.");
        }
        // Build gate noise
        if (std::find(gates2q.begin(), gates2q.end(), gate_op) != gates2q.end())
        { // if 2-qubit gate
            if (q2 < 0)
            {
                throw std::invalid_argument("Check your second qubit. Qubit index must be non-negative.");
            }
            const int qubit_dim = 4;
            // const int sp_dim = qubit_dim*qubit_dim;
            std::complex<double> ideal_gate[qubit_dim][qubit_dim] = {};
            set_CX(&ideal_gate[0][0]);

            std::complex<double> gate_sp[16][16] = {};

            if (!Config::ENABLE_NOISE)
            {
                std::complex<double> ideal_gate_conj[qubit_dim][qubit_dim] = {}; // conjugate of ideal gate
                conjbar(ideal_gate[0], ideal_gate_conj[0], qubit_dim);
                // std::complex<double> ideal_gate_sp[sp_dim][sp_dim] = {}; // super operator of ideal gate
                kronProd(ideal_gate_conj[0], ideal_gate[0], 1.0, gate_sp[0], qubit_dim);
            }
            else
            {
                const int sp_dim = qubit_dim * qubit_dim;
                double qubit_dim_double = (double)qubit_dim;
                double sp_dim_double = (double)sp_dim;
                // Read qubits properties
                double T1_1, T1_2, T2_1, T2_2, gate_len, err_rate;
                try
                {
                    T1_1 = noise_model.get_qubit_parameter(T1_KEY, q1);
                    T2_1 = noise_model.get_qubit_parameter(T2_KEY, q1);

                    T1_2 = noise_model.get_qubit_parameter(T1_KEY, q2);
                    T2_2 = noise_model.get_qubit_parameter(T2_KEY, q2);

                    gate_len = noise_model.get_gate_parameter(GATE_LENS_KEY, gate_name, q1, q2);
                    err_rate = noise_model.get_gate_parameter(GATE_ERRS_KEY, gate_name, q1, q2);
                }
                catch (...)
                {
                    // std::strin
                    throw std::invalid_argument("2-qubit gate: " +
                                                std::string(OP_NAMES[gate_op]) + std::to_string(q1) + std::to_string(q2) +
                                                " properties is not contained in the configuration file.");
                }

                std::complex<double> tr_sp[sp_dim][sp_dim] = {};
                addTRErr2Q(gate_len, T1_1, T2_1, gate_len, T1_2, T2_2, ideal_gate, tr_sp, true);

                // Now we construct depolarizing error
                std::complex<double> dep_sp[sp_dim][sp_dim] = {};
                addDepErr2Q(err_rate, ideal_gate, dep_sp, false);
                dotProd(tr_sp[0], dep_sp[0], gate_sp[0], sp_dim);
            }

            DMGate gate(OP::C4, q1, q2);
            gate.set_gm(gate_sp[0], 16);
            return gate;
        }
        else
        { // 1-qubit gate, which is basically the same code as in the previous if. There should be a smarter way.

            std::complex<double> gate_sp[4][4] = {};

            const int qubit_dim = 2;
            // const int sp_dim = qubit_dim*qubit_dim;
            std::complex<double> ideal_gate[qubit_dim][qubit_dim] = {};
            switch (gate_op)
            {
            case OP::X:
                set_X(&ideal_gate[0][0]);
                break;
            case OP::ID:
                set_ID(&ideal_gate[0][0]);
                break;
            case OP::DELAY:
                set_ID(&ideal_gate[0][0]);
                break;
            case OP::SX:
                set_SX(&ideal_gate[0][0]);
                break;
            case OP::RZ:
                set_RZ(&ideal_gate[0][0], theta);
                break;
            default:
                throw std::invalid_argument("Unsupported basis gate: " + std::string(OP_NAMES[gate_op]));
            }

            if (!Config::ENABLE_NOISE)
            {
                std::complex<double> ideal_gate_conj[qubit_dim][qubit_dim] = {}; // conjugate of ideal gate
                conjbar(ideal_gate[0], ideal_gate_conj[0], qubit_dim);
                // std::complex<double> ideal_gate_sp[sp_dim][sp_dim] = {}; // super operator of ideal gate
                kronProd(ideal_gate_conj[0], ideal_gate[0], 1.0, gate_sp[0], qubit_dim);
            }
            else
            {
                const int sp_dim = qubit_dim * qubit_dim;
                double qubit_dim_double = (double)qubit_dim;
                double sp_dim_double = (double)sp_dim;
                // Read qubits properties
                double T1, T2;
                try
                {
                    T1 = noise_model.get_qubit_parameter(T1_KEY, q1);
                    T2 = noise_model.get_qubit_parameter(T2_KEY, q1);
                }
                catch (...)
                {
                    throw std::invalid_argument("1-qubit gate properties is not contained in the configuration file.");
                }
                std::complex<double> tr_sp[sp_dim][sp_dim] = {};

                if (gate_op == OP::DELAY)
                {
                    // Delay gate
                    addTRErr1Q(theta, T1, T2, ideal_gate, tr_sp, false);
                    std::copy(&tr_sp[0][0], &tr_sp[0][0] + sp_dim * sp_dim, gate_sp[0]);
                }
                else
                {
                    double gate_len, err_rate;
                    try
                    {
                        gate_len = noise_model.get_gate_parameter(GATE_LENS_KEY, gate_name, q1);
                        err_rate = noise_model.get_gate_parameter(GATE_ERRS_KEY, gate_name, q1);
                    }
                    catch (...)
                    {
                        throw std::invalid_argument("1-qubit gate properties is not contained in the configuration file.");
                    }
                    // Other 1-qubit gates
                    addTRErr1Q(gate_len, T1, T2, ideal_gate, tr_sp, true);

                    // Now we construct depolarizing error
                    std::complex<double> dep_sp[sp_dim][sp_dim] = {};
                    addDepErr1Q(err_rate, ideal_gate, dep_sp, false);
                    dotProd(tr_sp[0], dep_sp[0], gate_sp[0], sp_dim);
                }
            }
            DMGate gate(OP::C2, q1, q2);
            gate.set_gm(gate_sp[0], 4);
            return gate;
        }
    }

    // Create measurement error superoperator
    /**
     * @brief Computer a measurement noise superoperator for a given qubit from a backend noise configuratio nfile
     *
     * @param noise_sp where the noise superoperator is saved
     * @param qubit_index the measured qubit
     * @param relaxation_noise if the relaxation noise should be added. Defult is false, since the parameter estimation usually already consider this part
     */
    void getMeasureSP(NoiseModel noise_model,
                      std::complex<double> *gate_sp,
                      int qubit_index,
                      bool relaxation_noise)
    {
        if (qubit_index < 0)
        {
            throw std::invalid_argument("Qubit index must be non-negative.");
        }
        // Read qubits properties
        double T1, T2, rd_len, p_1g0, p_0g1;
        try
        {
            T1 = noise_model.get_qubit_parameter(T1_KEY, qubit_index);
            T2 = noise_model.get_qubit_parameter(T2_KEY, qubit_index);
            rd_len = noise_model.get_qubit_parameter(READOUT_LENGTH_KEY, qubit_index);
            p_1g0 = noise_model.get_qubit_parameter(PROB_MEAS1_PREP0_KEY, qubit_index);
            p_0g1 = noise_model.get_qubit_parameter(PROB_MEAS0_PREP1_KEY, qubit_index);
        }
        catch (...)
        {
            throw std::invalid_argument("Qubit or gate properties is not contained in the configuration file.");
        }

        std::complex<double>(*noise_sp_mat)[4] = (std::complex<double>(*)[4])gate_sp;
        if (relaxation_noise)
        {
            addMeaErr1Q(p_1g0, p_0g1, noise_sp_mat, rd_len, T1, T2);
        }
        else
        {
            addMeaErr1Q(p_1g0, p_0g1, noise_sp_mat);
        }
    }
}
