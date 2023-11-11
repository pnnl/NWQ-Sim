// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Muqing Zheng and Ang Li
// Pacific Northwest National Laboratory(PNNL), U.S.
// GitHub repo: http://www.github.com/pnnl/SV-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
#ifndef DEVICE_NOISE_HPP_
#define DEVICE_NOISE_HPP_

#include <fstream>   // std::ifstream
#include <vector>    // std::vector, to tell if a gate is 1 or 2 qubits
#include <algorithm> // std::find
#include <complex>
#include "../nlohmann/json.hpp" // json parser
#include "noise_model.hpp"

#include "../../gate.hpp"
#include "../config.hpp"
#include "../sim_gate.hpp"

using json = nlohmann::json;

namespace NWQSim
{
    // Records all names of 2 qubit gates
    std::vector<OP> gates2q{OP::CX};

    //================ Set gates =================
    void set_X(std::complex<double> *res)
    {
        /******************************************
         * Pauli-X gate: bit-flip or NOT gate
         *  X = [0 1]
         *      [1 0]
         ******************************************/
        const std::complex<double> x_gate[2][2] = {{{0, 0}, {1, 0}},
                                                   {{1, 0}, {0, 0}}};
        std::copy(&x_gate[0][0], &x_gate[0][0] + 4, res);
    }
    void set_ID(std::complex<double> *res)
    {
        /******************************************
         * Identiy gate, this is meaningful
         * for noisy simulation
         ******************************************/
        const std::complex<double> id_gate[2][2] = {{{1, 0}, {0, 0}},
                                                    {{0, 0}, {1, 0}}};
        std::copy(&id_gate[0][0], &id_gate[0][0] + 4, res);
    }
    void set_SX(std::complex<double> *res)
    {
        /******************************************
         * sqrt(X) gate, basis gate for IBM-Q
         * SX = 1/2 [1+i 1-i]
         *          [1-i 1+i]
         ******************************************/
        const std::complex<double> sx_gate[2][2] = {{std::complex<double>(0.5, 0.5), std::complex<double>(0.5, -0.5)},
                                                    {std::complex<double>(0.5, -0.5), std::complex<double>(0.5, 0.5)}};
        std::copy(&sx_gate[0][0], &sx_gate[0][0] + 4, res);
    }
    void set_RZ(std::complex<double> *res, double theta)
    {
        /******************************************
         * Rotation around Z axis
         * RZ = [cos(a/2)-i*sin(a/2)  0]
         *      [0  cos(a/2)+i*sin(a/2)]
         ******************************************/
        const std::complex<double> rz_gate[2][2] = {
            {{cos(theta / 2), -sin(theta / 2)}, {0, 0}},
            {{0, 0}, {cos(theta / 2), sin(theta / 2)}}};
        std::copy(&rz_gate[0][0], &rz_gate[0][0] + 4, res);
    }
    void set_CX(std::complex<double> *res)
    {
        /******************************************
         * Controlled X gate (CNOT)
         * Apply X when the control qubit is 1
         ******************************************/

        const std::complex<double> cx_gate[4][4] = {{{1, 0}, {0, 0}, {0, 0}, {0, 0}},
                                                    {{0, 0}, {1, 0}, {0, 0}, {0, 0}},
                                                    {{0, 0}, {0, 0}, {0, 0}, {1, 0}},
                                                    {{0, 0}, {0, 0}, {1, 0}, {0, 0}}};
        std::copy(&cx_gate[0][0], &cx_gate[0][0] + 16, res);
    }

    // ---------------------------------------------- Convert Superoperator to Choi Matrix -------------------------------
    // I think this is not necessary since we can compare the noise channel with identity, where later one is unitary
    // So we just compute the trace of noise superoperator and divided it by dim^2 to get process fidelity

    // // Ref Eq.(4.3) in https://arxiv.org/pdf/1111.6950.pdf
    // void colReshuffle1q(int dim, int a, int b, int* new_a, int* new_b) {
    //     // To bits
    //     std::bitset<2> bit_a (a);
    //     std::bitset<2> bit_b (b);
    //     // Split to two halves
    //     std::bitset<1> a_right = std::bitset<1> ((bit_a >> 1).to_ulong());
    //     std::bitset<1> a_left = std::bitset<1> (((bit_a << 1) >> 1).to_ulong());
    //     std::bitset<1> b_right = std::bitset<1> ((bit_b >> 1).to_ulong());
    //     std::bitset<1> b_left = std::bitset<1> (((bit_b << 1) >> 1).to_ulong());
    //     // Reshuffle
    //     std::bitset<2> new_bit_a (b_left.to_string()+a_left.to_string());
    //     std::bitset<2> new_bit_b (b_right.to_string()+a_right.to_string());

    //     // Convert back to integer
    //     *new_a = static_cast<int>(new_bit_a.to_ulong());
    //     *new_b = static_cast<int>(new_bit_b.to_ulong());
    // }

    // void colReshuffle2q(int dim, int a, int b, int* new_a, int* new_b) {
    //     // To bits
    //     std::bitset<4> bit_a (a);
    //     std::bitset<4> bit_b (b);
    //     // Split to two halves
    //     std::bitset<2> a_right = std::bitset<2> ((bit_a >> 2).to_ulong());
    //     std::bitset<2> a_left = std::bitset<2> (((bit_a << 2) >> 2).to_ulong());
    //     std::bitset<2> b_right = std::bitset<2> ((bit_b >> 2).to_ulong());
    //     std::bitset<2> b_left = std::bitset<2> (((bit_b << 2) >> 2).to_ulong());
    //     // Reshuffle
    //     std::bitset<4> new_bit_a (b_left.to_string()+a_left.to_string());
    //     std::bitset<4> new_bit_b (b_right.to_string()+a_right.to_string());

    //     // Convert back to integer
    //     *new_a = static_cast<int>(new_bit_a.to_ulong());
    //     *new_b = static_cast<int>(new_bit_b.to_ulong());
    // }

    // // Superoperator to Choi
    // void sp2choi(std::complex<double>* noise_sp, int dim, std::complex<double>* noise_choi, bool if_normalize = true) {
    //     double rate = 1/dim;
    //     for (int i=0; i<dim; i++) {
    //         for (int j=0; j<dim; j++) {
    //             // Do reshuffle on coordinates
    //             int new_i;
    //             int new_j;
    //             if (dim == 4) {// 1-qubit
    //                 colReshuffle1q(4, i, j, &new_i, &new_j);
    //             } else if (dim == 16) {
    //                 colReshuffle2q(16, i, j, &new_i, &new_j);
    //             } else {
    //                 throw std::invalid_argument( "Only 1 or 2-qubit gates." );
    //             }
    //             // Do normalization if necessary
    //             if ( if_normalize) {
    //                 noise_choi[new_i*dim+new_j] = rate * noise_sp[i*dim+j];
    //             } else {
    //                 noise_choi[new_i*dim+new_j] = noise_sp[i*dim+j];
    //             }

    //         }
    //     }
    // }
    // ---------------------------------------------------------------------------------------------------------------

    /**
     * @brief Compute process fidelity of a given noisy superoperator (to the ideal situation, i.e., identity operator) in dimension dim
     *        Ref: https://qiskit.org/documentation/stubs/qiskit.quantum_info.process_fidelity.html#qiskit.quantum_info.process_fidelity
     *
     * @param noise_sp
     * @param dim
     * @return double
     */
    double procFid(std::complex<double> *noise_sp, int dim)
    {
        double trace = 0.0;
        double rate = 1.0 / dim;
        for (int i = 0; i < dim; i++)
        {
            trace += std::real(noise_sp[i * dim + i]); // the number should be real anyway
        }

        return trace * rate;
    }

    /**
     * @brief Compute average gate fidelity of a given noisy superoperator (to the ideal situation, i.e., identity operator) in dimension dim
     *        Ref: https://qiskit.org/documentation/stubs/qiskit.quantum_info.average_gate_fidelity.html#qiskit.quantum_info.average_gate_fidelity
     *           Equation: (dim*process_fidelity + 1)/(dim+1)   <---- this is a special case because we compare an operator to identity operator which is unitary
     *
     * @param noise_sp
     * @param dim
     * @return double
     */
    double aveGateFid(std::complex<double> *noise_sp, int dim)
    { // note, dim is dimension of superoperator, which is square of state_dim
        double state_dim = std::sqrt(dim);
        double process_fidelity = procFid(noise_sp, dim);
        double denom = 1.0 / (state_dim + 1.0);
        double ave_fidelity = (state_dim * process_fidelity + 1.0) * denom;

        return ave_fidelity;
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
    DMGate generateDMGate(OP gate_op, int q1, int q2 = -1, double theta = 0)
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
                std::string str_q1 = std::to_string(q1); // key must be string
                std::string str_q2 = std::to_string(q2); // key must be string

                double T1_1, T1_2, T2_1, T2_2, gate_len, err_rate;
                try
                {
                    T1_1 = Config::backend_config["T1"][str_q1];
                    T2_1 = Config::backend_config["T2"][str_q1];
                    T1_2 = Config::backend_config["T1"][str_q2];
                    T2_2 = Config::backend_config["T2"][str_q2];
                    gate_len = Config::backend_config["gate_lens"][gate_name + str_q1 + "_" + str_q2];
                    err_rate = Config::backend_config["gate_errs"][gate_name + str_q1 + "_" + str_q2];
                }
                catch (...)
                {
                    throw std::invalid_argument("2-qubit gate (%s,%s) properties is not contained in the configuration file.");
                }

                std::complex<double> tr_sp[sp_dim][sp_dim] = {};
                std::complex<double> tr_sp_noise[sp_dim][sp_dim] = {};
                addTRErr2Q(gate_len, T1_1, T2_1, gate_len, T1_2, T2_2, ideal_gate, tr_sp_noise, false);
                double tr_fid = aveGateFid(tr_sp_noise[0], 16);
                double tr_infid = 1.0 - tr_fid;
                addTRErr2Q(gate_len, T1_1, T2_1, gate_len, T1_2, T2_2, ideal_gate, tr_sp, true);
                // Depolarizing Error
                if (err_rate <= tr_infid)
                {
                    // If relaxation error is already too large, skip depolarizing error
                    std::copy(&tr_sp[0][0], &tr_sp[0][0] + sp_dim * sp_dim, gate_sp[0]);
                }
                else
                {                                                                      // add depolarizing error
                    double err_rate_max = qubit_dim_double / (qubit_dim_double + 1.0); // 2^n/(2^n+1), n = 1
                    double err_rate_fixed = err_rate;
                    if (err_rate > err_rate_max)
                    { // if when error rate is too large, I think this actually means MODEL FAILURE
                        err_rate_fixed = err_rate_max;
                    }
                    double dep_rate = qubit_dim_double * (err_rate_fixed - tr_infid) / (qubit_dim_double * tr_fid - 1.0);
                    double dep_rate_max = min(1.0, sp_dim_double / (sp_dim_double - 1.0)); // maximum depolarizing rate (note this is not the error probability)
                    double dep_rate_fixed = dep_rate;
                    if (dep_rate > dep_rate_max)
                    { // Again, I think this actually means MODEL FAILURE
                        dep_rate_fixed = dep_rate_max;
                    }
                    // Now we construct depolarizing error
                    std::complex<double> dep_sp[sp_dim][sp_dim] = {};
                    addDepErr2Q(dep_rate_fixed, ideal_gate, dep_sp, false);
                    dotProd(tr_sp[0], dep_sp[0], gate_sp[0], sp_dim);
                }
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
                throw std::invalid_argument("Unsupported basis gate!");
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
                std::string str_q1 = std::to_string(q1); // key must be string
                double T1, T2; 
                try
                {
                    T1 = Config::backend_config["T1"][str_q1];
                    T2 = Config::backend_config["T2"][str_q1];
                }
                catch (...)
                {
                    throw std::invalid_argument("1-qubit gate properties is not contained in the configuration file.");
                }
                std::complex<double> tr_sp[sp_dim][sp_dim] = {};

                if (gate_op == OP::DELAY)
                {
                    //Delay gate
                    addTRErr1Q(theta, T1, T2, ideal_gate, tr_sp, false);
                    std::copy(&tr_sp[0][0], &tr_sp[0][0] + sp_dim * sp_dim, gate_sp[0]);
                }
                else
                {
                    double gate_len, err_rate;
                    try 
                    {
                        gate_len = Config::backend_config["gate_lens"][gate_name + str_q1];
                        err_rate = Config::backend_config["gate_errs"][gate_name + str_q1];
                    }
                    catch (...)
                    {
                        throw std::invalid_argument("1-qubit gate properties is not contained in the configuration file.");
                    }
                    //Other 1-qubit gates
                    std::complex<double> tr_sp_noise[sp_dim][sp_dim] = {};
                    addTRErr1Q(gate_len, T1, T2, ideal_gate, tr_sp_noise, false);
                    double tr_fid = aveGateFid(tr_sp_noise[0], sp_dim);
                    double tr_infid = 1.0 - tr_fid;
                    addTRErr1Q(gate_len, T1, T2, ideal_gate, tr_sp, true);
                    // Depolarizing Error
                    if (err_rate <= tr_infid)
                    {
                        // If relaxation error is already too large, skip depolarizing error
                        std::copy(&tr_sp[0][0], &tr_sp[0][0] + sp_dim * sp_dim, gate_sp[0]);
                    }
                    else
                    {                                                                      // add depolarizing error
                        double err_rate_max = qubit_dim_double / (qubit_dim_double + 1.0); // 2^n/(2^n+1), n = 1
                        double err_rate_fixed = err_rate;

                        if (err_rate > err_rate_max)
                        { // if when error rate is too large, I think this actually means MODEL FAILURE
                            err_rate_fixed = err_rate_max;
                        }
                        double dep_rate = qubit_dim_double * (err_rate_fixed - tr_infid) / (qubit_dim_double * tr_fid - 1.0);
                        double dep_rate_max = min(1.0, sp_dim_double / (sp_dim_double - 1.0)); // maximum depolarizing rate (note this is not the error probability)
                        double dep_rate_fixed = dep_rate;

                        if (dep_rate > dep_rate_max)
                        { // Again, I think this actually means MODEL FAILURE
                            dep_rate_fixed = dep_rate_max;
                        }
                        // Now we construct depolarizing error
                        std::complex<double> dep_sp[sp_dim][sp_dim] = {};
                        addDepErr1Q(dep_rate_fixed, ideal_gate, dep_sp, false);
                        dotProd(tr_sp[0], dep_sp[0], gate_sp[0], sp_dim);
                    }
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
     * @param backend_config nlohmann::json object, stores information for backend noise parameters, can be load from function readConfigFile()
     * @param noise_sp where the noise superoperator is saved
     * @param qubit_index the measured qubit
     * @param relaxation_noise if the relaxation noise should be added. Defult is false, since the parameter estimation usually already consider this part
     */
    void getMeasureSP(std::complex<double> *gate_sp,
                      int qubit_index,
                      bool relaxation_noise = false)
    {
        if (qubit_index < 0)
        {
            throw std::invalid_argument("Qubit index must be non-negative.");
        }
        // Read qubits properties
        std::string str_q = std::to_string(qubit_index); // key must be string
        double T1, T2, rd_len, p_1g0, p_0g1;
        try
        {
            T1 = Config::backend_config["T1"][str_q];
            T2 = Config::backend_config["T2"][str_q];
            rd_len = Config::backend_config["readout_length"][str_q];
            p_1g0 = Config::backend_config["prob_meas1_prep0"][str_q];
            p_0g1 = Config::backend_config["prob_meas0_prep1"][str_q];
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

} // end of namespace DMSim
#endif
