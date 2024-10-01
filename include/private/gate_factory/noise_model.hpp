#ifndef NOISE_MODEL_HPP
#define NOISE_MODEL_HPP

#include <string>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip> // For setting precision
#include <cmath>   // For std::abs
#include <regex>   // For regex to clean up the scientific notation string

#include "../nlohmann/json.hpp"
#include "../../nwq_util.hpp"

namespace NWQSim
{
    class NoiseModel
    {

    public:
        NoiseModel(const std::string &device_filename, const std::string &layout_filename, const std::string &layout_str)
        {
            std::ifstream f(device_filename);
            if (f.fail())
                throw std::logic_error("Device config file not found at " + device_filename);
            backend_config = nlohmann::json::parse(f);

            if (!layout_str.empty())
            {
                loadLayoutString(layout_str);
            }
            else if (!layout_filename.empty())
            {
                loadLayoutFile(layout_filename);
            }
            else
            {
                IdxType n_qubits;
                backend_config["num_qubits"].get_to(n_qubits);
                for (IdxType n = 0; n < n_qubits; n++)
                {
                    layout[std::to_string(n)] = (int)n;
                }
            }

            f.close();

            default_backend_config = backend_config;
            default_layout = layout;
        }

        IdxType get_num_qubits()
        {
            return backend_config["num_qubits"];
        }

        double get_t1(IdxType q)
        {
            return backend_config["T1"][qindex(q)];
        }

        double get_t2(IdxType q)
        {
            return backend_config["T2"][qindex(q)];
        }

        double get_gate_len(const std::string &gate_name, IdxType q1, IdxType q2 = -1)
        {

            return backend_config["gate_lens"][gname(gate_name, q1, q2)];
        }

        double get_gate_err(const std::string &gate_name, IdxType q1, IdxType q2 = -1)
        {

            return backend_config["gate_errs"][gname(gate_name, q1, q2)];
        }

        double get_readout_length(IdxType q)
        {
            return backend_config["readout_length"][qindex(q)];
        }

        double get_readout_m1p0(IdxType q)
        {
            return backend_config["prob_meas1_prep0"][qindex(q)];
        }

        double get_readout_m0p1(IdxType q)
        {
            return backend_config["prob_meas0_prep1"][qindex(q)];
        }

        void modify_noise(std::string mod_op, std::string mod_noise, ValType value, std::vector<IdxType> qubit_list)
        {
            // convert mod_op and mod_noise to lowercase
            std::transform(mod_op.begin(), mod_op.end(), mod_op.begin(), ::tolower);
            std::transform(mod_noise.begin(), mod_noise.end(), mod_noise.begin(), ::tolower);

            std::string keyword;     // This will hold the relevant keyword for backend_config
            bool is_gate = false;    // Whether the modification is for a gate or not
            bool is_2q_gate = false; // Whether the gate is a 2-qubit gate or not
            std::string gate_name;   // To store the gate name for gate-related noise

            // Determine the correct keyword based on mod_noise
            if (mod_noise == "t1")
            {
                keyword = "T1";
            }
            else if (mod_noise == "t2")
            {
                keyword = "T2";
            }
            else if (mod_noise == "readout_len")
            {
                keyword = "readout_length";
            }
            else if (mod_noise == "readout_m0p1")
            {
                keyword = "prob_meas0_prep1";
            }
            else if (mod_noise == "readout_m1p0")
            {
                keyword = "prob_meas1_prep0";
            }
            else if (mod_noise.find("_len") != std::string::npos)
            {
                // Extract gate name before "_len" and set keyword to "gate_lens"
                gate_name = mod_noise.substr(0, mod_noise.find("_len"));
                keyword = "gate_lens";
                is_gate = true;
            }
            else if (mod_noise.find("_err") != std::string::npos)
            {
                // Extract gate name before "_err" and set keyword to "gate_errs"
                gate_name = mod_noise.substr(0, mod_noise.find("_err"));
                keyword = "gate_errs";
                is_gate = true;
            }
            else
            {
                std::cerr << "Unknown noise type: " << mod_noise << std::endl;
                return;
            }

            if (is_gate) // gate err/len
            {
                is_2q_gate = gate_name == "cx" || gate_name == "ecr"; // Check if the gate is a 2-qubit gate
                if (is_2q_gate)
                    assert(qubit_list.size() == 2);
            }

            // Update the noise model, 3 cases: two-qubit gate, single-qubit gate, and single-qubit property (T1, T2, READOUT, etc.)

            if (is_gate)
            {
                if (is_2q_gate) // case 1: two-qubit gate
                {
                    // Handle two-qubit gate operation (e.g., cx_len, cx_err)
                    IdxType q1 = qubit_list[0];
                    IdxType q2 = qubit_list[1];
                    std::string sub_keyword = gname(gate_name, q1, q2); // Generate index for gate

                    // Perform the operation
                    if (mod_op == "scale")
                    {
                        backend_config[keyword][sub_keyword] = static_cast<ValType>(backend_config[keyword][sub_keyword]) * value;
                    }
                    else if (mod_op == "set")
                    {
                        backend_config[keyword][sub_keyword] = value;
                    }
                    else if (mod_op == "reset")
                    {
                        backend_config[keyword][sub_keyword] = default_backend_config[keyword][sub_keyword];
                    }
                    else
                    {
                        std::cerr << "Unknown two-qubit gate type: " << mod_op << std::endl;
                    }
                }
                else // case 2: single-qubit gate
                {
                    // loop over all qubits in qubit_list
                    for (auto qubit : qubit_list)
                    {
                        std::string sub_keyword = gname(gate_name, qubit); // Generate gate name for single qubit

                        // Perform the operation
                        if (mod_op == "scale")
                        {
                            backend_config[keyword][sub_keyword] = static_cast<ValType>(backend_config[keyword][sub_keyword]) * value;
                        }
                        else if (mod_op == "set")
                        {
                            backend_config[keyword][sub_keyword] = value;
                        }
                        else if (mod_op == "reset")
                        {
                            backend_config[keyword][sub_keyword] = default_backend_config[keyword][sub_keyword];
                        }
                        else
                        {
                            std::cerr << "Unknown single qubit noise type: " << mod_op << std::endl;
                        }
                    }
                }
            }
            else
            {
                // case 3: single-qubit property
                for (auto qubit : qubit_list)
                {
                    std::string sub_keyword = qindex(qubit); // Generate index for single qubit

                    // Perform the operation
                    if (mod_op == "scale")
                    {
                        backend_config[keyword][sub_keyword] = static_cast<ValType>(backend_config[keyword][sub_keyword]) * value;
                    }
                    else if (mod_op == "set")
                    {
                        backend_config[keyword][sub_keyword] = value;
                    }
                    else if (mod_op == "reset")
                    {
                        backend_config[keyword][sub_keyword] = default_backend_config[keyword][sub_keyword];
                    }
                    else
                    {
                        std::cerr << "Unknown single qubit noise type: " << mod_op << std::endl;
                    }
                }
            }

            // print();
            // std::cout << "\n\n";
        }

        void loadLayoutFile(const std::string &filename, bool required = true)
        {
            std::ifstream i(filename);
            if (!i.is_open())
            {
                return;
                // throw std::runtime_error("Could not open " + filename);
            }
            layout = nlohmann::json::parse(i);
            i.close();
        }

        void loadLayoutString(std::string s)
        {
            size_t pos_start = 0;
            size_t pos_end;
            std::string token;
            while ((pos_end = s.find(",", pos_start)) != std::string::npos)
            {
                token = s.substr(pos_start, pos_end - pos_start);
                size_t src_end = token.find("=");
                if (src_end == std::string::npos)
                {
                    throw std::invalid_argument("Ill-formatted layout string\n");
                }
                std::string log_qb = token.substr(0, src_end);
                std::string phys_qb = token.substr(src_end + 1);
                layout[log_qb] = std::stoll(phys_qb);
                pos_start = pos_end + 1;
            }
            token = s.substr(pos_start);
            if (token.length())
            {
                size_t src_end = token.find("=");
                std::string log_qb = token.substr(0, src_end);
                std::string phys_qb = token.substr(src_end + 1);
                layout[log_qb] = std::stoll(phys_qb);
            }
        }

        void print(int precision = 3)
        {
            safe_print("Noise Profile:\n");

            // Print basic information
            safe_print("Device Name: ", backend_config["name"], "\n");
            safe_print("Version: ", backend_config["version"], "\n");
            safe_print("Number of Qubits: ", backend_config["num_qubits"], "\n");

            // Print basis gates
            safe_print("Basis Gates: ");
            for (const auto &gate : backend_config["basis_gates"])
            {
                safe_print(gate, " ");
            }
            safe_print("\n");

            // Print T1 values in scientific notation
            safe_print("T1 Times (in seconds):\n");
            for (auto &[qubit, value] : backend_config["T1"].items())
            {
                safe_print("  Qubit ", qubit, ": ", formatWithPrecision(value, precision), "\n");
            }

            // Print T2 values in scientific notation
            safe_print("T2 Times (in seconds):\n");
            for (auto &[qubit, value] : backend_config["T2"].items())
            {
                safe_print("  Qubit ", qubit, ": ", formatWithPrecision(value, precision), "\n");
            }

            // Print frequencies in scientific notation
            safe_print("Frequencies (in Hz):\n");
            for (auto &[qubit, value] : backend_config["freq"].items())
            {
                safe_print("  Qubit ", qubit, ": ", formatWithPrecision(value, precision), "\n");
            }

            // Print readout lengths in scientific notation
            safe_print("Readout Lengths (in seconds):\n");
            for (auto &[qubit, value] : backend_config["readout_length"].items())
            {
                safe_print("  Qubit ", qubit, ": ", formatWithPrecision(value, precision), "\n");
            }

            // Print measurement error probabilities in scientific notation
            safe_print("Measurement Error Probabilities (Meas0_Prep1):\n");
            for (auto &[qubit, value] : backend_config["prob_meas0_prep1"].items())
            {
                safe_print("  Qubit ", qubit, ": ", formatWithPrecision(value, precision), "\n");
            }

            safe_print("Measurement Error Probabilities (Meas1_Prep0):\n");
            for (auto &[qubit, value] : backend_config["prob_meas1_prep0"].items())
            {
                safe_print("  Qubit ", qubit, ": ", formatWithPrecision(value, precision), "\n");
            }

            // Print gate lengths in scientific notation
            safe_print("Gate Lengths (in seconds):\n");
            for (auto &[gate, value] : backend_config["gate_lens"].items())
            {
                safe_print("  ", gate, ": ", formatWithPrecision(value, precision), "\n");
            }

            // Print gate errors in scientific notation
            safe_print("Gate Errors:\n");
            for (auto &[gate, value] : backend_config["gate_errs"].items())
            {
                safe_print("  ", gate, ": ", formatWithPrecision(value, precision), "\n");
            }

            // Print coupling map
            safe_print("CX Coupling Map:\n");
            for (const auto &coupling : backend_config["cx_coupling"])
            {
                safe_print("  ", coupling, "\n");
            }

            safe_print("========================\n");
        }

        // void print(int precision = 3)
        // {
        //     std::cout << "Noise Profile:\n";

        //     // Print basic information
        //     std::cout << "Device Name: " << backend_config["name"] << "\n";
        //     std::cout << "Version: " << backend_config["version"] << "\n";
        //     std::cout << "Number of Qubits: " << backend_config["num_qubits"] << "\n";

        //     // Print basis gates
        //     std::cout << "Basis Gates: ";
        //     for (const auto &gate : backend_config["basis_gates"])
        //     {
        //         std::cout << gate << " ";
        //     }
        //     std::cout << "\n";

        //     // Print T1 values in scientific notation
        //     std::cout << "T1 Times (in seconds):\n";
        //     for (auto &[qubit, value] : backend_config["T1"].items())
        //     {
        //         std::cout << "  Qubit " << qubit << ": " << formatWithPrecision(value, precision) << "\n";
        //     }

        //     // Print T2 values in scientific notation
        //     std::cout << "T2 Times (in seconds):\n";
        //     for (auto &[qubit, value] : backend_config["T2"].items())
        //     {
        //         std::cout << "  Qubit " << qubit << ": " << formatWithPrecision(value, precision) << "\n";
        //     }

        //     // Print frequencies in scientific notation
        //     std::cout << "Frequencies (in Hz):\n";
        //     for (auto &[qubit, value] : backend_config["freq"].items())
        //     {
        //         std::cout << "  Qubit " << qubit << ": " << formatWithPrecision(value, precision) << "\n";
        //     }

        //     // Print readout lengths in scientific notation
        //     std::cout << "Readout Lengths (in seconds):\n";
        //     for (auto &[qubit, value] : backend_config["readout_length"].items())
        //     {
        //         std::cout << "  Qubit " << qubit << ": " << formatWithPrecision(value, precision) << "\n";
        //     }

        //     // Print measurement error probabilities in scientific notation
        //     std::cout << "Measurement Error Probabilities (Meas0_Prep1):\n";
        //     for (auto &[qubit, value] : backend_config["prob_meas0_prep1"].items())
        //     {
        //         std::cout << "  Qubit " << qubit << ": " << formatWithPrecision(value, precision) << "\n";
        //     }

        //     std::cout << "Measurement Error Probabilities (Meas1_Prep0):\n";
        //     for (auto &[qubit, value] : backend_config["prob_meas1_prep0"].items())
        //     {
        //         std::cout << "  Qubit " << qubit << ": " << formatWithPrecision(value, precision) << "\n";
        //     }

        //     // Print gate lengths in scientific notation
        //     std::cout << "Gate Lengths (in seconds):\n";
        //     for (auto &[gate, value] : backend_config["gate_lens"].items())
        //     {
        //         std::cout << "  " << gate << ": " << formatWithPrecision(value, precision) << "\n";
        //     }

        //     // Print gate errors in scientific notation
        //     std::cout << "Gate Errors:\n";
        //     for (auto &[gate, value] : backend_config["gate_errs"].items())
        //     {
        //         std::cout << "  " << gate << ": " << formatWithPrecision(value, precision) << "\n";
        //     }

        //     // Print coupling map
        //     std::cout << "CX Coupling Map:\n";
        //     for (const auto &coupling : backend_config["cx_coupling"])
        //     {
        //         std::cout << "  " << coupling << "\n";
        //     }

        //     std::cout << "========================\n";
        // }

    private:
        std::unordered_map<std::string, double> defaultNoiseParameters;

        nlohmann::json backend_config = {};
        nlohmann::json layout = {};

        nlohmann::json default_backend_config = {};
        nlohmann::json default_layout = {};

        std::string qindex(IdxType q)
        {
            return std::to_string(static_cast<IdxType>(
                layout[std::to_string(q)]));
        }

        std::string gname(const std::string &gate_name, IdxType q1, IdxType q2 = -1)
        {
            std::string name = gate_name + qindex(q1);

            if (q2 >= 0)
            {
                name += "_" + qindex(q2);
            }
            return name;
        }

        // Helper function to format a number with fixed precision
        std::string formatWithPrecision(double value, int precision)
        {
            if (std::abs(value) < 1e-10)
            {
                return "0";
            }

            std::ostringstream out;
            out << std::scientific << std::setprecision(precision) << value;
            return out.str();
        }
    };
}
#endif // NOISE_MODEL_HPP