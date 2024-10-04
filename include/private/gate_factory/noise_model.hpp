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

using json = nlohmann::json;
namespace NWQSim
{
    // Constants for qubit parameters
    constexpr const char *T1_KEY = "T1";
    constexpr const char *T2_KEY = "T2";
    constexpr const char *FREQ_KEY = "freq";
    constexpr const char *READOUT_LENGTH_KEY = "readout_length";
    constexpr const char *PROB_MEAS0_PREP1_KEY = "prob_meas0_prep1";
    constexpr const char *PROB_MEAS1_PREP0_KEY = "prob_meas1_prep0";

    // Constants for gate parameters
    constexpr const char *GATE_LENS_KEY = "gate_lens";
    constexpr const char *GATE_ERRS_KEY = "gate_errs";

    class NoiseModel
    {

    public:
        NoiseModel(const std::string &device_filename, const std::string &mapping_filename, const std::string &mapping_str)
        {
            load_from_file(device_filename);
            save_defaults();

            if (!mapping_str.empty())
            {
                loadMappingString(mapping_str);
            }
            else if (!mapping_filename.empty())
            {
                loadMappingFile(mapping_filename);
            }
            else
            {
                for (IdxType n = 0; n < num_qubits; n++)
                {
                    qubit_mapping[n] = n;
                }
            }
        }

        // Method to load the JSON file and set the current values
        void load_from_file(const std::string &file_path)
        {
            std::ifstream input_file(file_path);
            if (!input_file.is_open())
            {
                throw std::runtime_error("Could not open file");
            }

            json j;
            input_file >> j;

            name = j["name"];
            version = j["version"];
            num_qubits = j["num_qubits"];
            basis_gates = j["basis_gates"].get<std::vector<std::string>>();

            // Populate qubit-specific parameters
            qubit_parameters[T1_KEY] = j[T1_KEY].get<std::map<std::string, double>>();
            qubit_parameters[T2_KEY] = j[T2_KEY].get<std::map<std::string, double>>();
            qubit_parameters[FREQ_KEY] = j[FREQ_KEY].get<std::map<std::string, double>>();
            qubit_parameters[READOUT_LENGTH_KEY] = j[READOUT_LENGTH_KEY].get<std::map<std::string, double>>();
            qubit_parameters[PROB_MEAS0_PREP1_KEY] = j[PROB_MEAS0_PREP1_KEY].get<std::map<std::string, double>>();
            qubit_parameters[PROB_MEAS1_PREP0_KEY] = j[PROB_MEAS1_PREP0_KEY].get<std::map<std::string, double>>();

            // Populate gate-specific parameters
            gate_parameters[GATE_LENS_KEY] = j[GATE_LENS_KEY].get<std::map<std::string, double>>();
            gate_parameters[GATE_ERRS_KEY] = j[GATE_ERRS_KEY].get<std::map<std::string, double>>();

            cx_coupling = j["cx_coupling"].get<std::vector<std::string>>();
        }

        // Method to load qubit mapping from a file
        void loadMappingFile(const std::string &filename, bool required = true)
        {
            std::ifstream input_file(filename);
            if (!input_file.is_open())
            {
                if (required)
                {
                    throw std::runtime_error("Could not open mapping file: " + filename);
                }
                return; // If the file is not required, just return without loading
            }

            json j;
            input_file >> j;
            input_file.close();

            // Populate the qubit_mapping object from the JSON file
            for (const auto &item : j.items())
            {
                IdxType logical_qb = std::stoi(item.key());
                IdxType physical_qb = item.value().get<IdxType>();
                qubit_mapping[logical_qb] = physical_qb;
            }
        }

        // Method to load qubit mapping from a string
        void loadMappingString(const std::string &s)
        {
            size_t pos_start = 0;
            size_t pos_end;
            std::string token;

            // Parse the string "logical=physical,logical=physical"
            while ((pos_end = s.find(",", pos_start)) != std::string::npos)
            {
                token = s.substr(pos_start, pos_end - pos_start);
                size_t src_end = token.find("=");
                if (src_end == std::string::npos)
                {
                    throw std::invalid_argument("Ill-formatted mapping string");
                }
                IdxType log_qb = std::stoi(token.substr(0, src_end));
                IdxType phys_qb = std::stoi(token.substr(src_end + 1));
                qubit_mapping[log_qb] = phys_qb; // Store the mapping
                pos_start = pos_end + 1;
            }

            // Handle the last pair
            token = s.substr(pos_start);
            if (!token.empty())
            {
                size_t src_end = token.find("=");
                IdxType log_qb = std::stoi(token.substr(0, src_end));
                IdxType phys_qb = std::stoi(token.substr(src_end + 1));
                qubit_mapping[log_qb] = phys_qb;
            }
        }

        // Method to save the current values as defaults
        void save_defaults()
        {
            default_name = name;
            default_version = version;
            default_num_qubits = num_qubits;
            default_basis_gates = basis_gates;

            default_qubit_parameters = qubit_parameters;
            default_gate_parameters = gate_parameters;
            default_cx_coupling = cx_coupling;
        }

        // Method to reset the current values to defaults
        void reset_to_defaults()
        {
            name = default_name;
            version = default_version;
            num_qubits = default_num_qubits;
            basis_gates = default_basis_gates;

            qubit_parameters = default_qubit_parameters;
            gate_parameters = default_gate_parameters;
            cx_coupling = default_cx_coupling;
        }

        IdxType get_num_qubits()
        {
            return num_qubits;
        }

        // Method to safely get a qubit parameter
        double get_qubit_parameter(const std::string &param, IdxType qubit) const
        {
            std::string physical_qubit = get_physical_qubit(qubit);
            if (qubit_parameters.find(param) != qubit_parameters.end())
            {
                const auto &qubit_map = qubit_parameters.at(param);
                if (qubit_map.find(physical_qubit) != qubit_map.end())
                {
                    return qubit_map.at(physical_qubit);
                }
            }
            std::cerr << "Warning: Qubit parameter not found: " << param << ", Qubit " << qubit << "\n";
            return 0.0; // Return a default value or handle it as necessary
        }

        // Method to safely get a gate parameter
        double get_gate_parameter(const std::string &param, const std::string &gate, IdxType q1, IdxType q2 = -1) const
        {
            std::string gate_name = get_gate_name(gate, q1, q2);

            if (gate_parameters.find(param) != gate_parameters.end())
            {
                const auto &gate_map = gate_parameters.at(param);
                if (gate_map.find(gate_name) != gate_map.end())
                {
                    return gate_map.at(gate_name);
                }
            }
            std::cerr << "Warning: Gate parameter not found: " << param << ", Gate " << gate_name << "\n";
            return 0.0; // Return a default value or handle it as necessary
        }

        void modify_noise(std::string mod_op, std::string mod_noise, ValType value, std::vector<IdxType> qubit_list)
        {
            // Convert mod_op and mod_noise to lowercase
            std::transform(mod_op.begin(), mod_op.end(), mod_op.begin(), ::tolower);
            std::transform(mod_noise.begin(), mod_noise.end(), mod_noise.begin(), ::tolower);

            std::string keyword;     // This will hold the relevant keyword for qubit or gate parameters
            bool is_gate = false;    // Whether the modification is for a gate or not
            bool is_2q_gate = false; // Whether the gate is a 2-qubit gate or not
            std::string gate_name;   // To store the gate name for gate-related noise

            // Determine the correct keyword based on mod_noise
            if (mod_noise == "t1")
            {
                keyword = T1_KEY;
            }
            else if (mod_noise == "t2")
            {
                keyword = T2_KEY;
            }
            else if (mod_noise == "readout_len")
            {
                keyword = READOUT_LENGTH_KEY;
            }
            else if (mod_noise == "readout_m0p1")
            {
                keyword = PROB_MEAS0_PREP1_KEY;
            }
            else if (mod_noise == "readout_m1p0")
            {
                keyword = PROB_MEAS1_PREP0_KEY;
            }
            else if (mod_noise.find("_len") != std::string::npos)
            {
                // Extract gate name before "_len" and set keyword to "gate_lens"
                gate_name = mod_noise.substr(0, mod_noise.find("_len"));
                keyword = GATE_LENS_KEY;
                is_gate = true;
            }
            else if (mod_noise.find("_err") != std::string::npos)
            {
                // Extract gate name before "_err" and set keyword to "gate_errs"
                gate_name = mod_noise.substr(0, mod_noise.find("_err"));
                keyword = GATE_ERRS_KEY;
                is_gate = true;
            }
            else
            {
                std::cerr << "Unknown noise type: " << mod_noise << std::endl;
                return;
            }

            if (is_gate) // Gate noise (either 1-qubit or 2-qubit)
            {
                is_2q_gate = (gate_name == "cx" || gate_name == "ecr"); // Check if the gate is a 2-qubit gate
                if (is_2q_gate)
                    assert(qubit_list.size() == 2);

                if (is_2q_gate) // case 1: two-qubit gate
                {
                    IdxType q1 = qubit_list[1]; // target
                    IdxType q2 = qubit_list[0]; // ctrl
                    std::string gate_key = get_gate_name(gate_name, q1, q2);

                    // Perform the operation
                    if (mod_op == "scale")
                    {
                        gate_parameters[keyword][gate_key] *= value;
                    }
                    else if (mod_op == "set")
                    {
                        gate_parameters[keyword][gate_key] = value;
                    }
                    else if (mod_op == "reset")
                    {
                        gate_parameters[keyword][gate_key] = default_gate_parameters[keyword][gate_key];
                    }
                    else
                    {
                        std::cerr << "Unknown two-qubit gate modification operation: " << mod_op << std::endl;
                    }
                }
                else // case 2: single-qubit gate
                {
                    for (auto qubit : qubit_list)
                    {
                        std::string gate_key = get_gate_name(gate_name, qubit, -1);

                        // Perform the operation
                        if (mod_op == "scale")
                        {
                            gate_parameters[keyword][gate_key] *= value;
                        }
                        else if (mod_op == "set")
                        {
                            gate_parameters[keyword][gate_key] = value;
                        }
                        else if (mod_op == "reset")
                        {
                            gate_parameters[keyword][gate_key] = default_gate_parameters[keyword][gate_key];
                        }
                        else
                        {
                            std::cerr << "Unknown single-qubit gate modification operation: " << mod_op << std::endl;
                        }
                    }
                }
            }
            else // qubit properties (T1, T2, readout, etc.)
            {
                for (auto qubit : qubit_list)
                {
                    std::string physical_qubit = get_physical_qubit(qubit);
                    // Perform the operation
                    if (mod_op == "scale")
                    {
                        qubit_parameters[keyword][physical_qubit] *= value;
                    }
                    else if (mod_op == "set")
                    {
                        qubit_parameters[keyword][physical_qubit] = value;
                    }
                    else if (mod_op == "reset")
                    {
                        qubit_parameters[keyword][physical_qubit] = default_qubit_parameters[keyword][physical_qubit];
                    }
                    else
                    {
                        std::cerr << "Unknown qubit property modification operation: " << mod_op << std::endl;
                    }
                }
            }
        }

        // Method to print the contents of the noise profile
        void print(int precision = 3) const
        {
            std::cout << "Noise Profile:\n";

            // Print basic information
            std::cout << "Device Name: " << name << "\n";
            std::cout << "Version: " << version << "\n";
            std::cout << "Number of Qubits: " << num_qubits << "\n";

            // Print basis gates
            std::cout << "Basis Gates: ";
            for (const auto &gate : basis_gates)
            {
                std::cout << gate << " ";
            }
            std::cout << "\n";

            // Print T1 values using formatWithPrecision
            std::cout << "T1 Times (in seconds):\n";
            if (qubit_parameters.find("T1") != qubit_parameters.end())
            {
                for (const auto &[qubit, value] : qubit_parameters.at("T1"))
                {
                    std::cout << "  Qubit " << qubit << ": " << formatWithPrecision(value, precision) << "\n";
                }
            }

            // Print T2 values using formatWithPrecision
            std::cout << "T2 Times (in seconds):\n";
            if (qubit_parameters.find("T2") != qubit_parameters.end())
            {
                for (const auto &[qubit, value] : qubit_parameters.at("T2"))
                {
                    std::cout << "  Qubit " << qubit << ": " << formatWithPrecision(value, precision) << "\n";
                }
            }

            // Print frequencies using formatWithPrecision
            std::cout << "Frequencies (in Hz):\n";
            if (qubit_parameters.find("freq") != qubit_parameters.end())
            {
                for (const auto &[qubit, value] : qubit_parameters.at("freq"))
                {
                    std::cout << "  Qubit " << qubit << ": " << formatWithPrecision(value, precision) << "\n";
                }
            }

            // Print readout lengths using formatWithPrecision
            std::cout << "Readout Lengths (in seconds):\n";
            if (qubit_parameters.find("readout_length") != qubit_parameters.end())
            {
                for (const auto &[qubit, value] : qubit_parameters.at("readout_length"))
                {
                    std::cout << "  Qubit " << qubit << ": " << formatWithPrecision(value, precision) << "\n";
                }
            }

            // Print measurement error probabilities for Meas0_Prep1
            std::cout << "Measurement Error Probabilities (Meas0_Prep1):\n";
            if (qubit_parameters.find("prob_meas0_prep1") != qubit_parameters.end())
            {
                for (const auto &[qubit, value] : qubit_parameters.at("prob_meas0_prep1"))
                {
                    std::cout << "  Qubit " << qubit << ": " << formatWithPrecision(value, precision) << "\n";
                }
            }

            // Print measurement error probabilities for Meas1_Prep0
            std::cout << "Measurement Error Probabilities (Meas1_Prep0):\n";
            if (qubit_parameters.find("prob_meas1_prep0") != qubit_parameters.end())
            {
                for (const auto &[qubit, value] : qubit_parameters.at("prob_meas1_prep0"))
                {
                    std::cout << "  Qubit " << qubit << ": " << formatWithPrecision(value, precision) << "\n";
                }
            }

            // Print gate lengths using formatWithPrecision
            std::cout << "Gate Lengths (in seconds):\n";
            if (gate_parameters.find("gate_lens") != gate_parameters.end())
            {
                for (const auto &[gate, value] : gate_parameters.at("gate_lens"))
                {
                    std::cout << "  " << gate << ": " << formatWithPrecision(value, precision) << "\n";
                }
            }

            // Print gate errors using formatWithPrecision
            std::cout << "Gate Errors:\n";
            if (gate_parameters.find("gate_errs") != gate_parameters.end())
            {
                for (const auto &[gate, value] : gate_parameters.at("gate_errs"))
                {
                    std::cout << "  " << gate << ": " << formatWithPrecision(value, precision) << "\n";
                }
            }

            // Print coupling map
            std::cout << "CX Coupling Map:\n";
            for (const auto &coupling : cx_coupling)
            {
                std::cout << "  " << coupling << "\n";
            }

            std::cout << "========================\n";
        }

        // Method to print the current mapping
        void print_mapping() const
        {
            std::cout << "Qubit Mapping (Logical -> Physical):\n";
            for (const auto &[logical, physical] : qubit_mapping)
            {
                std::cout << "Logical Qubit: " << logical << " -> Physical Qubit: " << physical << "\n";
            }
        }

    private:
        // Current values
        std::string name;
        std::string version;
        IdxType num_qubits;
        std::vector<std::string> basis_gates;

        // Consolidated maps
        std::map<std::string, std::map<std::string, double>> qubit_parameters;
        std::map<std::string, std::map<std::string, double>> gate_parameters;
        std::vector<std::string> cx_coupling;

        // Backup for defaults
        std::string default_name;
        std::string default_version;
        IdxType default_num_qubits;
        std::vector<std::string> default_basis_gates;

        std::map<std::string, std::map<std::string, double>> default_qubit_parameters;
        std::map<std::string, std::map<std::string, double>> default_gate_parameters;
        std::vector<std::string> default_cx_coupling;

        // Mapping object: Logical qubit to physical qubit mapping

        std::map<IdxType, IdxType> qubit_mapping;

        // Method to get the physical qubit for a given logical qubit
        std::string get_physical_qubit(IdxType logical_qubit) const
        {
            if (qubit_mapping.find(logical_qubit) != qubit_mapping.end())
            {
                return std::to_string(qubit_mapping.at(logical_qubit));
            }
            else
            {
                std::cerr << "Warning: Logical qubit " << logical_qubit << " not found in mapping.\n";
                return ""; // Return an invalid value or handle it appropriately
            }
        }

        std::string get_gate_name(const std::string &gate_name, IdxType q1, IdxType q2) const
        {
            std::string name = gate_name;

            if (q2 >= 0)
            {
                name += get_physical_qubit(q2) + "_" + get_physical_qubit(q1);
            }
            else
            {
                name += get_physical_qubit(q1);
            }
            return name;
        }

        // Helper function to format a number with fixed precision
        std::string formatWithPrecision(double value, int precision) const
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