#pragma once

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <iomanip>  // For std::setw
#include <getopt.h> // For getopt_long
#include <memory>

// Class to represent a single configuration option
class ConfigOption
{
public:
    std::string long_name;     // Long flag (e.g., --shots)
    char short_name;           // Short flag (e.g., -q)
    std::string description;   // Description of the option
    std::string value;         // Current value of the option (for argument flags)
    std::string default_value; // Default value of the option (for argument flags)
    std::string desc_value;    // Description of the value (for help messages)
    bool requires_argument;    // Whether the option requires an argument
    bool flag_set;             // For boolean flags: true if flag is provided

    // Default constructor
    ConfigOption() : long_name(""), short_name(0), description(""), value(""), default_value(""), desc_value(""), requires_argument(false), flag_set(false) {}

    // Constructor for initializing a configuration option
    ConfigOption(const std::string &long_name, char short_name, const std::string &default_value, const std::string &description, const std::string &desc_value, bool requires_argument)
        : long_name(long_name), short_name(short_name), description(description), value(default_value), default_value(default_value), desc_value(desc_value), requires_argument(requires_argument), flag_set(false) {}

    // Set value from string (only for argument flags)
    void set_value_from_string(const std::string &input)
    {
        value = input;
    }

    // Set flag to true (only for boolean flags)
    void set_flag()
    {
        flag_set = true;
    }
};

// Class to hold and manage all configuration options
class ConfigParser
{
public:
    // Use shared pointers to point to the same ConfigOption object
    std::map<std::string, std::shared_ptr<ConfigOption>> long_options_map;
    std::map<char, std::shared_ptr<ConfigOption>> short_options_map;

    // Define the legacy options that need to be replaced
    const std::map<std::string, std::string> legacy_options = {
        {"-qs", "--qasm_string"},
        {"-js", "--json_string"},
        {"-metrics", "--metrics"},
        {"-initial", "--init_file"},
        {"-initial-format", "--init_format"},
        {"-dump", "--dump_file"},
        {"-layout", "--layout"},
        {"-layout-string", "--layout_str"},
        {"-device", "--device"},
        {"-backend", "--backend"},
        {"-backend_list", "--backend_list"},
        {"-shots", "--shots"},
        {"-sim", "--sim"},
        {"-basis", "--basis"},
        {"-fidelity", "--fidelity"}};

    // Constructor initializes default options
    ConfigParser()
    {
        // Initialize configuration options (LONG_FLAG, SHORT_FLAG, DEFAULT_VALUE, DESCRIPTION, DESC_VALUE, REQUIRES_ARGUMENT)

        // Program execution options
        add_option("qasm_file", 'q', "", "Execute simulation with the given QASM file", "FILE_PATH", true);
        add_option("qasm_string", 0, "", "Execute simulation with the provided QASM string", "STR", true); // Long flag only
        add_option("json_file", 'j', "", "Execute simulation with the given JSON file (Qiskit Qobj)", "FILE_PATH", true);
        add_option("json_string", 0, "", "Execute simulation with the provided JSON string (Qiskit Qobj)", "STR", true); // Long flag only
        add_option("test", 't', "", "Run testing benchmarks for the specified index", "INT", true);
        add_option("all_tests", 'a', "", "Run all available testing benchmarks", "", false); // Boolean flag, no argument

        // Circuit execution options
        add_option("shots", 's', "1024", "Specify the number of shots", "SHOTS", true);
        add_option("backend", 'b', "cpu", "Specify the simulation backend", "BACKEND", true);
        add_option("sim", 0, "sv", "Specify the simulation method", "METHOD", true);             // Long flag only
        add_option("basis", 0, "", "Run the test benchmark using basis gates", "", false);       // Boolean flag, no argument
        add_option("disable_fusion", 0, "", "Disable gate fusion ", "", false);                  // Boolean flag, no argument
        add_option("random_seed", 0, "", "Set the random seed for the simulation", "INT", true); // Long flag only

        // Noise model options
        add_option("device", 0, "", "Specify the device noise profile", "FILE_PATH", true);                       // Long flag only
        add_option("layout", 0, "", "Path to JSON mapping logical qubits to physical qubits", "FILE_PATH", true); // Long flag only
        add_option("layout_str", 0, "", "String format mapping logical qubits to physical qubits", "STR", true);  // Long flag only

        // Initial and resulting state file options
        add_option("init_file", 0, "", "Path to the initial statevector/density matrix file", "FILE_PATH", true);       // Long flag only
        add_option("init_format", 0, "", "Specify the format of the initial state", "FILE_PATH", true);                 // Long flag only
        add_option("dump_file", 0, "", "Path to dump the binary statevector/density matrix result", "FILE_PATH", true); // Long flag only

        // Helper options
        add_option("backend_list", 0, "", "Print the list of available simulation backends", "", false);    // Boolean flag, no argument
        add_option("metrics", 0, "", "Print the metrics of the executed circuit", "", false);               // Boolean flag, no argument
        add_option("verbose", 'v', "", "Enable verbose simulation trace", "", false);                       // Boolean flag, with short flag
        add_option("fidelity", 'f', "", "Run both DM-Sim and SV-Sim and report state fidelity", "", false); // Boolean flag, no argument
        add_option("help", 'h', "", "Print this help message", "", false);                                  // Boolean flag, with short flag

        // Hardware options (prefixed with hw_)
        add_option("hw_tensorcore", 0, "", "Enable the use of Tensor Cores", "", false);     // Boolean flag, no argument
        add_option("hw_threads", 0, "-1", "Specify the number of OMP threads", "INT", true); // Long flag only
        add_option("hw_avx512", 0, "", "Enable the use of AVX512", "", false);               // Boolean flag, no argument
        add_option("hw_matrixcore", 0, "", "Enable the use of MatrixCore", "", false);       // Boolean flag, no argument
    }

    // Method to add an option with both long and short flags
    void add_option(const std::string &long_name, char short_name, const std::string &default_value, const std::string &description, const std::string &desc_value, bool requires_argument)
    {
        // Create a shared pointer to the ConfigOption object
        std::shared_ptr<ConfigOption> option = std::make_shared<ConfigOption>(long_name, short_name, default_value, description, desc_value, requires_argument);

        // Add it to the long options map if long_name is provided
        if (!long_name.empty())
        {
            long_options_map[long_name] = option;
        }

        // Add the same object to the short options map if a short_name is provided
        if (short_name != 0)
        {
            short_options_map[short_name] = option;
        }
    }

    // Get the value of a configuration option (only for options that require arguments)
    std::string get_value(const std::string &name) const
    {
        // Check if it's a long flag
        if (long_options_map.find(name) != long_options_map.end())
        {
            return long_options_map.at(name)->value;
        }

        // Raise an error if the option is not found
        std::cerr << "Error: Option '" << name << "' not found.\n";
        return "";
    }

    std::string get_value(char short_flag) const
    {
        // Check if it's a short flag
        if (short_options_map.find(short_flag) != short_options_map.end())
        {
            return short_options_map.at(short_flag)->value;
        }

        // Raise an error if the option is not found
        std::cerr << "Error: Option '" << short_flag << "' not found.\n";
        return "";
    }

    // Check if a flag is set (only for boolean flags)
    bool is_flag_set(const std::string &name) const
    {
        // Check if it's a long flag
        if (long_options_map.find(name) != long_options_map.end())
        {
            return long_options_map.at(name)->flag_set;
        }

        // Raise an error if the option is not found
        std::cerr << "Error: Option '" << name << "' not found.\n";
        return false;
    }

    bool is_flag_set(char short_flag) const
    {
        // Check if it's a short flag
        if (short_options_map.find(short_flag) != short_options_map.end())
        {
            return short_options_map.at(short_flag)->flag_set;
        }

        // Raise an error if the option is not found
        std::cerr << "Error: Option '" << short_flag << "' not found.\n";
        return false;
    }

    // Function to set a flag or assign a value
    void set_flag_or_value(const std::string &name, const std::string &value = "")
    {
        // First check if it's a long flag
        if (long_options_map.find(name) != long_options_map.end())
        {
            auto opt = long_options_map[name];
            if (opt->requires_argument)
            {
                opt->set_value_from_string(value); // Set the value
            }
            else
            {
                opt->set_flag(); // Set the flag for boolean options
            }
        }
        // Check for short flag
        else if (name.size() == 1)
        {
            char short_name = name[0];
            if (short_options_map.find(short_name) != short_options_map.end())
            {
                auto opt = short_options_map[short_name];
                if (opt->requires_argument)
                {
                    opt->set_value_from_string(value); // Set the value
                }
                else
                {
                    opt->set_flag(); // Set the flag for boolean options
                }
            }
        }
    }

    void preprocess_single_dash_options(int argc, char *argv[])
    {
        // Iterate over the arguments, starting from argv[1] (skipping the program name)
        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];

            // Check if the argument is a legacy single-dash option
            if (legacy_options.find(arg) != legacy_options.end())
            {
                // Print a warning to the user
                std::cout << "Warning: '" << arg << "' is using a single dash. Please migrate to '" << legacy_options.at(arg) << "'.\n";

                // Replace the single-dash option with the double-dash version
                argv[i] = const_cast<char *>(legacy_options.at(arg).c_str());

                std::cout << "Replaced with: " << argv[i] << std::endl;
            }
        }
    }

    void parse_arguments(int argc, char *argv[])
    {

        // Preprocess single-dash options that were used in previous versions
        preprocess_single_dash_options(argc, argv);

        // Build the long option structure for getopt_long
        std::vector<struct option> long_options;
        std::string short_options;

        // Convert the config options into long options for getopt_long
        for (const auto &entry : long_options_map)
        {
            const std::shared_ptr<ConfigOption> &opt = entry.second;

            // Add long option to the long_options array
            struct option long_option;
            long_option.name = opt->long_name.c_str();                                      // Long option name
            long_option.has_arg = opt->requires_argument ? required_argument : no_argument; // Does this option require an argument?
            long_option.flag = nullptr;                                                     // This is null since we are not using flags directly
            long_option.val = opt->short_name != 0 ? opt->short_name : 0;                   // Short flag value, or 0 if none
            long_options.push_back(long_option);

            // Add short option if applicable
            if (opt->short_name != 0)
            {
                short_options += opt->short_name;
                if (opt->requires_argument)
                {
                    short_options += ":"; // `:` means this short option requires an argument
                }
            }
        }

        // Add a terminating zeroed-out option to indicate the end of options
        long_options.push_back({0, 0, 0, 0});

        // Parsing arguments using getopt_long
        int option_index = 0;
        int c;
        while ((c = getopt_long(argc, argv, short_options.c_str(), long_options.data(), &option_index)) != -1)
        {
            switch (c)
            {
            case 0: // This case is triggered for long flags without a short flag
            {
                const char *option_name = long_options[option_index].name;
                if (long_options_map.count(option_name))
                {
                    std::shared_ptr<ConfigOption> &opt = long_options_map[option_name];
                    if (opt->requires_argument)
                    {
                        opt->set_value_from_string(optarg); // Set the argument value
                    }
                    else
                    {
                        opt->set_flag(); // Set the boolean flag to true
                    }
                }
                break;
            }

            case '?': // Unrecognized option
                std::cerr << "Unknown option: " << argv[optind - 1] << std::endl;
                exit(EXIT_FAILURE);

            default: // Short flag or long flag with short equivalent
            {
                char short_flag = static_cast<char>(c);
                if (short_options_map.count(short_flag))
                {
                    std::shared_ptr<ConfigOption> &opt = short_options_map[short_flag];
                    if (opt->requires_argument)
                    {
                        opt->set_value_from_string(optarg); // Set the argument value
                    }
                    else
                    {
                        opt->set_flag(); // Set the boolean flag to true
                    }
                }
                break;
            }
            }
        }
    }

    void print_help() const
    {
        std::cout << "Usage: program [OPTIONS]\n\n";
        std::cout << "Available options:\n";

        // Define column widths for consistent alignment
        const int short_flag_width = 6;
        const int long_flag_width = 17;
        const int argument_width = 13;
        const int description_width = 80 - short_flag_width - long_flag_width - argument_width;

        for (const auto &entry : long_options_map)
        {
            auto opt = entry.second;

            // Print short flag if it exists
            if (opt->short_name != 0)
            {
                std::cout << "  -" << opt->short_name << ", ";
            }
            else
            {
                std::cout << "      "; // Maintain alignment if no short flag
            }

            // Print long flag
            std::cout << std::setw(long_flag_width) << std::left << "--" + opt->long_name;

            // Print required argument if applicable
            if (opt->requires_argument)
            {
                std::cout << std::setw(argument_width) << std::left << "<" + opt->desc_value + ">";
            }
            else
            {
                std::cout << std::setw(argument_width) << " "; // Maintain alignment if no argument
            }

            std::cout << std::setw(description_width) << std::left << opt->description << "\n";
        }
    }

    // Print the current configuration settings for debugging
    void print_configurations() const
    {
        safe_print("Current Configuration Settings:\n");

        // Loop through the long options map and print their values
        for (const auto &entry : long_options_map)
        {
            auto opt = entry.second;

            // Debugging: Print current and default values
            // std::cout << opt->long_name << ", Current Value: '" << opt->value << "', Default Value: '" << opt->default_value << "'\n";

            // Check if the current value is different from the default value
            bool is_non_default = (opt->value != opt->default_value);

            // Print the flag name and its current value
            safe_print("--%s = %s", opt->long_name.c_str(), opt->value.c_str());

            // Indicate non-default values with an asterisk
            if (is_non_default && !opt->value.empty())
            {
                safe_print(" *");
            }

            // For boolean flags, check if they were set
            if (!opt->requires_argument)
            {
                safe_print(" (Flag: %s)", opt->flag_set ? "ON" : "OFF");
            }

            safe_print("\n");
        }
    }
};
