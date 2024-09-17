#pragma once

#include <map>
#include <string>
#include <vector>
#include <getopt.h>
#include <iostream>
#include <functional>

#include <algorithm> // For std::max
#include <iomanip>   // For std::setw and std::left
#include <sstream>   // For stringstream

// Class to represent a single configuration option
class ConfigOption
{
public:
    std::string name;
    std::string description;
    std::string value;
    std::string desc_value;
    bool requires_argument;

    // Default constructor (required by std::map)
    ConfigOption() : name(""), value(""), description(""), desc_value(""), requires_argument(false) {}

    // Constructor for initializing a configuration option
    ConfigOption(const std::string &name, const std::string &default_value, const std::string &description, const std::string &desc_value, bool requires_argument)
        : name(name), description(description), value(default_value), desc_value(desc_value), requires_argument(requires_argument) {}

    // Set value from string, converting to correct type if necessary
    void set_value_from_string(const std::string &input)
    {
        value = input;
    }
};

// Class to hold and manage all configuration options
class ConfigParser
{
public:
    std::map<std::string, ConfigOption> options;

    // Constructor initializes default options
    ConfigParser()
    {
        // Initialize configuration options with name, description, default value, and whether it requires an argument

        // Program execution options (FLAG, DEFAULT_VALUE, DESCRIPTION, DESC_VALUE, REQUIRES_ARGUMENT)
        add_option("q", "", "Executes a simulation with the given QASM file", "FILE_PATH", true);
        add_option("qs", "", "Executes a simulation with the given QASM string", "QASM_STR", true);
        add_option("j", "", "Executes a simulation with the given json file with Qiskit experiment Qobj", "FILE_PATH", true);
        add_option("js", "", "Executes a simulation with the given json string", "QOBJ_STR", true);
        add_option("t", "", "Runs the testing benchmarks for the specific index provided", "INDEX", true);
        add_option("a", "", "Runs all testing benchmarks", "", false);

        // Circuit execution options
        add_option("shots", "1024", "Number of shots", "SHOTS", true);
        add_option("backend", "cpu", "Simulation backend to use", "BACKEND", true);
        add_option("sim_method", "sv", "Simulation method", "METHOD", true);
        add_option("basis", "", "Run test benchmark with basis gates", "", false);
        add_option("disable_fusion", "", "Disable gate fusion (enabled by default)", "", false);
        add_option("seed", "", "Set random seed", "INT", true);

        // Hardware options
        add_option("tensorcore", "", "Use tensor cores", "", false);
        add_option("threads", "-1", "Number of OMP threads", "NUM_THREADS", true);
        add_option("AVX512", "", "Use AVX512", "", false);
        add_option("matrixcore", "", "Use MatrixCore", "", false);

        // Noise model options
        add_option("noise_model", "", "Noise model", "FILE_PATH", true);
        add_option("layout", "", "Path to json mapping logical qubits to device qubits, used when constructing the DM-Sim noise gates. (default: \"\", no mapping)", "FILE_PATH", true);
        add_option("layout_str", "", "String denoting qubit layout. Format maps logical qubits (lq, 0...n_qubits) to physical qubits. Format is lq0=pq0;lq1=pq1...", "LAYOUT_STR", true);

        // Initial and resulting state file options
        add_option("init_file", "", "Initial statevector/density matrix to enable simulation reuse/checkpointing", "FILE_PATH", true);
        add_option("init_format", "", "Format of initial state, default is same as sim_method", "FILE_PATH", true);
        add_option("dump_file", "", "Path to dump binary statevector/density matrix result", "FILE_PATH", true);

        // Helper options
        add_option("backend_list", "", "Print the list of available simulation backends", "", false);
        add_option("metrics", "", "Print the metrics of the circuit", "", false);
        add_option("verbose", "", "Verbose simulation trace", "", false);
        add_option("fidelity", "", "Run both DM-Sim and SV-Sim and report the state fidelity", "", false);
        add_option("h", "", "Prints the help message", "", false);
    }

    // Add a configuration option
    void add_option(const std::string &name, const std::string &default_value, const std::string &description, const std::string &desc_value, bool requires_argument)
    {
        options[name] = ConfigOption(name, default_value, description, desc_value, requires_argument);
    }

    // Get the value of a configuration option
    std::string get_value(const std::string &name) const
    {
        // std::cout << "name: " << name << std::endl;
        return options.at(name).value;
    }

    // Parse the command-line arguments
    void parse_command_line_arguments(int argc, char *argv[])
    {
        // Create long options array
        std::vector<struct option> long_options;
        std::string short_options;

        for (const auto &pair : options)
        {
            const ConfigOption &opt = pair.second;
            long_options.push_back({opt.name.c_str(), opt.requires_argument ? required_argument : no_argument, nullptr, 0});
        }
        long_options.push_back({0, 0, 0, 0}); // End of options

        int opt;
        int option_index = 0;

        // Parse arguments using getopt_long
        while ((opt = getopt_long(argc, argv, "", long_options.data(), &option_index)) != -1)
        {

            const char *opt_name = long_options[option_index].name;

            if (opt_name)
            {
                std::string option_key = opt_name;

                // Handle cases where required arguments are missing
                if (options[option_key].requires_argument && (optarg == nullptr || (optarg && std::string(optarg).substr(0, 2) == "--")))
                {
                    std::cerr << "Missing required argument for option: " << option_key << std::endl;
                    exit(1); // Exit on missing required argument
                }

                // Set value for flags or options with arguments
                std::string value = optarg ? std::string(optarg) : "true"; // Use optarg if available, or "true" for flags

                if (options.count(option_key))
                {
                    options[option_key].set_value_from_string(value); // Update the ConfigParser option with the parsed value
                }
            }
        }
    }

    // Function to print out the configurations with dynamic alignment and a header
    void print_help() const
    {
        // Calculate the flag width dynamically
        size_t flag_width = calculate_flag_width(options, true);

        // Set the total line width (e.g., 80 characters) and calculate the description width
        const size_t total_width = 80;
        const size_t description_width = total_width - flag_width;

        // Print the header
        std::cout << std::left << std::setw(flag_width) << "Option"
                  << "Description" << "\n";
        std::cout << std::string(flag_width + description_width, '-') << "\n";

        // Helper lambda to print a single option with formatting
        auto print_option = [&](const ConfigOption &opt)
        {
            // Format the flag and argument part
            std::stringstream flag_output;
            flag_output << "--" << opt.name;
            if (opt.requires_argument)
            {
                flag_output << " <" << opt.desc_value << ">";
            }

            // Print the flag/argument part, and wrap the description
            std::cout << std::left << std::setw(flag_width) << flag_output.str()
                      << wrap_description(opt.description, description_width, flag_width) << "\n";
        };

        // Loop through all options and print them
        for (const auto &pair : options)
        {
            const ConfigOption &opt = pair.second;
            print_option(opt);
        }
    }

    // Function to print out the configurations in a tabular format
    void print_configs() const
    {
        // Calculate the flag width dynamically (for better alignment)
        size_t flag_width = calculate_flag_width(options, false);

        // Set the total line width (e.g., 60 characters) and calculate the value width
        const size_t total_width = 30;
        const size_t value_width = total_width - flag_width;

        // Print the header with centered text
        std::cout << center_text("Flag", flag_width)
                  << center_text("Value", value_width) << "\n";

        // Print the header underline
        std::cout << std::string(total_width, '-') << "\n";

        // Helper lambda to print a single option with formatting
        auto print_option = [&](const ConfigOption &opt)
        {
            // Format the flag and argument part
            std::stringstream flag_output;
            flag_output << opt.name;

            // Set value output (either a real value, "true/false" for flags)
            std::string value_output = opt.requires_argument ? opt.value : opt.value.empty() ? "false"
                                                                                             : "true";

            // Print the formatted row with centered columns
            std::cout << center_text(flag_output.str(), flag_width)
                      << center_text(value_output, value_width) << "\n";

            // Print a horizontal line after each row
            std::cout << std::string(total_width, '-') << "\n";
        };

        // Loop through all options and print them in a tabular format with horizontal lines
        for (const auto &pair : options)
        {
            const ConfigOption &opt = pair.second;
            print_option(opt);
        }
    }

private:
    // Helper function to wrap descriptions to a certain width with indentation
    std::string static wrap_description(const std::string &desc, size_t width, size_t indent)
    {
        std::stringstream wrapped;
        size_t pos = 0;
        size_t line_start = 0;

        while (line_start < desc.size())
        {
            size_t end_pos = line_start + width;

            // Try to break at the last space before the width limit, or at the width limit if no space
            if (end_pos < desc.size() && desc[end_pos] != ' ')
            {
                end_pos = desc.find_last_of(' ', end_pos);
                if (end_pos == std::string::npos || end_pos <= line_start)
                {
                    end_pos = std::min(line_start + width, desc.size());
                }
            }

            wrapped << desc.substr(line_start, end_pos - line_start) << "\n";

            // Move to the start of the next line and add indentation
            line_start = end_pos + 1;
            if (line_start < desc.size())
            {
                wrapped << std::string(indent, ' ');
            }
        }

        return wrapped.str();
    }
    // Function to dynamically calculate the maximum flag width
    size_t static calculate_flag_width(const std::map<std::string, ConfigOption> &options, bool include_value = true)
    {
        size_t max_width = 0;

        // Find the maximum width of "--flag <value>" combination
        for (const auto &pair : options)
        {
            const ConfigOption &opt = pair.second;
            std::stringstream flag_output;
            flag_output << "--" << opt.name;
            if (opt.requires_argument && include_value)
            {
                flag_output << " <" << opt.desc_value << ">";
            }
            max_width = std::max(max_width, flag_output.str().length());
        }

        return max_width + 2; // Add some padding
    }

    // Helper function to center text in a given width
    std::string center_text(const std::string &text, size_t width) const
    {
        if (text.length() >= width)
        {
            return text; // No need to center if text is longer than the width
        }
        size_t left_padding = (width - text.length()) / 2;
        size_t right_padding = width - text.length() - left_padding;
        return std::string(left_padding, ' ') + text + std::string(right_padding, ' ');
    }
};