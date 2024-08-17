#pragma once

#include "nlohmann/json.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

namespace NWQSim::Config
{
    inline bool PRINT_SIM_TRACE = true;
    inline bool ENABLE_NOISE = false;
    inline bool ENABLE_FUSION = true;
    inline bool ENABLE_TENSOR_CORE = false;
    inline bool DEVICE_READ = false;
    using IdxType = long long;
    using ValType = double;
    inline std::string DEVICE_CONFIG_PATH = "~/NWQ-Sim/data/device/";
    inline std::string DEVICE_CONFIG_FILE = "dummy_ibmq12";

    inline nlohmann::json backend_config = {};
    inline nlohmann::json layout = {};

    inline void LoadConfigFromFile(const std::string &filename, bool required = true)
    {
        std::ifstream i(filename);
        if (!i.is_open())
        {
            std::cout << "Could not open " << filename << std::endl;
            std::cout << "Using Default Configurations" << std::endl;
            return;
            // throw std::runtime_error("Could not open " + filename);
        }

        nlohmann::json j;
        try
        {
            i >> j;
        }
        catch (nlohmann::json::parse_error &e)
        {
            throw std::runtime_error(std::string("Error parsing JSON: ") + e.what());
        }

        if (required)
        {
            j.at("PRINT_SIM_TRACE").get_to(PRINT_SIM_TRACE);
            j.at("ENABLE_NOISE").get_to(ENABLE_NOISE);
            j.at("ENABLE_FUSION").get_to(ENABLE_FUSION);
            j.at("ENABLE_TENSOR_CORE").get_to(ENABLE_TENSOR_CORE);

            if (ENABLE_NOISE)
            {
                j.at("DEVICE_CONFIG_PATH").get_to(DEVICE_CONFIG_PATH);
                j.at("DEVICE_CONFIG_FILE").get_to(DEVICE_CONFIG_FILE);
            }
        }
        else
        {
            if (j.find("PRINT_SIM_TRACE") != j.end())
                j.at("PRINT_SIM_TRACE").get_to(PRINT_SIM_TRACE);

            if (j.find("ENABLE_NOISE") != j.end())
                j.at("ENABLE_NOISE").get_to(ENABLE_NOISE);

            if (j.find("ENABLE_FUSION") != j.end())
                j.at("ENABLE_FUSION").get_to(ENABLE_FUSION);

            if (j.find("ENABLE_TENSOR_CORE") != j.end())
                j.at("ENABLE_TENSOR_CORE").get_to(ENABLE_TENSOR_CORE);

            if (j.find("DEVICE_CONFIG_PATH") != j.end())
                j.at("DEVICE_CONFIG_PATH").get_to(DEVICE_CONFIG_PATH);

            if (j.find("DEVICE_CONFIG_FILE") != j.end())
                j.at("DEVICE_CONFIG_FILE").get_to(DEVICE_CONFIG_FILE);
        }
    }
    inline void LoadLayoutFromFile(const std::string &filename, bool required = true)
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

    inline void printConfig(IdxType i_proc, const std::string &sim_backend)
    {
        if (i_proc == 0)
        {
            std::cout << "\033[1;34m\nRunning simulation with the following configuration:\033[0m" << std::endl;
            std::cout << std::left << std::setw(40) << "\033[1;33mSim Trace:\033[0m"
                      << "\033[1;32m" << (PRINT_SIM_TRACE ? "Enabled" : "Disabled") << "\033[0m" << std::endl;
            std::cout << std::left << std::setw(40) << "\033[1;33mFusion:\033[0m"
                      << "\033[1;32m" << (ENABLE_FUSION ? "Enabled" : "Disabled") << "\033[0m" << std::endl;
            std::cout << std::left << std::setw(40) << "\033[1;33mNoise:\033[0m"
                      << "\033[1;32m" << (ENABLE_NOISE ? "Enabled" : "Disabled") << "\033[0m" << std::endl;
            std::cout << std::left << std::setw(40) << "\033[1;33mTensor Core:\033[0m"
                      << "\033[1;32m" << (ENABLE_TENSOR_CORE ? "Enabled" : "Disabled") << "\033[0m" << std::endl;
            std::cout << std::left << std::setw(40) << "\033[1;33mDevice Config Path:\033[0m"
                      << "\033[1;32m" << (DEVICE_CONFIG_PATH.empty() ? "Not Set" : DEVICE_CONFIG_PATH) << "\033[0m" << std::endl;
            std::cout << std::left << std::setw(40) << "\033[1;33mDevice Config File:\033[0m"
                      << "\033[1;32m" << (DEVICE_CONFIG_FILE.empty() ? "Not Set" : DEVICE_CONFIG_FILE) << "\033[0m" << std::endl;
            std::cout << std::left << std::setw(40) << "\033[1;33mSimulation Backend:\033[0m"
                      << "\033[1;32m" << (sim_backend == "sv" ? "SVSim" : "DMSim") << "\033[0m" << std::endl;
            if (Config::PRINT_SIM_TRACE)
            {
                std::cout << std::endl;
            }
        }
    }

    /**
     * @brief Read json file that stores backend noise properties. Output from the function in the python script device_config.py
     *        See an example usage in backedn_noise_validation.cpp
     *
     * @param file_dir the path where json file locates. e.g., "./Data/"
     * @param backend_name the backend name, e.g., "ibmq_toronto"
     */
    inline void readConfigFile()
    {
        std::string appendix = ".json";
        std::string path = DEVICE_CONFIG_PATH + DEVICE_CONFIG_FILE + appendix;
        std::ifstream f(path);
        if (f.fail())
            throw std::logic_error("Device config file not found at " + path);
        backend_config = nlohmann::json::parse(f);
        
    }

    /**
     * @brief Read json file that stores backend noise properties from input path. Output from the function in the python script device_config.py
     *        See an example usage in backedn_noise_validation.cpp
     *
     * @param path the absolute path path where json file locates. e.g., "./Data/"
     */
    inline void readDeviceConfig(const std::string& path) {
        std::ifstream f(path);
        if (f.fail())
            throw std::logic_error("Device config file not found at " + path);
        backend_config = nlohmann::json::parse(f);
        DEVICE_READ = true;
    }

    inline void Load(const std::string &filename = "../default_config.json")
    {
        
        LoadConfigFromFile(filename, true);
        if (ENABLE_NOISE && !DEVICE_READ) {
            readConfigFile();
        }
        if (layout.empty()) {
            IdxType n_qubits;
            backend_config["num_qubits"].get_to(n_qubits);
            for (IdxType n = 0; n < n_qubits; n++) {
                layout[std::to_string(n)] = (int)n;
            }
        }
        
    }
    inline std::string qindex(IdxType q) {
        return std::to_string(static_cast<IdxType>(
            layout[std::to_string(q)]
        ));
    }
    inline void readLayoutString(std::string s) {
        size_t pos_start = 0;
        size_t pos_end;
        std::string token;
        while((pos_end = s.find(",", pos_start)) != std::string::npos) {
            token = s.substr (pos_start, pos_end - pos_start);
            size_t src_end = token.find("=");
            if (src_end == std::string::npos) {
                throw std::invalid_argument("Ill-formatted layout string\n");
            }
            std::string log_qb = token.substr(0, src_end);
            std::string phys_qb = token.substr(src_end+1);
            layout[log_qb] = std::stoll(phys_qb);
            pos_start = pos_end + 1;
        }
        token = s.substr (pos_start);
        if (token.length()) {
            size_t src_end = token.find("=");
            std::string log_qb = token.substr(0, src_end);
            std::string phys_qb = token.substr(src_end+1);
            layout[log_qb] = std::stoll(phys_qb);
        }
    }
    inline void Update(const std::string &filename)
    {
        LoadConfigFromFile(filename, false);

        if (ENABLE_NOISE && !DEVICE_READ)
            readConfigFile();
    }

}