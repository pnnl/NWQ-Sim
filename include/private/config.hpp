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
    inline bool CUSTOMIZED_SP = false; 

    inline std::string DEVICE_CONFIG_PATH = "/people/lixi149/NWQ-Sim/data/device/";
    inline std::string DEVICE_CONFIG_FILE = "dummy_ibmq12";
    inline std::string CUSTMOIZED_GATES_FILE = "customized_basis_gates.json";

    inline nlohmann::json backend_config = {};

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

            if(CUSTOMIZED_SP){
                j.at("CUSTMOIZED_GATES_FILE").get_to(CUSTMOIZED_GATES_FILE);
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
           
            if (j.find("CUSTOMIZED_SP") != j.end())
                j.at("CUSTOMIZED_SP").get_to(CUSTMOIZED_GATES_FILE);
           
            if (j.find("DEVICE_CONFIG_PATH") != j.end())
                j.at("DEVICE_CONFIG_PATH").get_to(DEVICE_CONFIG_PATH);

            if (j.find("DEVICE_CONFIG_FILE") != j.end())
                j.at("DEVICE_CONFIG_FILE").get_to(DEVICE_CONFIG_FILE);

            if(j.find("CUSTMOIZED_GATES_FILE") != j.end())
                j.at("CUSTMOIZED_GATES_FILE").get_to(CUSTMOIZED_GATES_FILE);

        }
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

    inline void Load(const std::string &filename = "../default_config.json")
    {
        LoadConfigFromFile(filename, true);
        if (ENABLE_NOISE)
            readConfigFile();
    }

    inline void Update(const std::string &filename)
    {
        LoadConfigFromFile(filename, false);

        if (ENABLE_NOISE)
            readConfigFile();
    }

}