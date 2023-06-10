#pragma once

#include "nlohmann/json.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

namespace NWQSim::Config
{
    inline bool PRINT_SIM_TRACE;
    inline bool ENABLE_NOISE;
    inline bool ENABLE_FUSION;
    inline std::string DEVICE_CONFIG_PATH;
    inline std::string DEVICE_CONFIG_FILE;

    inline nlohmann::json backend_config = {};

    inline void LoadConfigFromFile(const std::string &filename, bool required = true)
    {
        std::ifstream i(filename);
        if (!i.is_open())
        {
            throw std::runtime_error("Could not open " + filename);
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

            if (j.find("DEVICE_CONFIG_PATH") != j.end())
                j.at("DEVICE_CONFIG_PATH").get_to(DEVICE_CONFIG_PATH);

            if (j.find("DEVICE_CONFIG_FILE") != j.end())
                j.at("DEVICE_CONFIG_FILE").get_to(DEVICE_CONFIG_FILE);
        }
    }

    inline void printConfig(IdxType i_proc = 0)
    {
        if (i_proc == 0)
        {
            std::cout << "Running with the following configuration:" << std::endl;
            std::cout << "PRINT_SIM_TRACE: " << (PRINT_SIM_TRACE ? "True" : "False") << std::endl;
            std::cout << "ENABLE_NOISE: " << (ENABLE_NOISE ? "True" : "False") << std::endl;
            std::cout << "ENABLE_FUSION: " << (ENABLE_FUSION ? "True" : "False") << std::endl;

            std::cout << "DEVICE_CONFIG_PATH: " << DEVICE_CONFIG_PATH << std::endl;
            std::cout << "DEVICE_CONFIG_FILE: " << DEVICE_CONFIG_FILE << std::endl;

            std::cout << std::endl
                      << std::endl;
        }
    }

    /**
     * @brief Read json file that stores backend noise properties. Output from the function in the python script device_config.py
     *        See an example usage in backedn_noise_validation.cpp
     *
     * @param file_dir the path where json file locates. e.g., "./Data/"
     * @param backend_name the backend name, e.g., "ibmq_toronto"
     */
    void readConfigFile()
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