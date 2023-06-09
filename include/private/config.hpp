#pragma once

#include "nlohmann/json.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

namespace NWQSim::Config
{
    inline bool PRINT_SIM_TRACE;

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
        }
        else
        {
            if (j.find("PRINT_SIM_TRACE") != j.end())
                j.at("PRINT_SIM_TRACE").get_to(PRINT_SIM_TRACE);
        }
    }

    inline void printConfig()
    {
        std::cout << "Running with the following configuration:" << std::endl;
        std::cout << "PRINT_SIM_TRACE: " << (PRINT_SIM_TRACE ? "True" : "False") << std::endl;

        std::cout << std::endl
                  << std::endl;
    }

    inline void Load(const std::string &filename = "../default_config.json")
    {
        LoadConfigFromFile(filename, true);
        printConfig();
    }

    inline void Update(const std::string &filename)
    {
        LoadConfigFromFile(filename, false);
    }

}
