#pragma once

#include "private/nlohmann/json.hpp"
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
    inline bool ENABLE_MATRIX_CORE = false;
    inline bool ENABLE_AVX512 = false; // AVX512 is not supported yet
    inline int OMP_NUM_THREADS = -1;   // -1 means using the default number of threads

    inline int RANDOM_SEED = 5489;

    using IdxType = long long;
    using ValType = double;

    inline nlohmann::json backend_config = {};
    inline nlohmann::json layout = {};

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
            std::cout << std::left << std::setw(40) << "\033[1;33mSimulation Backend:\033[0m"
                      << "\033[1;32m" << (sim_backend == "sv" ? "SVSim" : "DMSim") << "\033[0m" << std::endl;
            if (Config::PRINT_SIM_TRACE)
            {
                std::cout << std::endl;
            }
        }
    }

    /**
     * @brief Read json file that stores backend noise properties from input path. Output from the function in the python script device_config.py
     *        See an example usage in backedn_noise_validation.cpp
     *
     * @param path the absolute path path where json file locates. e.g., "./Data/"
     */
    inline void readDeviceConfig(const std::string &path)
    {
        std::ifstream f(path);
        if (f.fail())
            throw std::logic_error("Device config file not found at " + path);
        backend_config = nlohmann::json::parse(f);
        if (layout.empty())
        {
            IdxType n_qubits;
            backend_config["num_qubits"].get_to(n_qubits);
            for (IdxType n = 0; n < n_qubits; n++)
            {
                layout[std::to_string(n)] = (int)n;
            }
        }
    }

    inline std::string qindex(IdxType q)
    {
        return std::to_string(static_cast<IdxType>(
            layout[std::to_string(q)]));
    }

    inline void loadLayoutFile(const std::string &filename, bool required = true)
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

    inline void loadLayoutString(std::string s)
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

}