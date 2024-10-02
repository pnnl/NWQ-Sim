#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <time.h>

namespace NWQSim::Config
{
    inline bool PRINT_SIM_TRACE = true;
    inline bool ENABLE_NOISE = false;
    inline bool ENABLE_FUSION = true;
    inline bool ENABLE_TENSOR_CORE = false;
    inline bool ENABLE_MATRIX_CORE = false;
    inline bool ENABLE_AVX512 = false; // AVX512 is not supported yet
    inline int OMP_NUM_THREADS = -1;   // -1 means using the default number of threads

    inline int RANDOM_SEED = time(NULL); // 5489

    inline std::string device_noise_file = "";
    inline std::string device_layout_file = "";
    inline std::string device_layout_str = "";

    using IdxType = long long;
    using ValType = double;

    static void printConfig(IdxType i_proc, const std::string &sim_backend)
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

}