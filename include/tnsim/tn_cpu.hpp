#pragma once

#include "../state.hpp"

#include "../nwq_util.hpp"
#include "../gate.hpp"
#include "../circuit.hpp"
#include "../config.hpp"
#include "private/exp_gate_declarations_host.hpp"

#include "../circuit_pass/fusion.hpp"
#include "../private/macros.hpp"
#include "../private/sim_gate.hpp"

#include <random>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <vector>

#include "itensor/all.h"
#include <cmath>

namespace NWQSim
{
    class TN_CPU : public QuantumState
    {

    public:
        TN_CPU(IdxType _n_qubits) : QuantumState(SimType::TN)
        {
            // Initialize CPU side
            n_qubits = _n_qubits;

            n_cpu = 1;

            rng.seed(Config::RANDOM_SEED);
        }
        
        void reset_state() override
        {
            throw std::runtime_error("Not implemented");
	}

        void set_seed(IdxType seed) override
        {
            rng.seed(seed);
        }

        virtual void set_initial(std::string fpath, std::string format) override
        {
            throw std::runtime_error("Not implemented");
	}

        virtual void dump_res_state(std::string outpath) override
        {
            throw std::runtime_error("Not implemented");
	};

        void sim(std::shared_ptr<NWQSim::Circuit> circuit) override
        {
            throw std::runtime_error("Not implemented");
	}

        IdxType *get_results() override
        {
            throw std::runtime_error("Not implemented");
        }

        IdxType measure(IdxType qubit) override
        {
            throw std::runtime_error("Not implemented");
        }

        IdxType *measure_all(IdxType repetition) override
        {
            throw std::runtime_error("Not implemented");
        }

        virtual ValType get_exp_z(const std::vector<size_t> &in_bits) override
        {
            throw std::runtime_error("Not implemented");
        }

        virtual ValType get_exp_z() override
        {
            throw std::runtime_error("Not implemented");
        }

        void print_res_state() override
        {
            throw std::runtime_error("Not implemented");
        }

    protected:
        // n_qubits is the number of qubits
        IdxType n_qubits;
        IdxType n_cpu;

        IdxType *results = NULL;

        // Random
        std::mt19937 rng;
        std::uniform_real_distribution<ValType> uni_dist;

        // CPU memory usage
        ValType cpu_mem;
    };

} // namespace NWQSim
