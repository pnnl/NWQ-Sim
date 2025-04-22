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

	    sites_ = itensor::SpinHalf(n_qubits)
	    psi_full_ = itensor::ITensor(sites_);

        }
        
        void reset_state() override
        {
	    psi_full_fill(0.0);

	    std::vector<long> idx(n_qubits,1);
	    psi_full_.set(idx, 1.0);
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
    	}

        void sim(std::shared_ptr<NWQSim::Circuit> circuit) override
        {
    	    // initialize the state
    	    reset_state();
    
    	    // temporarily use the fuse gates from sv_sim ToDo: Make a tn_fuse circuits that returns gates
    	    
    	    assert(circuit->num_qubits()==n_qubits);
    	    auto gates = fuse_circuit_sv(circuit);
    
    	    // apply all of the gates one by one 
    	    for (auto const& g : gates)
    	    {
    	        if(g.op_name==OP::C1)
    		        apply_one_site(g);
    		    else if(g.op_name==OP::C2)
    		        apply_two_site(g);
            }
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
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType)*repetition);

            const IdxType nstates = IdxType(1) << n_qubits;

            std::vector<double> probs(nstates);
            std::vector<itensor::IndexVal> iv(n_qubits);

            for(IdxType s = 0; s < nstates; ++s)
            {

            }
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

        void one_qubit_gate(const GateRecord& g)
        {
            // itensor starts at 1
            int site = g.qubit + 1
            auto I = sites_(site), Ip = prime(I);

            // itensor builds tensors by naming indices 
            // gate tensor
            itensor::ITensor G(I, Ip);
            
            // ToDo: Fill in the gate from the sv gate
            
            psi_full_ = psi_full_ * G;
            psi_full_.noPrime();
        }

        void two_qubit_gate(const GateRecord& g)
        {
            int i1 = g.ctrl + 1, i2 = g.qubit + 1;
            auto I1 = sites_(i1), I2 = sites_(i2);
            auto Ip1 = prime(I1), Ip2 = prime(I2);

            itensor::ITensor G(I1, I2, Ip1, Ip2);

            // ToDo: Fill in the gate from the sv gate
            
            psi_full_ = psi_full_ * G;
            psi_full_.noPrime();
        }
    };

} // namespace NWQSim
