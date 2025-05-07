#pragma once

#include "../state.hpp"

#include "../nwq_util.hpp"
#include "../gate.hpp"
#include "../circuit.hpp"
#include "../config.hpp"
#include "../private/exp_gate_declarations_host.hpp"

#include "../circuit_pass/fusion.hpp"
#include "../private/macros.hpp"
#include "../private/sim_gate.hpp"

#include <random>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <vector>

#include "tensor.h"

namespace NWQSim
{
    class TN_CPU : public QuantumState
    {

    public:
        TN_CPU(IdxType _n_qubits) : QuantumState(SimType::TN)
        {
            // Initialize CPU side
            n_qubits = _n_qubits;

            dim = (IdxType)1 << (n_qubits);
            half_dim = (IdxType)1 << (n_qubits - 1);
            sv_size = dim * (IdxType)sizeof(ValType);
            n_cpu = 1;

            // CPU side initialization
            SAFE_ALOC_HOST(sv_real, sv_size);
            SAFE_ALOC_HOST(sv_imag, sv_size);
            memset(sv_real, 0, sv_size);
            memset(sv_imag, 0, sv_size);

            // State-vector initial state [0..0] = 1
            sv_real[0] = 1.;
            cpu_mem += sv_size * 4;

            SAFE_ALOC_HOST(m_real, sv_size + sizeof(ValType));
            memset(m_real, 0, sv_size + sizeof(ValType));
            rng.seed(Config::RANDOM_SEED);

            // ITensor MPS Initialization
	        auto sites = itensor::SpinHalf(int(n_qubits),{"ConserveQNs=", false});
    	    auto state = itensor::InitState(sites,"Up");
	        auto network = itensor::MPS(state);
            auto all_zeros = network;

        }

        ~TN_CPU()
        {
            // Release for CPU side
            SAFE_FREE_HOST(sv_real);
            SAFE_FREE_HOST(sv_imag);
            SAFE_FREE_HOST(m_real);
            SAFE_FREE_HOST(results);
        }

        void reset_state() override
        {
            // Reset CPU input & output
            memset(sv_real, 0, sv_size);
            memset(sv_imag, 0, sv_size);
            memset(m_real, 0, sv_size + sizeof(ValType));
            // State Vector initial state [0..0] = 1
            sv_real[0] = 1.;
	    
            // MPS initial state |00..0>
            sites = itensor::SpinHalf(int(n_qubits),{"ConserveQNs=", false});
            auto state = itensor::InitState(sites,"Up");
            network = itensor::MPS(state);
            all_zeros = network;
            network.position(1);
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
            IdxType origional_gates = circuit->num_gates();
            std::vector<SVGate> gates = fuse_circuit_sv(circuit);
            IdxType n_gates = gates.size();
            assert(circuit->num_qubits() == n_qubits);
            double sim_time;
            cpu_timer sim_timer;
            sim_timer.start_timer();
	    
            // Set Gauge of MPS to right-canonical
            network.position(1);
         

            simulation_kernel(gates);

            sim_timer.stop_timer();
            sim_time = sim_timer.measure();

            if (Config::PRINT_SIM_TRACE)
            {
                printf("\n============== TN-Sim ===============\n");
                printf("n_qubits:%lld, n_gates:%lld, sim_gates:%lld, ncpus:%lld, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_cpu:%.3lf MB\n",
                       n_qubits, origional_gates, n_gates, n_cpu, sim_time, 0.,
                       sim_time, cpu_mem / 1024 / 1024, cpu_mem / 1024 / 1024);
                printf("=====================================\n");
            }

            //=========================================
        }

        IdxType *get_results() override
        {
            return results;
        }

        IdxType measure(IdxType qubit) override
        {
            throw std::runtime_error("Not implemented");
            //M_GATE(qubit);
            //return results[0];
        }

        IdxType *measure_all(IdxType repetition) override
        {
            MA_GATE(repetition);
	    std::cout<<"measure_all was called"<<std::endl;
            return results;
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
        IdxType sv_size;
        IdxType dim;
        IdxType half_dim;
        IdxType n_cpu;

        // CPU arrays
        ValType *sv_real;
        ValType *sv_imag;
        ValType *m_real;

        //
        itensor::SpinHalf sites;
        itensor::MPS network,all_zeros;

        IdxType *results = NULL;

        // Random
        std::mt19937 rng;
        std::uniform_real_distribution<ValType> uni_dist;

        // CPU memory usage
        ValType cpu_mem;

        virtual void simulation_kernel(const std::vector<SVGate> &gates)
        {
            auto start = std::chrono::steady_clock::now();
            int n_gates = gates.size();
            for (int i = 0; i < n_gates; i++)
            {
                auto g = gates[i];

                if (g.op_name == OP::C1)
                {
                    C1_GATE(g.gm_real, g.gm_imag, g.qubit);
                }
                else if (g.op_name == OP::C2)
                {
                    C2_GATE(g.gm_real, g.gm_imag, g.ctrl, g.qubit);
                }
                else if (g.op_name == OP::RESET)
                {
                    RESET_GATE(g.qubit);
                }
                else if (g.op_name == OP::M)
                {
                    M_GATE(g.qubit);
                }
                else if (g.op_name == OP::MA)
                {
                    MA_GATE(g.qubit);
                }
                else if (g.op_name == OP::EXPECT)
                {
                    ObservableList *o = (ObservableList *)(g.data);
                    EXPECT_GATE(o);
                }
                else
                {
                    std::cout << "Unrecognized gates" << std::endl
                              << OP_NAMES[g.op_name] << std::endl;
                    std::logic_error("Invalid gate type");
                }

#ifdef PURITY_CHECK
                Purity_Check(g, i);
#endif
            }
            if (Config::PRINT_SIM_TRACE)
            {
                std::cout << std::endl;
            }
        }

        //============== C1 Gate ================
        // Arbitrary 1-qubit gate
        virtual void C1_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit)
        {
            // std::cout<<"C1"<<std::endl;
            int site = qubit + 1;

            //Initialize empty C1 Gate tensor
            auto j = sites(site);
            auto gate = itensor::ITensor(prime(j),j);

            // Set values of C1 gate tensor
            gate.set(1,1,std::complex<double>(gm_real[0],gm_imag[0]));
            gate.set(1,2,std::complex<double>(gm_real[1],gm_imag[1]));
            gate.set(2,1,std::complex<double>(gm_real[2],gm_imag[2]));
            gate.set(2,2,std::complex<double>(gm_real[3],gm_imag[3]));

            // Move center of orthongality to qubit that will be contracted with gate
            network.position(site);
            // PrintData(network);
            auto temp = gate * network(site);
            // PrintData(temp);
            temp.noPrime();
            // PrintData(gate);
            
            // Set the MPS site to new values
            network.set(site,temp);
            network.orthogonalize();
            
            
            // PrintData(network);
            

        }

        //============== C2 Gate ================
        // Arbitrary 2-qubit gate
        virtual void C2_GATE(const ValType *gm_real, const ValType *gm_imag,
                             const IdxType qubit0, const IdxType qubit1)
        {

            assert(qubit0 != qubit1); // Non-cloning

            itensor::ITensor full_mat;
            full_mat = network(1);
            for(auto i = 1 ; i <= n_qubits-1 ; i++){
                full_mat *= network(i+1);
            }
            // PrintData(full_mat);

            int site0 = qubit0 + 1;
            int site1 = qubit1 + 1;
            auto i = sites(site0);
            auto j = sites(site1);
            auto gate = itensor::ITensor(dag(i),dag(j),prime(i),prime(j));

            gate.set(1,1,1,1,std::complex<double>(gm_real[0],gm_imag[0]));
            gate.set(1,2,1,1,std::complex<double>(gm_real[1],gm_imag[1]));
            gate.set(2,1,1,1,std::complex<double>(gm_real[2],gm_imag[2]));
            gate.set(2,2,1,1,std::complex<double>(gm_real[3],gm_imag[3]));
                                                                        
            gate.set(1,1,1,2,std::complex<double>(gm_real[4],gm_imag[4]));
            gate.set(1,2,1,2,std::complex<double>(gm_real[5],gm_imag[5]));
            gate.set(2,1,1,2,std::complex<double>(gm_real[6],gm_imag[6]));
            gate.set(2,2,1,2,std::complex<double>(gm_real[7],gm_imag[7]));
                                                                    
            gate.set(1,1,2,1,std::complex<double>(gm_real[8],gm_imag[8]));
            gate.set(1,2,2,1,std::complex<double>(gm_real[9],gm_imag[9]));
            gate.set(2,1,2,1,std::complex<double>(gm_real[10],gm_imag[10]));
            gate.set(2,2,2,1,std::complex<double>(gm_real[11],gm_imag[11]));
                                                                        
            gate.set(1,1,2,2,std::complex<double>(gm_real[12],gm_imag[12]));
            gate.set(1,2,2,2,std::complex<double>(gm_real[13],gm_imag[13]));
            gate.set(2,1,2,2,std::complex<double>(gm_real[14],gm_imag[14]));
            gate.set(2,2,2,2,std::complex<double>(gm_real[15],gm_imag[15]));

            auto q_dist = std::abs(qubit0 - qubit1);
            // PrintData(gate);
            //
            // std::cout<<"C2"<<std::endl;
            // PrintData(network);
            // Non-adjacent 2-qubit gates using MPOs
            //
            auto method = true;
            if(std::abs(qubit0 - qubit1) != 1){
            // if(1){

            
                // std::cout<<"non local"<<std::endl;
                

                auto gate_MPO = itensor::MPO(n_qubits);
                auto [u,s,v] = itensor::svd(gate,{i,prime(i)},{j,prime(j)});
                auto sv = s*v;
                // PrintData(u);
                // PrintData(s);
                // PrintData(v);

                // Sequential Contraction Method based on Quantinuum pytket
                if(method){
                    itensor::ITensor propagating_bond;
                    itensor::Index lindex;
                    network.position(site0);
                    // PrintData(network);
                    auto site0_contract = u*network(site0);
                    // PrintData(site0_contract);
                    auto [U,S,V] = itensor::svd(site0_contract,{i,prime(i),leftLinkIndex(network,site0)});
                    network.set(site0,U);
                    network.position(site0);
                    propagating_bond =S*V;
                
                    // PrintData(U);
                    // PrintData(propagating_bond);

                    lindex = commonInds(U,S)[0];
                    // std::cout<<"leftlim "<<leftLim(network)<<std::endl;
                    // std::cout<<"rightlim "<<rightLim(network)<<std::endl;
            

                    for(auto i = site0+1 ; i < site1; i++ ){
                        // PrintData(lindex);
                        auto i_contract = propagating_bond * network(i);
                        // PrintData(network(i));
                        // PrintData(i_contract);
                        
                        auto [u,s,v] = itensor::svd(i_contract,{sites(i),lindex});
                        lindex = commonInds(u,s)[0];
                        // PrintData(u);
                        network.set(i,u);
                        network.position(i);
                        propagating_bond = s*v;
                        // PrintData(propagating_bond);
                        // std::cout<<"leftlim "<<leftLim(network)<<std::endl;
                        // std::cout<<"rightlim "<<rightLim(network)<<std::endl;
                
        
                    }
                    // PrintData(sv);
                    // PrintData(network(site1));
                    auto site1_contract = sv*propagating_bond*network(site1);
    
                    network.set(site1,site1_contract);
                    
                    // PrintData(network);
                    
                    // std::cout<<"leftlim "<<leftLim(network)<<std::endl;
                    // std::cout<<"rightlim "<<rightLim(network)<<std::endl;
            



                }
                else{
                    // MPO Method, Almost Works?
                    //
                    // WARNING Normalization breaks!!
                    auto lusv = commonInds(u,s)[0];

                    itensor::Index lright;
                    itensor::Index lleft;
                    network.position(1);

                    for (auto i = 1; i<=n_qubits;i++){
                        if(not((i == site0) || (i == site1))){
                            if(i==1){
                                lright = itensor::Index(4,"Link");
                                auto diag1 = itensor::delta(sites(i),prime(sites(i)));
                                auto diag2 = itensor::delta(lright);
                                auto diag3 = toDense(diag2);
                                auto temp = diag1 * diag3;
                                gate_MPO.set(i,temp);
                            }
                            else if( i==n_qubits){
                                lleft = lright;
                                auto diag1 = itensor::delta(sites(i),prime(sites(i)));
                                auto diag2 = itensor::delta(lleft);
                                auto diag3 = toDense(diag2);
                                auto temp = diag1 * diag3;
                                gate_MPO.set(i,temp);
                            }
                            else{
                                lleft = lright;
                                lright = itensor::Index(4,"Link");
                                auto diag1 = itensor::delta(sites(i),prime(sites(i)));
                                auto diag2 = itensor::delta(lleft,lright);
                                auto diag3 = toDense(diag2);
                                auto temp = diag1 * diag3;
                                gate_MPO.set(i,temp);
                            }
                        }
                        else if(i==site0){
                            if(not(i==1 or i==n_qubits)){
                                lleft = lright;
                                lright = lusv;
                                gate_MPO.set(site0, u*itensor::delta(lleft));
                            }
                            else{
                                lright = lusv;
                                gate_MPO.set(site0, u);
                            }
                        }
                        else if(i==site1){
                            if(not(i==1 or i==n_qubits)){
                                lleft = lright;
                                auto lsv = itensor::delta(lleft,lusv);
                                lright = itensor::Index(4,"Link");
                                gate_MPO.set(site1, s*v*lsv*itensor::delta(lright));
                            }
                            else{
    
                                lleft = lright;
                                auto lsv = itensor::delta(lleft,lusv);
                                gate_MPO.set(site1, s*v*lsv);
                            }
                        }
                    }
                    auto full_gate = gate_MPO(1);
                    for(auto i = 1 ; i <= n_qubits-1 ; i++){
                        full_gate *= gate_MPO(i+1);
                    }
                    // PrintData(full_gate);
                    // PrintData(full_mat);
                    // PrintData(gate_MPO);
                    network = applyMPO(gate_MPO,network);
                    network.normalize();
                    // PrintData(network);
                }
                network.noPrime();
                network.orthogonalize();
                // PrintData(network);
                // std::cout<<"leftlim "<<leftLim(network)<<std::endl;
                // std::cout<<"rightlim "<<rightLim(network)<<std::endl;

                full_mat = network(1);
                for(auto i = 1 ; i <= n_qubits-1 ; i++){
                    full_mat *= network(i+1);
                }
                // PrintData(full_mat);
            }
            else{
                // std::cout<<"local"<<std::endl;
                network.position(site0);
                auto contract_location = network(site0)*network(site1);
                auto temp = gate*contract_location;
                    temp.noPrime();
                auto [u,s,v] = itensor::svd(temp,itensor::inds(network(site0)),{"Cutoff=", 0.0, "MaxDim=", 10});	    
                network.set(site0, u);
                network.set(site1, s*v);
                network.orthogonalize();
            }
        }

        //============== C4 Gate ================
        // Arbitrary 4-qubit gate
        virtual void C4_GATE(const ValType *gm_real, const ValType *gm_imag,
                             const IdxType qubit0, const IdxType qubit1,
                             const IdxType qubit2, const IdxType qubit3)
        {
            throw std::runtime_error("Not implemented");
        }

        virtual void M_GATE(const IdxType qubit)
        {
            throw std::runtime_error("Not implemented");
        }

        //============== MA Gate (Measure all qubits in Pauli-Z) ================
        virtual void MA_GATE(const IdxType repetition)
        {
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
            memset(results, 0, sizeof(IdxType) * repetition);

            
 
            // Move to right orthogonal gauge for improved scaling before sampling
            network.position(1);
            // PrintData(network);
            // std::cout<<"Net norm "<<norm(network)<<std::endl;
         
            // std::cout<<"0000 "<<innerC(network,all_zeros)<<std::endl;
            // itensor::ITensor full_mat;
            // full_mat = network(1);
            // for(auto i = 1 ; i <= n_qubits-1 ; i++){
            //     full_mat *= network(i+1);
            // }
            // PrintData(full_mat);
            // std::cout<<"0000 "<<innerC(network,all_zeros)<<std::endl;
            
            for (IdxType i = 0; i < repetition; i++)
            {
        
                auto network_ =  network;
                IdxType res = 0;

                for(IdxType j = 1;j <= n_qubits; j++){

                    ValType r = uni_dist(rng);

                    if(j==n_qubits){

     
                        auto rdm = network_(1) * prime(dag(network_(1)));
      
                        auto p_si = std::real(eltC(rdm,1,1));
                    
                        if (r >= p_si){
                            res |= static_cast<IdxType>(1) << (j-1);
                            }
                        results[i] = res;
                        }
                    else{
                 
                        auto rdm = prime(dag(network_(n_qubits-j+1)),"Link")*network_(n_qubits-j+1);
                        
                      
                        
                        for (auto k = n_qubits-j; k >= 1 ; k--){
                            
                            rdm *= network_(k);
                            if (k==1){
                                rdm *= prime(dag(network_(k)));
                            }
                            else{
                                rdm *= prime(dag(network_(k)),"Link");
                            }
                        }
                        
                        auto p_si = std::real(eltC(rdm,1,1));

                        auto site = sites(j);
                        auto si = itensor::ITensor(site);

                        if (r <= p_si){
       
                            si.set(1,1.);
                            si.set(2,0.);

                            auto temp = si * network_(1);

                            temp *= network_(2);


                            auto temp_net = itensor::MPS(n_qubits-j);
                            for (int i = 1;i<=(n_qubits-j);i++){
                            
                                temp_net.ref(i) = network_(i+1);}
    
                            temp_net.set(1,temp);
                            network_ = temp_net;

                            if (p_si > 0){
                                network_ /= std::sqrt(p_si);
                                }


                        }
                        else{
                           
                            si.set(1,0.);
                            si.set(2,1.);
      
                            auto temp = si * network_(1);
                
              
                            temp *= network_(2);
                 



                            auto temp_net = itensor::MPS(n_qubits-j);
              
                            for (int i = 1;i<=(n_qubits-j);i++){
                             
                                temp_net.ref(i) = network_(i+1);}
                    
                            temp_net.set(1,temp);
                          
                            network_ = temp_net;
                  
                    

                            res |= static_cast<IdxType>(1) << (j-1);

                            if (p_si != 1){
                            network_ /= std::sqrt(1-p_si);
                            }
                    }
                    }
		        }
            }
        }
        virtual double EXPECT_C4_GATE(const ValType *gm_real, const ValType *gm_imag, IdxType qubit0, IdxType qubit1, IdxType qubit2, IdxType qubit3, IdxType mask)
        {
            throw std::runtime_error("Not implemented");
        }

        virtual double EXPECT_C2_GATE(const ValType *gm_real, const ValType *gm_imag, IdxType qubit0, IdxType qubit1, IdxType mask)
        {
            throw std::runtime_error("Not implemented");
        }

        virtual double Expect_C0(IdxType mask)
        {
            throw std::runtime_error("Not implemented");
        }


        virtual void EXPECT_GATE(ObservableList *o)
        {
            throw std::runtime_error("Not implemented");
        }
        //============== Reset ================
        virtual void RESET_GATE(const IdxType qubit)
        {
            throw std::runtime_error("Not implemented");
        }

        //============== Purity Check  ================
        void Purity_Check(SVGate g, const IdxType t)
        {
            ValType purity = 0;
            for (IdxType i = 0; i < (((IdxType)1 << n_qubits)); i++)
                purity += ((sv_real[i] * sv_real[i]) + (sv_imag[i] * sv_imag[i]));
            if (fabs(purity - 1.0) > ERROR_BAR)
            {
                printf("Purity Check fails after Gate-%lld=>%s(ctrl:%lld,qubit:%lld) with %lf\n", t, OP_NAMES[g.op_name], g.ctrl, g.qubit, purity);
            }
        }

        virtual ValType *get_real() const override { return sv_real; };
        virtual ValType *get_imag() const override { return sv_imag; };
    };

} // namespace NWQSim
