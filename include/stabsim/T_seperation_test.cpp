#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../backendManager.hpp"
#include "../state.hpp"
#include "../circuit.hpp"
#include "../nwq_util.hpp"

//#include "T_seperation.hpp"

namespace NWQSim
{

    void circuit_reverse(std::shared_ptr<Circuit>& circuit)
    {
        auto new_circ = std::make_shared<Circuit>(circuit->num_qubits());
        std::vector<Gate> gates = circuit->get_gates();
        IdxType num_gates = gates.size();
        int a; //qubit
        int b; //control

        for(int i = num_gates-1; i > -1; i--)
        {
            auto gate = gates[i];
            a = gate.qubit;
            b = gate.ctrl;
            if(gate.op_name == OP::S)
            {
                new_circ->S(a);
            }
            else if(gate.op_name == OP::H)
            {
                new_circ->H(a);
            }
            else if(gate.op_name == OP::CX)
            {
                new_circ->CX(b, a);
            }
            else
                std::cerr << "Unsupported gate in T transpilation!" << std::endl;
        }

        circuit = new_circ;
    }

    //Main compilation function for T pushthrough using tableau simulation
    void T_passthrough(std::shared_ptr<Circuit>& circuit, std::shared_ptr<QuantumState>& T_tab, std::shared_ptr<Circuit>& M_circ, int reps = 1)
    {
        IdxType n = circuit->num_qubits();

        std::vector<Gate> gates = circuit->get_gates();
        IdxType num_gates = gates.size();

        int a; //qubit
        int b; //control
        //As we loop through the circuit, keep track of the status of T
        for(int i = num_gates-1; i > -1; i--)
        {
            auto gate = gates[i];
            a = gate.qubit;
            b = gate.ctrl;
            if(gate.op_name == OP::T)
            {
                //Create a new Z stabilizer for the T gate
                std::string new_row;
                for(int i = 0; i < n; i++)
                {
                    new_row.push_back('I');
                }
                new_row[a] = 'Z';

                T_tab->add_stabilizer(new_row);

                std::cout << "T stab " << new_row << std::endl;

            }
            else if(gate.op_name == OP::S)
            {
                T_tab->apply_gate("S", a);
                M_circ->S(a);
            }
            else if(gate.op_name == OP::H)
            {
                T_tab->apply_gate("H", a);
                M_circ->H(a);
            }
            else if(gate.op_name == OP::CX)
            {
                T_tab->apply_gate("CX", a, b);
                M_circ->CX(b, a);
            }
            else
                std::cerr << "Unsupported gate in T transpilation!" << std::endl;
        }//End for loop of original circuit
        
        
        std::vector<std::shared_ptr<QuantumState>> P; //Empty vector of tableaus
        std::string tempStab;
        
        //Initialize P with the last row of the T tableau
        std::string backend = "CPU";
        std::string sim_method = "stab";
        double timer = 0;
        auto temp = BackendManager::create_state(backend, n, sim_method); //Temp tableau
        P.push_back(temp);
        P[0]->delete_all_rows();
        P[0]->remove_destabilizers();
        int T_rows = T_tab->get_num_rows();
        P[0]->add_stabilizer((T_tab->get_stabilizer_line(T_rows-1)).first);

        bool commutes = false;

        /*Process T tableaus*/

        std::cout << "---- T tableau before processing. Rows: " << T_rows << " ----" << std::endl;
        T_tab->print_res_state();

        //Sort through the T tableau in reverse order to undo our first reverse parse of the Cliff+T circuit
        for(int i = T_rows-2; i > -1; i--)
        {
            tempStab = (T_tab->get_stabilizer_line(i)).first;

            //Check all of the tableaus we have right now and see if any of them commute with the T line, starting with the first tableau
            for(int i = 0; i < P.size(); i++)
            {
                //Check if it commutes with a tableau
                //If it does, append it to that existing tableau
                if(P[i]->check_commutation(tempStab))
                {
                    P[i]->add_stabilizer(tempStab);
                    commutes = true;
                }   
                std::cout << "Temp stab " << tempStab << " commutes with P" << i << std::endl;
            }
            //If the stabilizer doesn't commute with any of the existing tableaus, make a new tableau
            if(!commutes)
            {

                auto newTab = BackendManager::create_state(backend, n, sim_method);
                P.push_back(newTab);
                P.back()->delete_all_rows();
                P.back()->remove_destabilizers();
                P.back()->add_stabilizer(tempStab);
                std::cout << "Temp stab " << tempStab << " does not commute. New P state ---" << std::endl;
                P.back()->print_res_state();
            }
            commutes = false; //reset commutation flag   
        }
        std::cout << "T tableau splitting finished.\n" << std::endl;
        std::cout << "T tableaus after splitting:" << std::endl;
        for(int i = P.size()-1; i > -1; i--)
        {
            std::cout << "P[" << i << "]: " << std::endl;

            P[i]->print_res_state();
        }


        //Push through any repeating stabilizers that may result in Clifford gates
        //Start from the last P, and push right. i.e. P1, P2, P3, check P3 for clifford gates first, then P2, then P1. If found in P3, push it through P2 and P1
        for(int i = P.size()-1; i > -1; i--)
        {
            //Fill in the stabilizer map
            std::unordered_map<std::string, int> stab_map;
            stab_map = P[i]->stabilizer_count();
            //For each recurring stabilizer
            for (const auto& pair : stab_map)
            {
                //Every 2 T gates is an S gate
                int num_s_gates = pair.second/2;
                if(pair.second > 1)
                {
                    std::cout << pair.first << " OCCURS " << pair.second << " TIMES " << std::endl;
                    P[i]->remove_repetitions(pair.first, num_s_gates);
                    
                    //push S gates through next tableaus
                    for(int j = i; j > -1; j--)
                    {
                        int P_rows = P[j]->get_num_rows();
                        P[j]->add_stabilizer(pair.first);

                        //Check the commutation with every line of the tableau
                        //If a line doesn't commute, rowsum that line with the temp stabilizer
                        for(int k = 0; k < P_rows; k++)
                        {
                            if(!(P[j]->check_row_commutation(pair.first, k)))
                            {
                                P[j]->rowsum(k, P_rows);
                            }
                        }
                        P[j]->remove_stabilizer(P_rows);
                    }
                }
                //Once all of the repeating stabilizers in a P tableau are pushed through,
                //apply S to the measurement tableau
                //The original T gate was applied at the qubit where z = 1
                int target_qubit = pair.first.find('Z');
                for(int num = 0; num < num_s_gates; num++)
                {
                    M_circ->S(target_qubit);
                }
                
            } //Done one P tableau
        } //Done all P tableaus

        T_tab->delete_all_rows();
        //Fill in the T tableau with the newly reduced P tableaus
        for(int i = 0; i < P.size(); i++)
        {
            std::vector<std::string> stabs = P[i]->get_stabilizers();
            std::cout << "First stabilizer at " << i << ": " << stabs[0] << std::endl;
            for(int j = 0; j < stabs.size(); j++)
            {
                T_tab->add_stabilizer(stabs[j]);
            }
        }

        /*Process is done, M and T have been seperated and returned*/

    } //End T_passthrough
}//End namespace


int main(){

    // Create a circuit with 2 qubits
    int n_qubits = 4;
    int shots = 10;

    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

    circuit->H(3);
    circuit->CX(2,0);
    circuit->CX(0,1);
    circuit->H(3);
    circuit->CX(0,2);
    circuit->CX(0,2);
    circuit->T(1);
    circuit->CX(3,2);
    circuit->CX(2,0);
    circuit->H(3);
    circuit->T(2);
    circuit->CX(0,1);
    circuit->CX(1,0);
    circuit->CX(3,2);
    circuit->CX(1,0);

    //Measurement circuit will be filled in the passthrough function
    auto M_circ = std::make_shared<NWQSim::Circuit>(n_qubits);

    std::string backend = "CPU";
    std::string sim_method = "stab";
    double timer = 0;
    
    /*Create T and Measurement Tableaus with only stabilizers. T starts empty, M starts as identity.*/
    auto T_tab = BackendManager::create_state(backend, n_qubits, sim_method);
    T_tab->delete_all_rows();
    T_tab->remove_destabilizers();

    T_passthrough(circuit, T_tab, M_circ);

    /*T Tableau*/
    std::cout << "---- T tableau -----" << std::endl;
    T_tab->print_res_state();
    /*T Tableau*/

    /*Measurement Tableau*/
    circuit_reverse(M_circ);

    //After all M_tab gates and additional Clifford gates from optimization
    //Simulate the M circuit evolution after construction
    auto M_tab = BackendManager::create_state(backend, n_qubits, sim_method);
    M_tab->sim(M_circ, timer);
    std::cout << "---- M tableau -----" << std::endl;

    M_tab->print_res_state();
    std::cout << "M tableau measurement results: " << (M_tab->measure_all(10))[0] << std::endl;
    /*Measurement Tableau*/

    return 0;
}
 