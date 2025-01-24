#pragma once

#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../backendManager.hpp"
#include "../state.hpp"
#include "../circuit.hpp"
#include "../nwq_util.hpp"

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

    void T_optimize(std::shared_ptr<QuantumState>& T_tab, std::shared_ptr<Circuit>& M_circ, int reps)
    {
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
            std::cout << i << std::endl;
            tempStab = (T_tab->get_stabilizer_line(i)).first;

            //Check all of the tableaus we have right now and see if any of them commute with the T line, starting with the first tableau
            for(int i = P.size()-1; i > -1; i--)
            {
                //Check if it commutes with a tableau
                //If it does, append it to that existing tableau
                if(P[i]->check_commutation(tempStab))
                {
                    P[i]->add_stabilizer(tempStab);
                    commutes = true;
                    std::cout << "Temp stab " << tempStab << " commutes with P" << i << std::endl;
                }  
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
        T_tab->delete_all_rows();

        std::cout << "T tableaus after splitting:" << std::endl;
        for(int i = 0; i < P.size()-1; i++)
        {
            std::cout << "P[" << i << "]: " << std::endl;

            P[i]->print_res_state();
        }

        //Push through any repeating stabilizers that may result in Clifford gates
        //Start from the last P, and push right. i.e. for (P1, P2, P3) check P3 for clifford gates first, then P2, then P1. If found in P1, push through P2 and P3.
        for(int i = P.size()-1; i > -1; i--)
        {
            //Fill in the stabilizer map
            std::unordered_map<std::string, int> stab_map;
            P[i]->stabilizer_count(stab_map);
            //For each recurring stabilizer
            for(const auto& pair : stab_map)
            {
                //Every 2 T gates is an S gate
                int num_quarter_gates = pair.second/2;

                std::cout << pair.first<< " OCCURS " << pair.second << " TIMES " << std::endl;

                //Positive rotation gates come out of multiple T's
                if(num_quarter_gates > 0)
                {
                    std::cout << "In num quarter gates" << std::endl;

                    P[i]->remove_repetitions(pair.first, num_quarter_gates * 2);

                    //4 S gates is a full rotation around the bloch sphere and cancel out
                    num_quarter_gates = num_quarter_gates % 4;
                    if(num_quarter_gates != 0)
                    {
                        //Perform the push through routine for every S gate that doesn't cancel
                        for(int gate_num = 0; gate_num < num_quarter_gates; gate_num++)
                        {
                            //Push an S gate through remaining tableaus in P order (to the right/to the end of the circuit)
                            for(int j = i; j < P.size(); j++)
                            {
                                int P_rows = P[j]->get_num_rows();
                                P[j]->add_stabilizer(pair.first);

                                //Check the commutation with every line of the tableau
                                //If a line doesn't commute, rowsum that line with the temp 'S' stabilizer
                                for(int k = 0; k < P_rows; k++)
                                {
                                    if(!(P[j]->check_row_commutation(pair.first, k)))
                                    {
                                        P[j]->rowsum(k, P_rows);
                                    }
                                }
                                P[j]->remove_stabilizer(P_rows);
                            }
                            //Once all of the repeating stabilizers in a P tableau are pushed through,
                            //apply S to the measurement tableau.
                            //The original T gate was applied where z = 1 in the 'S' stabilizer.
                            for(int col= 0; col < pair.first.size(); col++)
                            {
                                if(pair.first[col] == 'Z')
                                {
                                    M_circ->S(target_qubit);
                                }
                                else if(pair.first[col] == 'Y')
                                {
                                    M_circ->RY(PI/2, target_qubit);
                                }
                                else if(pair.first[col] == 'X')
                                {
                                    M_circ->RX(PI/2, target_qubit);
                                }
                            }
                        }
                    }
                }
                //Negative quarter rotation gates that come out of repeated TDG or similar
                else if(num_quarter_gates < 0)
                {
                    std::cout << "In num sdg gates" << std::endl;

                    P[i]->remove_repetitions(pair.first, num_quarter_gates * -2);

                    //4 S gates is a full rotation around the bloch sphere and cancel out
                    num_quarter_gates = num_quarter_gates % 4;
                    if(num_quarter_gates != 0)
                    {
                        //Perform the push through routine for every S gate that doesn't cancel
                        for(int s_gates = 0; s_gates > num_quarter_gates; s_gates--)
                        {
                            //Push an SDG gate through remaining tableaus in P order (to the right/to the end of the circuit)
                            for(int j = i; j < P.size(); j++)
                            {
                                int P_rows = P[j]->get_num_rows();
                                //SDG
                                P[j]->add_stabilizer(pair.first, 1);

                                //Check the commutation with every line of the tableau
                                //If a line doesn't commute, rowsum that line with the temp 'S' stabilizer
                                for(int k = 0; k < P_rows; k++)
                                {
                                    if(!(P[j]->check_row_commutation(pair.first, k)))
                                    {
                                        P[j]->rowsum(k, P_rows);
                                    }
                                }
                                P[j]->remove_stabilizer(P_rows);
                            }
                            //Once all of the repeating stabilizers in a P tableau are pushed through,
                            //apply S to the measurement tableau.
                            //The original T gate was applied where z = 1 in the 'S' stabilizer.
                            for(int col= 0; col < pair.first.size(); col++)
                            {
                                if(pair.first[col] == 'Z')
                                {
                                    M_circ->SDG(target_qubit);
                                }
                                else if(pair.first[col] == 'Y')
                                {
                                    M_circ->RY(-PI/2, target_qubit);
                                }
                                else if(pair.first[col] == 'X')
                                {
                                    M_circ->RX(-PI/2, target_qubit);
                                }
                            }
                        }
                    }
                }
            }
        }
        //Fill in the T tableau with the newly reduced P tableaus
        for(int i = 0; i < P.size(); i++)
        {
            std::cout << "P[" << i << "]: " << std::endl;

            P[i]->print_res_state();
            std::vector<std::string> stabs = P[i]->get_stabilizers();
            std::cout << stabs.size() << " -- stabs.size()" << std::endl;

            for(int j = 0; j < stabs.size(); j++)
            {
                std::cout << "First stabilizer at " << i << ": " << stabs[j] << std::endl;
                T_tab->add_stabilizer(stabs[j]);
            }
        }
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
            else if(gate.op_name == OP::TDG)
            {
                //Create a new Z stabilizer for the T gate
                std::string new_row;
                for(int i = 0; i < n; i++)
                {
                    new_row.push_back('I');
                }
                new_row[a] = 'Z';

                T_tab->add_stabilizer(new_row, 1);

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

        for(int i = 0; i < reps; i++)
        {
            T_optimize(T_tab, M_circ, reps);
        }
        
        //Put the M circuit back in forward time after the new Clifford gates have been appended
        circuit_reverse(M_circ);

        /*Process is done, M and T have been seperated and returned*/

    } //End T_passthrough
}//End namespace