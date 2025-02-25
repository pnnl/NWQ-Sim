#pragma once

#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>

#include "../../include/backendManager.hpp"
#include "../../include/state.hpp"
#include "../../include/circuit.hpp"
#include "../../include/nwq_util.hpp"

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
            if(gate.op_name == OP::S)
            {
                new_circ->S(a);
                // std::cout << "S(" << a << ")" << std::endl;
            }
            else if(gate.op_name == OP::SDG)
            {
                new_circ->SDG(a);
                // std::cout << "SDG(" << a << ")" << std::endl;
            } 
            else if(gate.op_name == OP::H)
            {
                new_circ->H(a);
                // std::cout << "H(" << a << ")" << std::endl;
            }
            else if(gate.op_name == OP::RX)
            {
                new_circ->RX(gate.theta, a);
                // std::cout << "RX(" << a << ")" << std::endl;
            }
            else if(gate.op_name == OP::RY)
            {
                new_circ->RY(gate.theta, a);
                // std::cout << "RY(" << a << ")" << std::endl;
            }
            else if(gate.op_name == OP::CX)
            {
                a = gate.ctrl;
                b = gate.qubit;
                new_circ->CX(a, b);
                // std::cout << "CX(" << a << ", " << b << ")" << std::endl;
            }
            else
                std::cerr << "Unsupported gate in circuit reverse!" << std::endl;
        }
        circuit = new_circ;
    }

    void T_combine(std::shared_ptr<QuantumState>& T_tab, std::vector<std::shared_ptr<QuantumState>>& P)
    {
        //Fill in the T tableau with the newly reduced P tableaus
        for(int i = 0; i < P.size(); i++)
        {
            // // std::cout <<"P[" << i << "]: " << std::endl;

            // P[i]->print_res_state();
            std::vector<std::string> stabs = P[i]->get_stabilizers();
            // // std::cout <<stabs.size() << " -- stabs.size()" << std::endl;

            for(int j = 0; j < stabs.size(); j++)
            {
                // // std::cout <<"First stabilizer at " << i << ": " << stabs[j] << std::endl;
                T_tab->add_stabilizer(stabs[j]);
            }
        }

    }

    void T_seperate(std::shared_ptr<QuantumState>& T_tab, std::vector<std::shared_ptr<QuantumState>>& P)
    {
        std::string tempStab;
        int phase;
        
        //Initialize P with the last row of the T tableau
        std::string backend = "CPU";
        std::string sim_method = "stab";
        double timer = 0;
        auto temp = BackendManager::create_state(backend, n, sim_method); //Temp tableau
        P.push_back(temp);
        P[0]->delete_all_rows();
        P[0]->remove_destabilizers();
        int T_rows = T_tab->get_num_rows();
        P[0]->add_stabilizer((T_tab->get_stabilizer_line(T_rows-1).first), (T_tab->get_stabilizer_line(T_rows-1).second));

        bool commutes = false;

        /*Process T tableaus*/

        // std::cout <<"---- T tableau before processing. Rows: " << T_rows << " ----" << std::endl;
        // T_tab->print_res_state();

        //Sort through the T tableau in reverse order to undo our first reverse parse of the Cliff+T circuit
        for(int i = T_rows-2; i > -1; i--) //Start at T_rows-2 since we manually started at T-1
        {
            // // std::cout <<i << std::endl;
            tempStab = (T_tab->get_stabilizer_line(i)).first;
            phase = (T_tab->get_stabilizer_line(i)).second;

            //Check all of the tableaus we have right now and see if any of them commute with the T line
            //We start with the last tableau and find the earliest one that commutes moving backwards
            for(int j = P.size()-1; j > -1; j--)
            {
                //Check if it commutes with a tableau
                if(P[j]->check_commutation(tempStab))
                {
                    if(j == 0)
                        P[0]->add_stabilizer(tempStab, phase);
                }
                //If it doesn't add it to the last tableau it commuted with
                else
                {
                    //If it didn't commute with the first tableau we checked, we have to make a new one
                    if(j == P.size()-1)
                    {
                        auto newTab = BackendManager::create_state(backend, n, sim_method);
                        P.push_back(newTab);
                        P.back()->delete_all_rows();
                        P.back()->remove_destabilizers();
                        P.back()->add_stabilizer(tempStab, phase);
                        break;
                    }
                    //If it didn't commute and we're not on the last tableau, add it to the previous tableau that
                    //it commuted with
                    else
                        P[j+1]->add_stabilizer(tempStab, phase);
                }
            }
        }
        T_tab->delete_all_rows(); //T tableau is processed so it can be deleted
    }

    void T_optimize(std::vector<std::shared_ptr<QuantumState>>& P, std::shared_ptr<Circuit>& M_circ, int n)
    {

        // // std::cout <<"T tableau splitting finished.\n" << std::endl;
        

        // std::cout <<"T tableaus after splitting:" << std::endl;
        // for(int i = 0; i < P.size(); i++)
        // {
        //     // std::cout <<"P[" << i << "]: " << std::endl;

        //     P[i]->print_res_state();
        // }

        //Push through any repeating stabilizers that may result in Clifford gates
        //Start from the last P, and push right. i.e. for (P1, P2, P3) check P3 for clifford gates first, then P2, then P1. If found in P1, push through P2 and P3.
        for(int i = P.size()-1; i > -1; i--)
        {
            //Fill in the stabilizer map for repeating eigth rotations
            std::unordered_map<std::string, int> stab_map;
            P[i]->stabilizer_count(stab_map);
            //For each recurring stabilizer
            for(const auto& pair : stab_map)
            {
                std::string repeated = pair.first;

                //Adds/subtracts the positive and negative rotations
                int total_counts = P[i]->stabilizer_reps(repeated);

                //Every 2 T rotations is a quarter gate
                int num_quarter_gates = total_counts/2;

                // // std::cout <<pair.first<< " OCCURS " << pair.second << " TIMES " << std::endl;

                //Remove the repeating stabilizers, but add one back if there is an odd number of gates. 
                //If there is an extra negative gate, makes sure to add the negative phase.
                P[i]->remove_repetitions(repeated, total_counts%2);

                //4 S gates is a full rotation around the bloch sphere and cancel out so we might be done
                num_quarter_gates = num_quarter_gates % 4;
                //Positive rotation gates come out of more positive T rotations than negative
                if(num_quarter_gates > 0)
                {
                    //Perform the push through routine for every S gate that didn't cancel earlier (up to 3 gates)
                    for(int gate_num = 0; gate_num < num_quarter_gates; gate_num++)
                    {
                        //Push an S gate through remaining tableaus in P order (to the right/to the end of the circuit)
                        for(int j = i; j < P.size(); j++)
                        {
                            int P_rows = P[j]->get_num_rows();
                            P[j]->add_stabilizer(repeated);

                            //Check the commutation with every line of the tableau
                            //If a line doesn't commute, rowsum that line with the temp 'S' stabilizer
                            for(int k = 0; k < P_rows; k++)
                            {
                                if(!(P[j]->check_row_commutation(repeated, k)))
                                {
                                    P[j]->rowsum(k, P_rows);
                                }
                            }
                            P[j]->remove_stabilizer(P_rows);
                            // std::cout <<"Stabilizer removed" << std::endl;
                        }
                        //Once all of the repeating stabilizers in a P tableau are pushed through,
                        //apply S to the measurement tableau.
                        //The original T gate was applied where z = 1 in the 'S' stabilizer.
                        for(int col = 0; col < n; col++)
                        {
                            if(pair.first[col] == 'Z')
                            {
                                M_circ->S(col);
                            }
                            else if(pair.first[col] == 'Y')
                            {
                                M_circ->RY(PI/2, col);
                            }
                            else if(pair.first[col] == 'X')
                            {
                                M_circ->RX(PI/2, col);
                            }
                        }
                    }
                }
                //Negative quarter rotation gates that come out of repeated TDG or similar
                else if(num_quarter_gates < 0)
                {
                    // std::cout <<"In num sdg gates" << std::endl;
                    
                    //Perform the push through routine for every SDG gate that doesn't cancel (1-3 SDG gates)
                    for(int s_gates = 0; s_gates > num_quarter_gates; s_gates--)
                    {
                        //Push an SDG gate through remaining tableaus in P order (to the right/to the end of the circuit)
                        for(int j = i; j < P.size(); j++)
                        {
                            int P_rows = P[j]->get_num_rows();
                            //SDG
                            P[j]->add_stabilizer(repeated, 1);

                            //Check the commutation with every line of the tableau
                            //If a line doesn't commute, rowsum that line with the temp 'S' stabilizer
                            for(int k = 0; k < P_rows; k++)
                            {
                                if(!(P[j]->check_row_commutation(repeated, k)))
                                {
                                    P[j]->rowsum(k, P_rows);
                                }
                            }
                            P[j]->remove_stabilizer(P_rows);
                        }
                        //Once all of the repeating stabilizers in a P tableau are pushed through,
                        //apply S to the measurement tableau.
                        //The original T gate was applied where z = 1 in the 'S' stabilizer.
                        for(int col = 0; col < n; col++)
                        {
                            if(pair.first[col] == 'Z')
                            {
                                M_circ->SDG(col);
                            }
                            else if(pair.first[col] == 'Y')
                            {
                                //Backwards so we apply S first for -pi/2
                                M_circ->RY(-PI/2, col);
                            }
                            else if(pair.first[col] == 'X')
                            {
                                M_circ->RX(-PI/2, col);
                            }
                        }
                    }
                }
            }
        }

    }


    //Main compilation function for T pushthrough using tableau simulation
    void T_passthrough(std::shared_ptr<Circuit>& circuit, std::shared_ptr<QuantumState>& T_tab, std::shared_ptr<Circuit>& M_circ, std::string outFile, int reps = 1)
    {

        std::ofstream outfile(outFile); // Open a file for writing

        
        IdxType n = circuit->num_qubits();

        std::vector<Gate> gates = circuit->get_gates();
        IdxType num_gates = gates.size();

        int a; //qubit
        int b; //control
        int T_count = 0;
        int T_rows = 0;

        auto proc_start = std::chrono::high_resolution_clock::now();
        //As we loop through the circuit, keep track of the status of T
        for(int i = num_gates-1; i > -1; i--)
        {
            auto gate = gates[i];
            a = gate.qubit;
            if(gate.op_name == OP::T)
            {
                T_count++;
                //Create a new Z stabilizer for the T gate
                std::string new_row = "";
                for(int i = 0; i < n; i++)
                {
                    new_row.push_back('I');
                }
                new_row[a] = 'Z';

                T_tab->add_stabilizer(new_row);

                std::cout <<"T stab " << new_row << std::endl;

            }
            else if(gate.op_name == OP::TDG)
            {
                T_count++;
                //Create a new Z stabilizer for the T gate
                std::string new_row = "";
                for(int i = 0; i < n; i++)
                {
                    new_row.push_back('I');
                }
                new_row[a] = 'Z';

                T_tab->add_stabilizer(new_row, 1);

                std::cout <<"TDG stab " << new_row << std::endl;

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
                a = gate.ctrl;
                b = gate.qubit;
                T_tab->apply_gate("CX", a, b);
                M_circ->CX(a, b);
            }
            else
                std::cerr << "Unsupported gate in T transpilation!" << std::endl;
        }//End for loop of original circuit
        std::cout << "------\n\n\n TCount: " << T_count << " \n\n\n------" << std::endl;

        T_seperate
        
        auto proc_end = std::chrono::high_resolution_clock::now();
        auto proc_time = std::chrono::duration_cast<std::chrono::microseconds>(proc_end - proc_start);
        
        auto opt_start = std::chrono::high_resolution_clock::now();

        int rep_print = reps;

        for(int i = 0; i < reps; i++)
        {
            T_rows = T_tab->get_num_rows();
            std::cout << "------\n\n" << T_rows << "\n\n------" << std::endl;
            if(T_rows)
                T_optimize(T_tab, M_circ, n);
            else
            {
                std::cout << "No T gates to optimize." << std::endl;
                rep_print = i;
                break;
            }
            //If the number of rows is the same after optimizing, end the optimization
            if(T_rows == (T_tab->get_num_rows()))
            {
                std::cout << "T optimization done at " << i << " reps." << std::endl;
                rep_print = i;
                break;
            }

        }

        auto opt_end = std::chrono::high_resolution_clock::now();
        auto opt_time = std::chrono::duration_cast<std::chrono::microseconds>(opt_end - opt_start);

        outfile << proc_time.count() << std::endl;
        outfile << opt_time.count() << std::endl;
        outfile << T_rows << std::endl;
        outfile << T_tab->get_num_rows() << std::endl;
        outfile << rep_print << std::endl;
        
        
        //Put the M circuit back in forward time after the new Clifford gates have been appended
        circuit_reverse(M_circ);

        /*Process is done, M and T have been seperated and returned*/

        outfile.close(); // Close the file
    } //End T_passthrough
}//End namespace