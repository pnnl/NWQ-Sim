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

    void T_separate(std::shared_ptr<QuantumState>& T_tab, std::vector<std::shared_ptr<QuantumState>>& P, int n)
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
                    {
                        P[0]->add_stabilizer(tempStab, phase);
                        break;
                    }
                    //else it commutes and we continue to the next tableau
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
                    {
                        P[j+1]->add_stabilizer(tempStab, phase);
                        break;
                    }
                }
            }
        }

        // // std::cout <<"T tableau splitting finished.\n" << std::endl;
        

        // std::cout <<"T tableaus after splitting:" << std::endl;
        // for(int i = 0; i < P.size(); i++)
        // {
        //     // std::cout <<"P[" << i << "]: " << std::endl;

        //     P[i]->print_res_state();
        // }
    }

    void T_optimize(std::vector<std::shared_ptr<QuantumState>>& P, std::shared_ptr<QuantumState>& M, int n)
    {
        //Push through any repeating stabilizers that may result in Clifford gates
        //Start from the last P, and push right. i.e. for (P1, P2, P3) check P3 for clifford gates first, then P2, then P1. If found in P1, push through P2 and P3.
        for(int i = P.size()-1; i > -1; i--)
        // for(int i = 0; i < P.size(); i++)
        {
            //Fill in the stabilizer map for repeating eigth rotations and remove repitions
            std::unordered_map<std::string, std::pair<int, int>> stab_map;

            P[i]->stabilizer_count(stab_map);
            // std::cout << "-----" << std::endl;
            // std::cout << "P[" << i << "]" << std::endl;
            // std::cout << "-----" << std::endl;

            //For each recurring stabilizer
            for(const auto& pair : stab_map)
            {
                int num_rotations = pair.second.second;
                
                //Check first that there is more than one stabilizer and they didn't all cancel out
                if((abs(num_rotations) > 1))
                {
                    std::string stabilizer = pair.first;
                    //If there were an odd number of stabilizers remaining, add back the one that doesn't cancel out

                    // std::cout << num_rotations << " repetitions of " << stabilizer << std::endl;
                    //8 1/8 gates cancel out
                    //Don't push through the odd number gate that was added back (if there were an odd number of stabilizers)
                    int rotations = num_rotations % 8;
                    rotations = (rotations/2);
                    // std::cout << "Rotations to push through for " << stabilizer << ": " << rotations << std::endl;

                    //Positive rotation gates come out of more positive T rotations than negative
                    if(rotations > 0)
                    {
                        //Apply rotation to the measurement tableau for every T gate that didn't cancel earlier
                        for(int num = 0; num < rotations; num++)
                        {
                            //Push an S gate through remaining tableaus in P order (to the right/to the end of the circuit)
                            for(int j = i+1; j < P.size(); j++)
                            // for(int j = i-1; j > -1 ; j--)
                            {
                                int P_rows = P[j]->get_num_rows();
                                P[j]->add_stabilizer(stabilizer);
                                // std::cout << "Added: " << stabilizer << std::endl;

                                //Check the commutation with every line of the tableau
                                //If a line doesn't commute, rowsum that line with the temp 'S' stabilizer
                                for(int k = 0; k < P_rows; k++)
                                {
                                    if(!(P[j]->check_row_commutation(stabilizer, k)))
                                    {
                                        P[j]->i_rowsum(k, P_rows);
                                        // std::cout << "Rowsum done at " << k << std::endl;
                                    }
                                }
                                P[j]->remove_stabilizer(P_rows);
                                // std::cout << "Removed: " << stabilizer << std::endl;
                                // std::cout <<"Stabilizer removed" << std::endl;
                            }

                            //Push S gate through the measurement circuit
                            int M_rows = M->get_num_rows();
                            M->add_stabilizer(stabilizer);

                            //Check the commutation with every line of the tableau
                            //If a line doesn't commute, rowsum that line with the temp 'S' stabilizer
                            for(int k = 0; k < M_rows; k++)
                            {
                                if(!(M->check_row_commutation(stabilizer, k)))
                                {
                                    M->i_rowsum(k, M_rows);
                                    // std::cout << "Applying rowsum on M at row " << k << std::endl;
                                }
                            }
                            M->remove_stabilizer(M_rows);
                        }
                    }
                    //Negative quarter rotation gates that come out of repeated TDG or similar
                    else if(rotations < 0)
                    {
                        //Apply rotation to the measurement tableau for every S gate that didn't cancel earlier
                        for(int num = rotations; num < 0; num++)
                        {
                            //Push an S gate through remaining tableaus in P order (to the right/to the end of the circuit)
                            for(int j = i+1; j < P.size(); j++)
                            {
                                int P_rows = P[j]->get_num_rows();
                                P[j]->add_stabilizer(stabilizer, 1);

                                //Check the commutation with every line of the tableau
                                //If a line doesn't commute, rowsum that line with the temp 'S' stabilizer
                                for(int k = 0; k < P_rows; k++)
                                {
                                    if(!(P[j]->check_row_commutation(stabilizer, k)))
                                    {
                                        P[j]->i_rowsum(k, P_rows);
                                    }
                                }
                                P[j]->remove_stabilizer(P_rows);
                                // std::cout <<"Stabilizer removed" << std::endl;
                            }

                            //Push S gate through the measurement circuit
                            int M_rows = M->get_num_rows();
                            M->add_stabilizer(stabilizer, 1);

                            //Check the commutation with every line of the tableau
                            //If a line doesn't commute, rowsum that line with the temp 'S' stabilizer
                            for(int k = 0; k < M_rows; k++)
                            {
                                if(!(M->check_row_commutation(stabilizer, k)))
                                {
                                    M->i_rowsum(k, M_rows);
                                    // std::cout << "Applying rowsum on M at row " << k << std::endl;
                                }
                            }
                            M->remove_stabilizer(M_rows);
                        }
                    } //If there are gates that need to be pushed
                } //If statement to check if there is more than 1 stabilizer
            } //End for loop of all stabilizers

            //If removing Clifford gates empties the tableau, then delete the tableau
            if(P[i]->get_num_rows() == 0)
            {
                P.erase(P.begin()+i);
            }
        }
    }

    void T_process(std::shared_ptr<Circuit>& circuit, std::shared_ptr<Circuit>& M_circ, std::shared_ptr<QuantumState>& T_tab, std::chrono::duration<long long, std::ratio<1, 1000000>>& proc_time)
    {
        IdxType n = M_circ->num_qubits();
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
            switch(gate.op_name)
            {
                case(OP::T):
                {
                    T_count++;
                    //Create a new Z stabilizer for the T gate
                    std::string new_row(n, 'I');
                    new_row[a] = 'Z';

                    T_tab->add_stabilizer(new_row);

                    // std::cout <<"T stab " << new_row << std::endl;
                    break;

                }
                case(OP::TDG):
                {
                    T_count++;
                    //Create a new Z stabilizer for the TDG gate
                    std::string new_row(n, 'I');
                    new_row[a] = 'Z';

                    T_tab->add_stabilizer(new_row, 1); //Phase for TDG

                    // std::cout <<"TDG stab " << new_row << std::endl;
                    break;

                }
                case(OP::S):
                {
                    T_tab->apply_gate("S", a);
                    M_circ->S(a);
                    break;
                }
                case(OP::H):
                {
                    T_tab->apply_gate("H", a);
                    M_circ->H(a);
                    break;
                }
                case(OP::CX):
                {
                    a = gate.ctrl;
                    b = gate.qubit;
                    T_tab->apply_gate("CX", a, b);
                    M_circ->CX(a, b);
                    break;
                }
                default:
                    std::cerr << "Unsupported gate in T transpilation!" << std::endl;
                    break;
            }//End Switch
        }//End for loop of original circuit
        // std::cout << "------\n\n\n TCount: " << T_count << " \n\n\n------" << std::endl;
        auto proc_end = std::chrono::high_resolution_clock::now();
        proc_time = std::chrono::duration_cast<std::chrono::microseconds>(proc_end - proc_start);
    }


    //Main compilation function for T pushthrough using tableau simulation
    void T_passthrough(std::shared_ptr<QuantumState>& T_tab, std::shared_ptr<QuantumState>& M_tab, std::ofstream& outfile, std::chrono::duration<long long, std::ratio<1, 1000000>> proc_time, int reps = 1)
    {
        
        IdxType n = M_tab->get_qubits();

        std::vector<std::shared_ptr<QuantumState>> P;
        int initial_T_count = T_tab->get_num_rows();
        if(initial_T_count)
            T_separate(T_tab, P, n);
        int P_rows = 0;
        int P_temp;
        for(int j = 0; j < P.size(); j++)
        {
            P_rows+= P[j]->get_num_rows();
        }
        int rep_print = reps;

        T_tab->delete_all_rows(); //T tableau is processed so it can be deleted for use later


        //Start the optimization defined by the number of reps
        std::chrono::microseconds opt_time(0);
        for(int i = 0; i < reps; i++)
        {
            P_temp = 0;
            // std::cout << "------\n\n Rep: " << i+1 << "\n\n------" << std::endl;
            if(P_rows)
            {
                auto opt_start = std::chrono::high_resolution_clock::now();
                T_optimize(P, M_tab, n);
                auto opt_end = std::chrono::high_resolution_clock::now();
                opt_time += std::chrono::duration_cast<std::chrono::microseconds>(opt_end - opt_start);
            }
            else
            {
                std::cout << "No T gates left to optimize." << std::endl;
                rep_print = i;
                break;
            }
            for(int j = 0; j < P.size(); j++)
            {
                P_temp += P[j]->get_num_rows();
            }
            //If the number of rows is the same after optimizing, end the optimization
            if(P_rows == P_temp)
            {
                std::cout << "T optimization done after " << i+1 << " reps." << std::endl;
                rep_print = i+1;
                break;
            }
            P_rows = P_temp;
        }



        //Processing is done, so we can combine into one tableau if needed
        // std::cout << "---- T tableau -----" << std::endl;
        for(int i = 0; i < P.size(); i++)
        {
            // std::cout << "---- P Tableau: " << i << " -----" << std::endl;
            // P[i]->print_res_state();
        }
        // std::cout << "---- End T tableau -----" << std::endl;
        T_combine(T_tab, P);

        outfile << proc_time.count() / 1000000.0 << std::endl;
        outfile << opt_time.count() / 1000000.0 << std::endl;
        outfile << initial_T_count << std::endl;
        outfile << T_tab->get_num_rows() << std::endl;
        outfile << rep_print << std::endl;

        /*Process is done, M and T have been seperated and returned*/

    } //End T_passthrough
}//End namespace