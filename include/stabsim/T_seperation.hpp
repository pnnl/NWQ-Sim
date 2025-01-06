#include "../backendManager.hpp"
#include "../state.hpp"
#include "../circuit.hpp"
#include "../nwq_util.hpp"

namespace NWQSim
{

    void T_optimize(std::shared_ptr<Circuit>& T_circ, int reps = 1)
    {

    }

    //Main compilation function for T pushthrough using tableau simulation
    void T_passthrough(std::shared_ptr<Circuit>& circuit, std::shared_ptr<Circuit>& T_circ, std::shared_ptr<Circuit>& M_circ, int reps = 1)
    {
        IdxType n = circuit->num_qubits();

        std::string backend = "CPU";
        std::string sim_method = "stab";
        double timer = 0;
        
        /*Create T and Measurement Tableaus with only stabilizers. T starts empty, M starts as identity.*/
        auto T_tab = BackendManager::create_state(backend, n, sim_method);
        T_tab->delete_all_rows();
        T_tab->has_destabilizers = false;

        auto M_tab = BackendManager::create_state(backend, n, sim_method);
        M_tab->remove_destabilizers();
        /**/
    

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
                //Create a 
                std::string new_row;
                for(int i = 0; i < n; i++)
                {
                    new_row.push_back('I');
                }
                new_row[a] = 'Z';
                T_tab->add_stabilizer(new_row);
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
        //Simulate the M circuit evolution after construction
        M_tab->sim(M_circ, timer);
        
        std::vector<std::shared_ptr<QuantumState>> P; //Empty vector of tableaus
        std::string tempStab;
        
        //Initialize P with the last row of the T tableau
        auto temp = BackendManager::create_state(backend, n, sim_method); //Temp tableau
        P.push_back(temp);
        P[0]->delete_all_rows();
        P[0]->remove_destabilizers();
        int T_rows = T_tab->get_num_rows();
        P[0]->add_stabilizer((T_tab->get_stabilizer_line(T_rows-1)).first);

        bool commutes = false;

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
            }
            //If the stabilizer doesn't commute with any of the existing tableaus, make a new tableau
            if(!commutes)
            {
                auto newTab = BackendManager::create_state(backend, n, sim_method);
                P.push_back(newTab);
            }
            commutes = false; //reset commutation flag   
        }

        std::unordered_map<std::string, int> stab_map;

        //Push through any repeating stabilizers that may result in Clifford gates
        //Start from the last P, and push right. i.e. P1, P2, P3, check P3 for clifford gates first, then P2, then P1. If found in P3, push it through P2 and P1
        for(int i = P.size()-1; i > -1; i--)
        {
            //Fill in the stabilizer map
            stab_map = P[i]->stabilizer_count();
            //For each recurring stabilizer
            for (const auto& pair : stab_map)
            {
                //Every 2 T gates is an S gate
                int num_s_gates = pair.second/2;
                if(pair.second > 1)
                {
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
        //Simulate the final Clifford Circuit
        M_tab->sim(M_circ, timer);

        T_tab->delete_all_rows();
        //Fill in the T tableau with the newly reduced P tableaus
        for(int i = 0; i < P.size(); i++)
        {
            std::vector<std::string> stabs = P[i]->get_stabilizers();
            for(int j = 0; j < stabs.size(); j++)
            {
                T_tab->add_stabilizer(stabs[j]);
            }
        }

        //Process is done, M and T have been seperated and returned
    }
}