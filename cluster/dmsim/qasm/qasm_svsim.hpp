#ifndef QASM_ASSEMBLER
#define QASM_ASSEMBLER

#include <algorithm>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <regex>
#include <bits/stdc++.h>
#include <map>


#include "../src/util.h"
#include "../src/svsim_nvgpu_sin.cuh"
#include "qasm_util.hpp"

using namespace std;
using namespace NWQSim;

map<string, reg> list_qregs;
vector<reg> list_cregs;

vector<gate> list_custom_gates;
vector<string> buffered_measure;
void load_gate(string op, Simulation &sim, bool buffer_measure);

/**
 * @brief Load defined registers
 *
 * @param reg_line Register defination string
 * @return reg
 */
reg load_reg(string reg_line, int *offset)
{
    reg cur_reg;
    vector<string> reg_fields = split(split(reg_line, ' ')[1], '[');

    cur_reg.name = reg_fields[0];

    string reg_witdth = regex_replace(reg_fields[1], std::regex(" +$"), "");

    reg_witdth.pop_back(); // remove ;
    reg_witdth.pop_back(); // remove ]

    cur_reg.width = stoi(reg_witdth);

    if (reg_line.substr(0, 1) == "q")
    {
        cur_reg.offset = *offset;
        *offset += cur_reg.width;
    }
    else
    {
        cur_reg.offset = 0;
    }

    // printf("Reg name:%5s, Reg width:%2d, Reg offset:%2d\n", cur_reg.name.c_str(), cur_reg.width, cur_reg.offset);
    return cur_reg;
}

/**
 * @brief Load custom gate definations
 *
 * @param gate_line Header line of the defination
 * @param file File handler to load gate body
 * @return gate
 */
gate load_gate_defination(string gate_line, ifstream *file)
{
    // Initialize gate object
    gate cur_gate;

    // Spilit gate defination into name, parameter, and args
    vector<string> header_fields = split(gate_line, ' ');
    vector<string> name_paramters = split(header_fields[1], '(');
    cur_gate.name = name_paramters[0];
    if (name_paramters.size() > 1)
    {
        name_paramters[1].pop_back();
        cur_gate.params = split(name_paramters[1], ',');
    }

    cur_gate.args = split(header_fields[2], ',');

    // Load gate body
    string line;
    while (getline(*file, line))
    {
        if (line.substr(0, 1) == "{")
        {
            // Starting line of current gate defination
            continue;
        }
        else if (line.substr(0, 1) == "}")
        {
            // Reacing the end of current gate defination
            return cur_gate;
        }
        else
        {
            line.erase(remove(line.begin(), line.end(), '\r'), line.end()); // remove '\r'
            line.erase(remove(line.begin(), line.end(), ';'), line.end());  // remove ';'

            // Trim the leading/tailing white space of the operation string and append it to the list
            cur_gate.ops.push_back(regex_replace(line, std::regex("^ +| +$|( ) +"), "$1"));
        }
    }

    return cur_gate;
}

void custom_ccx(vector<IdxType> qubits, vector<ValType> params, Simulation &sim)
{
    sim.H(qubits[2]);
    sim.CX(qubits[1], qubits[2]);
    sim.TDG(qubits[2]);

    sim.CX(qubits[0], qubits[2]);
    sim.T(qubits[2]);

    sim.CX(qubits[1], qubits[2]);
    sim.T(qubits[1]);
    sim.TDG(qubits[2]);

    sim.CX(qubits[0], qubits[2]);
    sim.T(qubits[2]);

    sim.CX(qubits[0], qubits[1]);
    sim.T(qubits[0]);
    sim.TDG(qubits[1]);
    sim.H(qubits[2]);

    sim.CX(qubits[0], qubits[1]);
}

void append_gate(string gate_name, vector<IdxType> qubits, vector<ValType> params, Simulation &sim)
{
    transform(gate_name.begin(), gate_name.end(), gate_name.begin(), ::toupper);

    if (gate_name == "U")
        sim.U(params[0], params[1], params[2], qubits[0]);
    else if (gate_name == "U1")
        sim.U(0, 0, params[0], qubits[0]);
    else if (gate_name == "U2")
        sim.U(pi / 2, params[0], params[1], qubits[0]);
    else if (gate_name == "U3")
        sim.U(params[0], params[1], params[2], qubits[0]);
    else if (gate_name == "X")
        sim.X(qubits[0]);
    else if (gate_name == "Y")
        sim.Y(qubits[0]);
    else if (gate_name == "Z")
        sim.Z(qubits[0]);
    else if (gate_name == "H")
        sim.H(qubits[0]);
    else if (gate_name == "S")
        sim.S(qubits[0]);
    else if (gate_name == "SDG")
        sim.SDG(qubits[0]);
    else if (gate_name == "T")
        sim.T(qubits[0]);
    else if (gate_name == "TDG")
        sim.TDG(qubits[0]);
    else if (gate_name == "RX")
        sim.RX(params[0], qubits[0]);
    else if (gate_name == "RY")
        sim.RY(params[0], qubits[0]);
    else if (gate_name == "RZ")
        sim.RZ(params[0], qubits[0]);
    else if (gate_name == "CX")
        sim.CX(qubits[0], qubits[1]);
    else if (gate_name == "CY")
        sim.CY(qubits[0], qubits[1]);
    else if (gate_name == "CZ")
        sim.CZ(qubits[0], qubits[1]);
    else if (gate_name == "CH")
        sim.CH(qubits[0], qubits[1]);
    else if (gate_name == "CCX")
        custom_ccx(qubits, params, sim);
    else if (gate_name == "CRX")
        sim.CRX(params[0], qubits[0], qubits[1]);
    else if (gate_name == "CRY")
        sim.CRY(params[0], qubits[0], qubits[1]);
    else if (gate_name == "CRZ")
        sim.CRZ(params[0], qubits[0], qubits[1]);
    else if (gate_name == "CU")
        sim.CU(params[0], params[1], params[2], params[3], qubits[0], qubits[1]);
    else if (gate_name == "CU1")
        sim.CU(0, 0, params[0], 0, qubits[0], qubits[1]);
    else if (gate_name == "CU3")
        sim.CU(params[0], params[1], params[2], 0, qubits[0], qubits[1]);
    else if (gate_name == "RESET")
        sim.RESET(qubits[0]);
    else if (gate_name == "SWAP")
        sim.SWAP(qubits[0], qubits[1]);
    else if (gate_name == "SX")
        sim.SX(qubits[0]);

    else if (gate_name == "RI")
        sim.RI(params[0], qubits[0]);
    else if (gate_name == "P")
        sim.P(params[0], qubits[0]);

    else if (gate_name == "CS")
        sim.CS(qubits[0], qubits[1]);
    else if (gate_name == "CSDG")
        sim.CSDG(qubits[0], qubits[1]);
    else if (gate_name == "CT")
        sim.CT(qubits[0], qubits[1]);
    else if (gate_name == "CTDG")
        sim.CTDG(qubits[0], qubits[1]);

    else if (gate_name == "CSX")
        sim.CSX(params[0], qubits[0], qubits[1]);
    else if (gate_name == "CP")
        sim.CP(params[0], qubits[0], qubits[1]);

    else if (gate_name == "CSWAP")
    {
        sim.CX(qubits[2], qubits[1]);
        custom_ccx(qubits, params, sim);
        sim.CX(qubits[2], qubits[1]);
    }

    else if (gate_name == "ID")
        sim.ID(qubits[0]);
    else if (gate_name == "RZZ")
    {
        sim.CX(qubits[0], qubits[1]);
        sim.RZ(params[0], qubits[1]);
        sim.CX(qubits[0], qubits[1]);
    }
    else if (gate_name == "RXX")
    {
        sim.H(qubits[0]);
        sim.H(qubits[1]);
        sim.CX(qubits[0], qubits[1]);
        sim.RZ(params[0], qubits[1]);
        sim.CX(qubits[0], qubits[1]);
        sim.H(qubits[0]);
        sim.H(qubits[1]);
    }
    else if (gate_name == "RYY")
    {
        sim.RX(pi / 2, qubits[0]);
        sim.RX(pi / 2, qubits[1]);
        sim.CX(qubits[0], qubits[1]);
        sim.RZ(params[0], qubits[1]);
        sim.CX(qubits[0], qubits[1]);
        sim.RX(pi / 2, qubits[0]);
        sim.RX(pi / 2, qubits[1]);
    }

    // sim.reload_file();
}

/**
 * @brief Add gate to the circuit
 *
 * @param sim Simulator handle
 */
void apply_buffer(Simulation &sim)
{
    // Apply the buffered measurements if any
    for (string buffered_op : buffered_measure)
        load_gate(buffered_op, sim, false);
    buffered_measure.clear();
}

/**
 * @brief Add gate to the circuit
 *
 * @param op_prefix Gate name and parameters
 * @param qubit Qubit registers it apply to
 * @param sim Simulator handle
 */
void add_to_circuit(vector<string> op_prefix, vector<IdxType> qubits, Simulation &sim)
{
    apply_buffer(sim);

    string gate_name = op_prefix[0];
    vector<ValType> params;

    if (op_prefix.size() > 1)
    {
        op_prefix[1].pop_back(); // remove )
        vector<string> params_str = split(op_prefix[1], ',');
        for (string s : params_str)
        {
            params.push_back(get_param_value(s));
        }
    }

    append_gate(gate_name, qubits, params, sim);
}

/**
 * @brief Load the predefined gates from QASM file
 *
 * @param op_prefix Prefix that contains gate's name and parameters
 * @param qubit The qubit registers the gate is applying to
 * @param list_qregs Previouly loaded list of registers (both qreg and creg)
 * @return string Converted Representation of the Gate
 */
void load_default_gate(vector<string> op_prefix, string qubit, Simulation &sim)
{
    qubit.erase(remove(qubit.begin(), qubit.end(), ';'), qubit.end()); // remove ';'

    vector<vector<string>> regs;

    vector<string> given_regs = split(qubit, ',');

    int op_count = 0;

    for (int i = 0; i < given_regs.size(); i++)
    {
        vector<string> reg_strs;

        vector<string> reg_name_fields = split(given_regs[i], '[');

        reg cur_reg = list_qregs.at(reg_name_fields[0]);

        if (reg_name_fields.size() > 1)
        {
            int reg_index = stoi(split(reg_name_fields[1], ']')[0]);

            // reg_strs.push_back(given_regs[i]);
            reg_strs.push_back(to_string(cur_reg.offset + reg_index));
        }
        else
        {
            for (int j = 0; j < cur_reg.width; j++)
            {
                reg_strs.push_back(to_string(cur_reg.offset + j));
            }
        }
        op_count = reg_strs.size() > op_count ? reg_strs.size() : op_count;

        regs.push_back(reg_strs);
    }

    for (int i = 0; i < op_count; i++)
    {
        stringstream target_qubits;

        vector<IdxType> target_qubit_indices;

        for (int j = 0; j < regs.size(); j++)
        {
            if (regs[j].size() == 1)
            {
                target_qubits << regs[j][0];

                target_qubit_indices.push_back(stoi(regs[j][0]));
            }
            else
            {
                target_qubits << regs[j][i];

                target_qubit_indices.push_back(stoi(regs[j][i]));
            }
            target_qubits << " ";
        }

        add_to_circuit(op_prefix, target_qubit_indices, sim);
    }
}

/**
 * @brief Load custom gates
 *
 * @param op_prefix Gate prefix that includes gate's name and gate's parameters
 * @param args List of arguments for the gate
 * @return string Converted Representation of the Gate
 */
void load_custom_gate(vector<string> op_prefix, string args, Simulation &sim)
{

    for (gate gate : list_custom_gates)
    {
        if (gate.name == op_prefix[0])
        {
            vector<string> param_vals;
            vector<string> arg_vals;

            if (op_prefix.size() > 1)
            {
                op_prefix[1].pop_back(); // remove )
                param_vals = split(op_prefix[1], ',');
            }

            args.erase(remove(args.begin(), args.end(), ';'), args.end()); // remove ';'
            arg_vals = split(args, ',');

            for (string s : gate.ops)
            {
                stringstream new_op;

                vector<string> opParam_args = split(s, ' ');
                vector<string> op_params = split(opParam_args[0], '(');

                new_op << op_params[0];

                // update parameter values
                if (op_params.size() > 1)
                {
                    new_op << "(";

                    op_params[1].erase(remove(op_params[1].begin(), op_params[1].end(), ')'), op_params[1].end()); // remove ')'

                    vector<string> origional_params = split(op_params[1], ',');

                    bool append_comma = false;

                    for (string param : origional_params)
                    {
                        if (append_comma)
                            new_op << ",";

                        new_op << param_vals[get_index(gate.params, param)];
                        append_comma = true;
                    }
                    new_op << ")";
                }

                new_op << " ";

                vector<string> origional_args = split(opParam_args[1], ',');

                bool append_comma = false;
                for (string arg : origional_args)
                {
                    if (append_comma)
                        new_op << ",";

                    new_op << arg_vals[get_index(gate.args, arg)];
                    append_comma = true;
                }

                string test = new_op.str();

                load_gate(new_op.str(), sim, true);
            }
            return;
        }
    }
    // printf("%s\n", op_prefix[0].c_str());
}

void execute_if(string condition, string conditional_op, Simulation &sim)
{
    apply_buffer(sim);

    condition.pop_back(); // remove )
    vector<string> condition_terms = split(condition, '=');

    vector<string> reg_name_fields = split(condition_terms[0], '[');
    reg cur_reg = list_qregs.at(reg_name_fields[0]);

    int reg_index;
    if (reg_name_fields.size() > 1)
        reg_index = stoi(split(reg_name_fields[1], ']')[0]);
    else
        reg_index = 0;

    IdxType measured_results = sim.measure(cur_reg.offset + reg_index);

    sim.clear_circuit();

    if (measured_results == stoi(condition_terms[1]))
        load_gate(conditional_op, sim, true);
}

/**
 * @brief Load individual line of code in qasm file
 *
 * @param op QASM Code
 * @return string Converted Representation of the Gate
 */
void load_gate(string op, Simulation &sim, bool buffer_measure)
{
    vector<string> op_fields = split(op, ' ');

    if (op_fields.size() > 1)
    {
        transform(op_fields[0].begin(), op_fields[0].end(), op_fields[0].begin(), ::tolower);

        vector<string> op_prefix = split(op_fields[0], '(');

        if (op_prefix[0] == "if")
        {
            stringstream conditional_op;
            for (int i = 1; i < op_fields.size(); i++)
            {
                if (i != 1)
                    conditional_op << " ";
                conditional_op << op_fields[i];
            }
            execute_if(op_prefix[1], conditional_op.str(), sim);
        }
        else if (find(begin(DEFAULT_GATES), end(DEFAULT_GATES), op_prefix[0]) != end(DEFAULT_GATES))
        {
            load_default_gate(op_prefix, op_fields[1], sim);
        }
        else if (op_prefix[0] == "measure")
        {
            if (buffer_measure)
            {
                buffered_measure.push_back(op);
            }
            else
            {
                op_fields[3].erase(remove(op_fields[3].begin(), op_fields[3].end(), ';'), op_fields[3].end()); // remove '\r'
                vector<string> reg_name_fields = split(op_fields[1], '[');
                reg cur_reg = list_qregs.at(reg_name_fields[0]);

                int reg_index;
                if (reg_name_fields.size() > 1)
                    reg_index = stoi(split(reg_name_fields[1], ']')[0]);
                else
                    reg_index = 0;

                sim.M(cur_reg.offset + reg_index);
            }
        }
        else
        {
            load_custom_gate(op_prefix, op_fields[1], sim);
        }
    }
    // else
    // {
    //     printf("%s\n", op.c_str());
    // }
}

int get_global_qreg_index(string qreg_str)
{
    vector<string> qreg_fields = split(qreg_str, '[');

    int qreg_index = 0;
    if (qreg_fields.size() > 1)
    {
        qreg_index = stoi(split(qreg_fields[1], ']')[0]);
    }
    return list_qregs.at(qreg_fields[0]).offset + qreg_index;
}

string convert_outcome(IdxType original_out, map<string, vector<int>> creg_indicies)
{
    stringstream ss;

    int cur_index = 0;
    for (auto const &[key, val] : creg_indicies)
    {
        if (cur_index != 0)
            ss << " ";

        for (int i = val.size() - 1; i >= 0; i--)
        {
            int index = val[i];

            if (index == -1)
                ss << 0;
            else
                ss << ((original_out >> index) & 1);
        }

        cur_index++;
    }
    return ss.str();
}

map<string, int> *convert_dictionary(map<IdxType, int> dict)
{
    map<string, vector<int>> creg_indicies;

    for (reg reg : list_cregs)
    {
        vector<int> indicies;
        for (int i = 0; i < reg.width; i++)
            indicies.push_back(-1);

        creg_indicies.insert({reg.name, indicies});
    }

    vector<string> measured_cregs;
    for (string s : buffered_measure)
    {
        vector<string> op_fields = split(s, ' ');

        int qreg_offset = get_global_qreg_index(op_fields[1]);

        string creg = op_fields[3];

        creg.pop_back();

        vector<string> creg_fields = split(creg, '[');

        if (creg_fields.size() > 1)
        {
            int creg_bit_index = stoi(split(creg_fields[1], ']')[0]);

            creg_indicies.at(creg_fields[0])[creg_bit_index] = qreg_offset;
        }
        else
        {
            for (int i = 0; i < creg_indicies.at(creg_fields[0]).size(); i++)
            {
                creg_indicies.at(creg_fields[0])[i] = qreg_offset + i;
            }
        }
        measured_cregs.push_back(creg_fields[0]);
    }

    // Eliminate unmeasured cregs
    for (reg reg : list_cregs)
    {
        if (find(measured_cregs.begin(), measured_cregs.end(), reg.name) == measured_cregs.end())
            creg_indicies.erase(reg.name);
    }

    map<string, int> *converted_counts = new map<string, int>;
    for (auto const &[key, val] : dict)
    {
        string converted_key = convert_outcome(key, creg_indicies);
        if (converted_counts->count(converted_key))
        {
            converted_counts->at(converted_key) += val;
        }
        else
        {
            converted_counts->insert({converted_key, val});
        }
    }
    return converted_counts;

    // convert_outcome(IdxType original_out, map<string, vector<int>> creg_indicies)
}

map<string, int> *to_binary_dictionary(IdxType num_qubits, map<IdxType, int> counts)
{
    map<string, int> *binary_counts = new map<string, int>;

    stringstream ss;
    for (auto const &[key, val] : counts)
    {
        ss.str(string());
        for (int i = num_qubits - 1; i >= 0; i--)
        {
            ss << ((key >> i) & 1);
        }
        binary_counts->insert({ss.str(), val});
    }
    return binary_counts;
}

/**
 * @brief Parse and execute the given QASM file
 *
 * @param file The path to the qasm file.
 */
map<string, int> *execute_qasm(const char *file, Simulation &sim, IdxType repetition = DEFAULT_REPETITIONS)
{

    ifstream qasmFile(file);
    if (!qasmFile)
    {
        printf("%s %s\n", "Could not open ", file);
        return NULL;
    }

    list_qregs.clear();
    list_custom_gates.clear();
    buffered_measure.clear();

    string outFileName = "testing_regs.txt";

    ofstream outRegsFile(outFileName);

    outRegsFile << "QUBIT ALLOCATION\n";

    string line;
    int offset = 0;
    int loaded_regs = 0;
    int num_qubits = 0;

    while (getline(qasmFile, line))
    {
        if (line.size() > 1)
        {
            line.erase(remove(line.begin(), line.end(), '\r'), line.end()); // remove '\r'

            if (line.substr(0, 4) == "gate")
            {
                list_custom_gates.push_back(load_gate_defination(line, &qasmFile));
            }
            else if (line.substr(1, 3) == "reg")
            {
                reg temp_reg = load_reg(line, &offset);

                if (line[0] == 'q')
                {
                    list_qregs.insert({temp_reg.name, temp_reg});
                    num_qubits += temp_reg.width;
                }
                else
                    list_cregs.push_back(temp_reg);

                outRegsFile << temp_reg.name << "\n"
                            << "NUM: " << temp_reg.width << " OFFSET: " << temp_reg.offset << "\n";

                loaded_regs++;
            }
            else
            {
                load_gate(line, sim, true);
            }
        }
    }
    // MEASURE AND RETURN RESULTS;
    outRegsFile.close();

    IdxType *results = sim.measure_all(repetition);

    map<IdxType, int> result_dict;

    for (int i = 0; i < repetition; i++)
        if (result_dict.count(results[i]))
            result_dict[results[i]] = result_dict.at(results[i]) + 1;
        else
            result_dict[results[i]] = 1;

    if (buffered_measure.size() == 0 || list_cregs.size() == 0)

        return to_binary_dictionary(num_qubits, result_dict);
    else
        return convert_dictionary(result_dict);
}

/**
 * @brief Parse and get the number of qubits
 *
 * @param file The path to the qasm file.
 * */
int get_num_qubits(const char *file)
{
    ifstream qasmFile(file);
    if (!qasmFile)
    {
        printf("%s\n", "not opened!");
        return -1;
    }
    else
    {
        printf("Executing %s\n", file);
    }

    int num_qubits = 0;
    int offset = 0;
    string line;

    while (getline(qasmFile, line))
    {
        if (line.size() > 1)
        {
            line.erase(remove(line.begin(), line.end(), '\r'), line.end()); // remove '\r'

            if (line.substr(0, 4) == "qreg")
            {
                // printf("%s\n", line.c_str());
                reg temp_reg = load_reg(line, &offset);

                num_qubits += temp_reg.width;
            }
        }
    }

    qasmFile.close();

    return num_qubits;
}

#endif
