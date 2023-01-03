#ifndef QASM_PARSER
#define QASM_PARSER

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

/********** UPDATE THE INCLUDE PATH FOR LOCAL HEADER FILES HERE ************/
#include "../src/util.h"
#include "../src/svsim_nvgpu_mpi.cuh"
#include "parser_util.hpp"
/***************************************************************************/

using namespace std;
using namespace NWQSim;

const IdxType UN_DEF = -1;

string DEFAULT_GATES[] = {
    "U", "U3", "U2", "U1", "X", "Y", "Z", "H",
    "S", "SDG", "T", "TDG", "SX",
    "RX", "RY", "RZ",
    "CZ", "CX", "CY", "CH",
    "CCX", "CRX", "CRY", "CRZ", "CU1", "CU3",
    "RESET", "SWAP", "CSWAP",
    "ID", "RI", "P", "CS", "CSDG", "CT", "CTDG", "CSX", "CP",
    "RZZ", "RXX", "RYY"};

const string REG("REG");
const string GATE("GATE");
const string IF("IF");
const string MEASURE("MEASURE");

struct gate
{
    // Common gate fields
    string name;
    vector<ValType> params;
    vector<IdxType> qubits;

    // Fields used for measurement operations
    string creg_name;
    IdxType creg_index;
    bool final_measurements;

    // Fields used for if operations
    IdxType if_offset;
    IdxType if_creg_val;
};

struct custom_gate
{
    string name;
    vector<string> params;
    vector<string> args;

    vector<string> ops;
};

struct qreg
{
    string name;
    IdxType width;
    IdxType offset;
};

struct creg
{
    string name;
    IdxType width;
    vector<IdxType> qubit_indices;

    IdxType val = 0;
};

class qasm_parser
{
private:
    /* data */
    map<string, qreg> list_qregs;
    map<string, creg> list_cregs;

    vector<custom_gate> list_defined_gates;
    vector<gate> *list_gates = NULL;
    vector<gate> *list_conditional_gates = NULL;
    vector<gate> list_buffered_measure;

    IdxType _n_qubits = 0;
    IdxType global_qubit_offset = 0;

    string fileline;

    string outRegsFileName = "temp_files/regs_class.txt";
    ofstream outRegsFile;

    string outGateFileName = "temp_files/gate_class.txt";
    ofstream outGateFile;

    IdxType prior_if_index = UN_DEF;

    bool contains_IF = false;

    /* Local Helper Functions */
    void parse_reg();
    void parse_gate_defination(ifstream *file);
    void parse_gate(string gateline);

    void parse_native_gate(vector<string> op_prefix, string qubit);
    void parse_custom_gate(vector<string> op_prefix, string args);

    void classify_measurements();

    IdxType *sub_execute(Simulation &sim, IdxType repetition = DEFAULT_REPETITIONS);
    map<string, int> *convert_dictionary(map<IdxType, int> dict);
    string convert_outcome(IdxType original_out);
    map<string, int> *to_binary_dictionary(IdxType num_qubits, map<IdxType, int> counts);

    void print_cregs();

public:
    qasm_parser(const char *file);
    IdxType num_qubits();
    map<string, int> *execute(Simulation &sim, IdxType repetition = DEFAULT_REPETITIONS, bool measure_all = false, bool repeat_per_shot = false);

    ~qasm_parser();
};

string replace_pi(string s)
{
    vector<int> start_indices;
    vector<int> end_indices;
    for (int i = 0; i < s.size() - 1;)
    {
        if (s.substr(i, 2) == "PI")
        {
            bool nested_braket = false;
            bool visited_comma = false;

            int start_i, end_i;

            if (i > 0 && s[i - 1] == '(')
            {
                start_i = i - 1;
                for (int j = i; j < s.size(); j++)
                    if (s[j] == ',')
                        visited_comma = true;
                    else if (s[j] == ')' && !visited_comma)
                    {
                        end_i = j;
                        nested_braket = true;
                        break;
                    }
            }
            else if (i < s.size() - 2 && s[i + 2] == ')')
            {
                end_i = i + 2;
                for (int j = i; j > 0; j--)
                    if (s[j] == ',')
                        visited_comma = true;
                    else if (s[j] == '(' && !visited_comma)
                    {
                        start_i = j;
                        nested_braket = true;
                        break;
                    }
            }
            if (nested_braket)
            {
                start_indices.push_back(start_i);
                end_indices.push_back(end_i);
            }

            i += 2;
        }
        else
            i++;
    }

    for (int i = start_indices.size() - 1; i >= 0; i--)
    {
        int start_i = start_indices[i];
        int end_i = end_indices[i];

        if (s[end_i + 1] == ' ')
            continue;

        // cout << "Origianl " << s << endl;
        s.replace(start_i, end_i - start_i + 1, to_string(get_param_value(s.substr(start_i + 1, end_i - start_i - 1))));
        // cout << "Updated " << s << endl;
    }
    return s;
}

string cleanup_str(string s)
{
    // s.erase(remove_if(s.begin(), s.end(), invalid_char), s.end());

    // Remove '\r' from string
    s.erase(remove(s.begin(), s.end(), '\r'), s.end());

    // Remove leading and tailing spaces
    s = regex_replace(s, std::regex("^ +| +$|( ) +"), "$1");

    transform(s.begin(), s.end(), s.begin(), ::toupper);

    return replace_pi(s);
}

qasm_parser::qasm_parser(const char *file)
{
    list_gates = new vector<gate>;
    ifstream qasmFile(file);
    if (!qasmFile)
    {
        printf("Can not open %s\n", file);
        exit(-1);
    }

    outRegsFile = ofstream(outRegsFileName);
    outGateFile.open(outGateFileName);

    outRegsFile << "QUBIT ALLOCATION\n";

    IdxType index = 0;

    while (getline(qasmFile, fileline))
    {
        fileline = cleanup_str(fileline);

        if (fileline.size() > 1)
        {
            if (fileline.substr(0, 4) == GATE)
            {
                parse_gate_defination(&qasmFile);
            }
            else
            {
                fileline.erase(remove(fileline.begin(), fileline.end(), ';'), fileline.end());

                if (fileline.substr(1, 3) == REG)
                {
                    parse_reg();
                }
                else
                {
                    parse_gate(fileline);
                }
            }

            if (prior_if_index != UN_DEF)
            {
                list_gates->at(prior_if_index).if_offset = list_gates->size() - prior_if_index;
                prior_if_index = UN_DEF;
            }
        }

        index++;
    }

    qasmFile.close();
    outRegsFile.close();

    classify_measurements();
}

void qasm_parser::parse_reg()
{

    vector<string> reg_fields = split(split(fileline, ' ')[1], '[');

    string reg_name = reg_fields[0];
    string reg_witdth_str = regex_replace(reg_fields[1], std::regex(" +$"), "");
    reg_witdth_str.pop_back(); // remove ]
    IdxType reg_width = stoi(reg_witdth_str);

    if (fileline.at(0) == 'Q')
    {
        qreg qreg;

        qreg.name = reg_name;
        qreg.width = reg_width;
        qreg.offset = global_qubit_offset;

        global_qubit_offset += qreg.width;
        _n_qubits += qreg.width;

        list_qregs.insert({qreg.name, qreg});

        outRegsFile << qreg.name << "\n"
                    << "NUM: " << qreg.width << " OFFSET: " << qreg.offset << "\n";
    }
    else
    {
        creg creg;
        creg.name = reg_name;
        creg.width = reg_width;

        creg.qubit_indices.insert(creg.qubit_indices.end(), creg.width, UN_DEF);

        list_cregs.insert({creg.name, creg});

        outRegsFile << creg.name << "\n"
                    << "Classical NUM: " << creg.width << "\n";
    }
}

void qasm_parser::parse_gate_defination(ifstream *file)
{
    // Initialize gate object
    custom_gate cur_gate;

    stringstream ss;

    ss << fileline;

    if (fileline.find('}') == string::npos)
    {
        while (getline(*file, fileline))
        {
            fileline = cleanup_str(fileline);

            ss << fileline;

            if (fileline.find('}') != string::npos)
            {
                break;
            }
        }
    }
    string gate_def = ss.str();

    vector<string> gate_fields = split(gate_def, '{');

    gate_fields[1].erase(remove(gate_fields[1].begin(), gate_fields[1].end(), '}'), gate_fields[1].end()); // remove }
    vector<string> gate_ops = split(gate_fields[1], ';');

    // Spilit gate defination into name, parameter, and args
    vector<string> header_fields = split(gate_fields[0], ' ');
    vector<string> name_paramters = split(header_fields[1], '(');
    cur_gate.name = name_paramters[0];

    if (name_paramters.size() > 1)
    {
        name_paramters[1].pop_back();
        cur_gate.params = split(name_paramters[1], ',');
    }

    cur_gate.args = split(header_fields[2], ',');

    for (string op : gate_ops)
    {
        op = std::regex_replace(op, std::regex("^ +"), "");
        if (op.size() > 0)
        {
            //cout << op << endl;
            cur_gate.ops.push_back(op);
        }
    }

    list_defined_gates.push_back(cur_gate);

    outGateFile << "CUSTOME GATE DEF:\n";
    outGateFile << cur_gate.name << "\nPARAMS:\n";
    for (string s : cur_gate.params)

        outGateFile << s << " ";

    outGateFile << "\n\nArgs:\n";
    for (string s : cur_gate.args)
        outGateFile << s << " ";
    outGateFile << "\n\nOps:\n";
    for (string s : cur_gate.ops)
        outGateFile << s << endl;
    outGateFile << "-------------------\n\n";
}

void qasm_parser::parse_gate(string gateline)
{
    vector<string> op_fields = split(gateline, ' ');

    if (op_fields.size() > 1)
    {
        vector<string> op_prefix = split(op_fields[0], '(');

        if (find(begin(DEFAULT_GATES), end(DEFAULT_GATES), op_prefix[0]) != end(DEFAULT_GATES))
        {
            // Native Gate
            parse_native_gate(op_prefix, op_fields[1]);
        }
        else if (op_prefix[0] == IF)
        {
            contains_IF = true;

            op_prefix[1].pop_back(); // remove )

            gate gate;
            gate.name = IF;

            vector<string> creg_condition = split(op_prefix[1], '=');

            gate.creg_name = creg_condition[0];
            gate.if_creg_val = stoi(creg_condition[2]);

            prior_if_index = list_gates->size();
            list_gates->push_back(gate);

            stringstream ss;
            for (IdxType i = 1; i < op_fields.size(); i++)
                ss << op_fields[i] << " ";

            parse_gate(ss.str());
        }
        else if (op_prefix[0] == MEASURE)
        {
            vector<string> qreg_name_fields = split(op_fields[1], '[');
            vector<string> creg_name_fields = split(op_fields[3], '[');

            qreg cur_qreg = list_qregs.at(qreg_name_fields[0]);

            if (qreg_name_fields.size() > 1)
            // indexed measurements
            {
                gate gate;
                gate.name = MEASURE;

                gate.qubits.push_back(cur_qreg.offset + stoi(split(qreg_name_fields[1], ']')[0]));

                gate.creg_name = creg_name_fields[0];
                gate.creg_index = stoi(split(creg_name_fields[1], ']')[0]);

                list_gates->push_back(gate);
            }
            else
            // unindexed measurements
            {
                for (IdxType i = 0; i < cur_qreg.width; i++)
                {

                    gate gate;
                    gate.name = MEASURE;

                    gate.qubits.push_back(cur_qreg.offset + i);

                    gate.creg_name = creg_name_fields[0];
                    gate.creg_index = i;

                    list_gates->push_back(gate);
                }
            }
        }
        else
        {
            parse_custom_gate(op_prefix, op_fields[1]);
        }
    }
}

void qasm_parser::parse_native_gate(vector<string> op_prefix, string qubit)
{
    vector<vector<string>> regs;

    vector<string> given_regs = split(qubit, ',');

    int op_count = 0;

    for (int i = 0; i < given_regs.size(); i++)
    {
        vector<string> reg_strs;

        vector<string> reg_name_fields = split(given_regs[i], '[');

        qreg cur_qreg = list_qregs.at(reg_name_fields[0]);

        if (reg_name_fields.size() > 1)
        {
            IdxType reg_index = stoi(split(reg_name_fields[1], ']')[0]);

            reg_strs.push_back(to_string(cur_qreg.offset + reg_index));
        }
        else
        {
            for (int j = 0; j < cur_qreg.width; j++)
            {
                reg_strs.push_back(to_string(cur_qreg.offset + j));
            }
        }
        op_count = reg_strs.size() > op_count ? reg_strs.size() : op_count;

        regs.push_back(reg_strs);
    }

    for (int i = 0; i < op_count; i++)
    {

        gate gate;
        gate.name = op_prefix[0];

        for (int j = 0; j < regs.size(); j++)
            if (regs[j].size() == 1)
                gate.qubits.push_back(stoi(regs[j][0]));
            else
                gate.qubits.push_back(stoi(regs[j][i]));

        if (op_prefix.size() > 1)
        {
            op_prefix[1].pop_back(); // remove )
            vector<string> params_str = split(op_prefix[1], ',');
            for (string s : params_str)
            {
                gate.params.push_back(get_param_value(s));
            }
        }
        list_gates->push_back(gate);
    }
}

void qasm_parser::parse_custom_gate(vector<string> op_prefix, string args)
{
    for (custom_gate gate : list_defined_gates)
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

                        IdxType param_index = get_index(gate.params, param);

                        if (param_index == -1)
                            new_op << param;
                        else
                            new_op << param_vals[param_index];
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

                    IdxType arg_index = get_index(gate.args, arg);

                    if (arg_index == -1)
                        new_op << arg;
                    else
                        new_op << arg_vals[arg_index];

                    append_comma = true;
                }

                parse_gate(new_op.str());
            }
            return;
        }
    }
}

void qasm_parser::classify_measurements()
{
    bool final_measurements = true;

    for (IdxType i = list_gates->size() - 1; i >= 0; i--)
    {
        if (list_gates->at(i).name == MEASURE)
        {
            list_gates->at(i).final_measurements = final_measurements;
        }
        else
        {
            final_measurements = false;
        }
    }
}

void append_gate(Simulation &sim, string gate_name, vector<IdxType> qubits, vector<ValType> params)
{
    if (gate_name == "U")
        sim.U(params[0], params[1], params[2], qubits[0]);
    else if (gate_name == "U1")
        sim.U(0, 0, params[0], qubits[0]);
    else if (gate_name == "U2")
        sim.U(PI / 2, params[0], params[1], qubits[0]);
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
        sim.CCX(qubits[0], qubits[1], qubits[2]);
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
        sim.CSX(qubits[0], qubits[1]);
    else if (gate_name == "CP")
        sim.CP(params[0], qubits[0], qubits[1]);
    else if (gate_name == "CSWAP")
        sim.CSWAP(qubits[0], qubits[1], qubits[2]);
    else if (gate_name == "ID")
        sim.ID(qubits[0]);
    else if (gate_name == "RXX")
        sim.RXX(params[0],qubits[0], qubits[1]);
    else if (gate_name == "RYY")
        sim.RYY(params[0],qubits[0], qubits[1]);
    else if (gate_name == "RZZ")
        sim.RZZ(params[0],qubits[0], qubits[1]);
    else
        throw logic_error("Undefined gate is called!");
}

IdxType qasm_parser::num_qubits() { return _n_qubits; }

IdxType modifyBit(IdxType n, IdxType p, IdxType b)
{
    IdxType mask = 1 << p;
    return ((n & ~mask) | (b << p));
}

map<string, int> *qasm_parser::execute(Simulation &sim, IdxType repetition, bool measure_all, bool repeat_per_shot)
{
    IdxType *results;

    if (repeat_per_shot || contains_IF)
    {
        IdxType results_arr[repetition];

        results = results_arr;

        for (IdxType i = 0; i < repetition; i++)
        {
            IdxType *sub_result = sub_execute(sim, 1);
            results[i] = sub_result[0];
        }
    }
    else
    {
        results = sub_execute(sim, repetition);
    }

    map<IdxType, int> result_dict;

    for (int i = 0; i < repetition; i++)
        if (result_dict.count(results[i]))
            result_dict[results[i]] = result_dict.at(results[i]) + 1;
        else
            result_dict[results[i]] = 1;

    if (measure_all)
    {
        return to_binary_dictionary(_n_qubits, result_dict);
    }
    else
    {
        return convert_dictionary(result_dict);
    }
}

IdxType *qasm_parser::sub_execute(Simulation &sim, IdxType repetition)
{
    sim.clear_circuit();

    for (IdxType i = 0; i < list_gates->size();)
    {
        gate gate = list_gates->at(i);

        if (gate.name == IF)
        {
            creg creg = list_cregs.at(gate.creg_name);

            if (creg.val != gate.if_creg_val)
            {
                i += gate.if_offset;
            }
            else
            {
                i++;
            }
        }
        else
        {
            if (gate.name == MEASURE)
            {

                if (!gate.final_measurements)
                {
                    // Measure and update creg value for intermediate measurements
                    IdxType result = sim.measure(gate.qubits[0]);

                    list_cregs.at(gate.creg_name).val = modifyBit(list_cregs.at(gate.creg_name).val, gate.creg_index, result);
                    list_cregs.at(gate.creg_name).qubit_indices[gate.creg_index] = UN_DEF;
                }
                else
                {
                    // Update creg qubit indices for final measurements
                    list_cregs.at(gate.creg_name).qubit_indices[gate.creg_index] = gate.qubits[0];
                }
            }
            else
            {
                append_gate(sim, gate.name, gate.qubits, gate.params);
            }

            i++;
        }
    }

    return sim.measure_all(repetition);
}

map<string, int> *qasm_parser::convert_dictionary(map<IdxType, int> dict)
{
    map<string, int> *converted_counts = new map<string, int>;

    for (auto const &[key, val] : dict)
    {
        string converted_key = convert_outcome(key);

        if (converted_key.size() > 0)
        {
            if (converted_counts->count(converted_key))
                converted_counts->at(converted_key) += val;
            else
                converted_counts->insert({converted_key, val});
        }
    }
    return converted_counts;
}

void qasm_parser::print_cregs()
{
    for (auto const &[key, val] : list_cregs)
    {
        cout << key << " width: " << val.width;
        for (IdxType n : val.qubit_indices)
            cout << n << " ";
        cout << endl
             << endl;
    }
}

string qasm_parser::convert_outcome(IdxType original_out)
{
    stringstream ss;

    IdxType cur_index = 0;
    for (auto const &[key, val] : list_cregs)
    {
        if (cur_index != 0)
            ss << " ";

        vector<IdxType> creg_qubit_indices = val.qubit_indices;

        for (IdxType i = creg_qubit_indices.size() - 1; i >= 0; i--)
        {
            IdxType index = creg_qubit_indices[i];

            if (index == UN_DEF)
            {
                ss << 0;
            }
            else
                ss << ((original_out >> index) & 1);
        }

        cur_index++;
    }
    return ss.str();
}

map<string, int> *qasm_parser::to_binary_dictionary(IdxType num_qubits, map<IdxType, int> counts)
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

qasm_parser::~qasm_parser()
{
    if (list_gates != NULL)
    {
        delete list_gates;
    }
}

#endif
