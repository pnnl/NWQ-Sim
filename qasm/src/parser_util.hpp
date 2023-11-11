#pragma once

#include <stdlib.h>
#include <string>
#include <vector>

#include <stack>
#include <queue>
#include <iostream>
#include <cmath>
#include <regex>
#include <cassert>

#include "lexer.hpp"
#include "qasm_parser_expr.hpp"
#include "nwq_util.hpp"

using namespace std;
using namespace NWQSim;
using namespace lexertk;

#define INST_NAME 0
#define INST_QASM_VERSION 1
#define INST_REG_NAME 1
#define INST_REG_WIDTH 3

#define INST_GATE_NAME 1

#define INST_MEASURE_QREG_NAME 1
#define INST_MEASURE_QREG_BIT 3
#define INST_MEASURE_CREG_NAME 7
#define INST_MEASURE_CREG_BIT 9

#define INST_IF_CREG 2
#define INST_IF_VAL 4
#define INST_IF_INST_START 6

typedef struct inst_indicies
{
    IdxType param_start = -1;
    IdxType param_end = -1;

    IdxType qubit_start = -1;
    IdxType qubit_end = -1;
} inst_indicies;

const IdxType UN_DEF = -1;

string DEFAULT_GATES[] = {
    "U", "U3", "U2", "U1", "X", "Y", "Z", "H",
    "S", "SDG", "T", "TDG", "SX",
    "RX", "RY", "RZ",
    "CZ", "CX", "CY", "CH",
    "CCX", "CRX", "CRY", "CRZ", "CU", "CU1", "CU3",
    "RESET", "SWAP", "CSWAP",
    "ID", "DELAY", "RI", "P", "CS", "CSDG", "CT", "CTDG", "CSX", "CP",
    "RZZ", "RXX", "RYY"};

const string OPENQASM("OPENQASM");
const string QREG("QREG");
const string CREG("CREG");
const string GATE("GATE");
const string IF("IF");
const string MEASURE("MEASURE");
const string BARRIER("BARRIER");

struct qasm_gate
{
    // Common qasm_gate fields
    string name;
    vector<ValType> params;
    vector<IdxType> qubits;
    IdxType repetation = 1;

    // Fields used for measurement operations
    string creg_name;
    IdxType creg_index;
    IdxType measured_qubit_index;
    bool final_measurements;

    // Fields used for if operations
    IdxType if_offset;
    IdxType if_creg_val;
    vector<qasm_gate> *conditional_inst;
};

struct defined_gate
{
    string name;
    vector<string> params;
    vector<string> qubits;
    vector<vector<token>> instructions;
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

IdxType get_last_rbraket(vector<token> &inst, IdxType start, IdxType end)
{
    for (IdxType i = end - 1; i > start; i--)
        if (inst[i].type == token::e_rbracket)
            return i;
    return -1;
}
inst_indicies get_indices(vector<token> &inst, IdxType start, IdxType end)
{
    inst_indicies indices;

    IdxType last_rbraket = get_last_rbraket(inst, start, end);
    if (last_rbraket == -1)
    {
        indices.qubit_start = start + 1;
        indices.qubit_end = end;
    }
    else
    {
        indices.param_start = start + 2;
        indices.param_end = last_rbraket;

        indices.qubit_start = last_rbraket + 1;
        indices.qubit_end = end;
    }
    return indices;
}

vector<ValType> get_params(vector<token> &inst, IdxType start, IdxType end)
{
    vector<ValType> params;
    if (start == -1)
        return params;

    IdxType cur_start = start;

    for (IdxType i = start; i < end; ++i)
    {
        if (inst[i].type == token::e_comma)
        {
            params.push_back(parse_expr(inst, cur_start, i));
            cur_start = i + 1;
        }
    }
    params.push_back(parse_expr(inst, cur_start, end));
    return params;
}

pair<IdxType, vector<vector<IdxType>>> get_qubits(vector<token> &inst, IdxType start, IdxType end, map<string, qreg> &list_qregs)
{
    vector<vector<IdxType>> qubits;
    IdxType repetition = 1;

    if (start == -1 || end == -1)
        throw runtime_error("No Qubits Found");

    for (IdxType i = start; i < end;)
    {
        vector<IdxType> sub_qubits;

        if (inst[i].type == token::e_number)
        {
            sub_qubits.push_back(stoll(inst[i].value));
            i++;
        }
        else
        {
            IdxType sub_qubit_offset = list_qregs.at(inst[i].value).offset;

            if (inst[i + 1].type == token::e_lsqrbracket)
            {
                sub_qubits.push_back(sub_qubit_offset + stoll(inst[i + 2].value));

                i += 4;
            }
            else
            {
                for (IdxType j = 0; j < list_qregs.at(inst[i].value).width; j++)
                    sub_qubits.push_back(sub_qubit_offset + j);

                repetition = sub_qubits.size();
                i++;
            }
        }

        qubits.push_back(sub_qubits);

        if (inst[i].type == token::e_eof)
            break;
        else
            i++;
    }
    return make_pair(repetition, qubits);
    // return qubits;
}

template <typename T>
vector<T> slices(vector<T> const &v, int m, int n)
{
    auto first = v.cbegin() + m;
    auto last = v.cbegin() + n;

    std::vector<T> vec(first, last);
    return vec;
}

IdxType modifyBit(IdxType n, IdxType p, IdxType b)
{
    IdxType mask = 1 << p;
    return ((n & ~mask) | (b << p));
}

string convert_outcome(IdxType original_out, map<string, creg> &cregs)
{
    stringstream ss;

    IdxType cur_index = 0;
    for (const auto &entry : cregs)
    {
        auto key = entry.first;
        auto val = entry.second;
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

map<string, IdxType> *convert_dictionary(map<IdxType, IdxType> &counts, map<string, creg> &cregs)
{
    map<string, IdxType> *converted_counts = new map<string, IdxType>;

    for (const auto &entry : counts)
    {
        auto key = entry.first;
        auto val = entry.second;
        string converted_key = convert_outcome(key, cregs);

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

map<string, IdxType> *to_binary_dictionary(IdxType num_qubits, map<IdxType, IdxType> &counts)
{
    map<string, IdxType> *binary_counts = new map<string, IdxType>;

    stringstream ss;
    for (const auto &entry : counts)
    {
        auto key = entry.first;
        auto val = entry.second;
        ss.str(string());
        for (IdxType i = num_qubits - 1; i >= 0; i--)
        {
            ss << ((key >> i) & 1);
        }
        binary_counts->insert({ss.str(), val});
    }
    return binary_counts;
}

map<string, IdxType> *outcome_to_dict(IdxType *results, IdxType num_qubits, IdxType repetition)
{
    map<IdxType, IdxType> result_dict;

    for (IdxType i = 0; i < repetition; i++)
    {
        if (result_dict.find(results[i]) != result_dict.end())
            result_dict[results[i]] += 1;
        else
            result_dict.insert({results[i], 1});
    }
    return to_binary_dictionary(num_qubits, result_dict);
}

/**************UTILITY FUNCTION DEFINIATIONS********************/

/**
 * @brief Split the given string based on the provided delim
 *
 * @param s
 * @param delim
 * @return vector<string>
 */
vector<string> split(const string &s, char delim);

/**
 * @brief Get the index of the target in the vector
 *
 * @param vec The vector that contains the target object
 * @param target Target to look for
 * @return int The found index of the target in the list
 */
int get_index(vector<string> vec, string target);

/**
 * @brief Print measurement outcomes
 *
 * @param counts Pointer to the counts dictionary.
 * @param repetition Number of shots performed.
 */
void print_counts(map<string, IdxType> *counts, IdxType repetition);

char *getCmdOption(char **begin, char **end, const std::string &option);
bool cmdOptionExists(char **begin, char **end, const std::string &option);

/************************** IMPLEMENTATION OF UTILITY FUNCTIONS **************************/

vector<string> split(const string &s, char delim)
{
    vector<string> elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}

int get_index(vector<string> vec, string target)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        if (vec[i] == target)
        {
            return i;
        }
    }
    return -1;
}

void print_counts(map<string, IdxType> *counts, IdxType repetition)
{
    assert(counts != NULL);
    printf("\n===============  Measurement (tests=%lld) ================\n", repetition);

    for (const auto &entry : *counts)
    {
        printf("\"%s\" : %lld\n", entry.first.c_str(), entry.second);
        // printf("\"%s\" : %lf\n", key.c_str(), double(val)/double(repetition));
    }
}

char *getCmdOption(char **begin, char **end, const string &option)
{
    char **itr = find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char **begin, char **end, const string &option)
{
    return find(begin, end, option) != end;
}
