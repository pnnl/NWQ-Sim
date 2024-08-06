#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>

#include "../nwq_util.hpp"
#include <cstring>
#include <algorithm>
#include <map>



#include "../private/config.hpp"
#include "../gate.hpp"


#include "../private/sim_gate.hpp"
#include "../private/gate_factory/device_noise.hpp"



namespace NWQSim
{
    // A utility function to get the OP values for the input basis gates
    OP getGateOp(const std::string &gateName)
    {
        for (int i = 0; i < sizeof(OP_NAMES) / sizeof(OP_NAMES[0]); ++i)
        {
            if (gateName == OP_NAMES[i])
            {
                return static_cast<OP>(i);
            }
        }
        throw std::runtime_error("Unknown gate name: " + gateName);
    }
    std::map<OP, std::map<std::string, DMGate>> getBsDMGateSP()
    {
        Config::ENABLE_NOISE = true;
        Config::readConfigFile();
        std::map<OP, std::map<std::string, DMGate>> dm_qubit_gates;

        // TODO: Make these parameters can be defined by the users
        IdxType qubit = 0;
        // IdxType ctrl = -1;
        ValType theta = 0;
        // ValType gate_len = theta;

        std::vector<std::string> basisGates = Config::backend_config["basis_gates"];
        IdxType num_qubits = Config::backend_config["num_qubits"];
        for (auto getName : basisGates)
        {
            std::transform(getName.begin(), getName.end(), getName.begin(), ::toupper);
            OP gate_op = getGateOp(getName);
            Gate G(OP::ID, qubit);
            if(gate_op == OP::CX){
                std::vector<DMGate> dm_gates;
                std::vector<std::string> cx_coupling = Config::backend_config["cx_coupling"];
                for (const auto &coupling : cx_coupling)
                {

                    size_t pos = coupling.find('_');
                    IdxType qubit = std::stoi(coupling.substr(0, pos));
                    IdxType ctrl = std::stoi(coupling.substr(pos + 1));
                    std::string key = std::to_string(qubit) + "_" + std::to_string(ctrl);
                    dm_qubit_gates[gate_op].emplace(key, generateDMGate(gate_op, qubit, ctrl, 0));

                }
            }
            else if (gate_op == OP::RESET){
                dm_qubit_gates[gate_op].emplace("0",DMGate(OP::RESET,0, 0));  
            }
            else if (gate_op == OP::X || 
                    gate_op == OP::ID || 
                    gate_op == OP::DELAY || 
                    gate_op == OP::RX || 
                    gate_op == OP::SX || 
                        gate_op == OP::RZ)
            {
                std::vector<DMGate> dm_gates;
                for (IdxType qubit = 0; qubit < num_qubits; ++qubit)
                {
                std::string key = std::to_string(qubit);
                dm_qubit_gates[gate_op].emplace(key,generateDMGate(gate_op, qubit, -1, 0));                
                }
            }
            else
            {
                throw std::invalid_argument("Unsupported basis gate: " + getName);
            }
    
            // dm_gates.emplace(gate_op, generateDMGate(G.op_name, G.qubit, G.ctrl, G.theta));
            // dm_gates[gate_op] = generateDMGate(G.op_name, G.qubit, G.ctrl, G.theta);

        }
        return dm_qubit_gates;
    }

nlohmann::json DMGateToJson(const DMGate &gate)
{
    nlohmann::json j;
    //TODO: Add parameters
    // Dertermine 1bit gate or 2bit gate
    // std::string qubit_couple = std::to_string(gate.qubit) + "_" + std::to_string(gate.ctrl);
    IdxType dim = (gate.op_name == OP::C2) ? 4 : 16; 
    std::vector<std::vector<std::string>> matrix(dim, std::vector<std::string>(dim));
    for (IdxType i = 0; i < dim; ++i)
    {
        for (IdxType j = 0; j < dim; ++j)
        {
            std::complex<ValType> complex_element(gate.gm_real[i * dim + j], gate.gm_imag[i * dim + j]);
            matrix[i][j] = std::to_string(complex_element.real()) + "+" + std::to_string(complex_element.imag()) + "i";
        }
    }

    // j["dim"] = dim;
    j["superoperator"] = matrix;

    return j;
}

void saveMapAsJson(const std::map<OP, std::map<std::string, DMGate>> &dm_qubit_gates, const std::string &filename)
{
    nlohmann::json j;
    for (const auto &pair : dm_qubit_gates)
    {
        json gates_json;
        for (const auto &gate_pair : pair.second)
        {
            gates_json[gate_pair.first] = DMGateToJson(gate_pair.second);
        }
        j[OP_NAMES[pair.first]] = gates_json;
    }

    std::ofstream file(filename);
    if (file.is_open())
    {
        file << j.dump(4); // Pretty print with 4 spaces
        file.close();
    }
    else
    {
        throw std::runtime_error("Unable to open file for writing");
    }
}


}