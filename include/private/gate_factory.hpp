#pragma once

#include <functional>
#include <map>

#include "../public/util.hpp"
#include "../public/gate.hpp"

#include "sim_gate.hpp"

namespace NWQSim
{
    class GateFactory
    {
    public:
        using GateMatrixFunction = std::function<SVGate(const Gate &, SIM_TYPE)>;

        static GateFactory &getInstance()
        {
            static GateFactory getInstance;
            return getInstance;
        }

        void registerGate(OP op, GateMatrixFunction func)
        {
            gateMatrixFunctions[op] = func;
        }

        SVGate generateGate(const Gate &g, SIM_TYPE sim_type = SIM_TYPE::SV)
        {
            auto it = gateMatrixFunctions.find(g.op_name);

            if (it != gateMatrixFunctions.end())
            {
                SVGate sim_gate = it->second(g, sim_type);

                return sim_gate;
            }
            else if (g.op_name == OP::M)
            {
                if (sim_type == SIM_TYPE::SV)
                    return SVGate(g.op_name, g.qubit, g.ctrl);
            }
            else if (g.op_name == OP::MA)
            {
                if (sim_type == SIM_TYPE::SV)
                    return SVGate(g.op_name, g.repetition, g.ctrl);
                else
                    throw std::runtime_error("DM GATE NOT IMPLEMENTED");
            }
            else if (g.op_name == OP::RESET)
            {
                if (sim_type == SIM_TYPE::SV)
                    return SVGate(g.op_name, g.qubit, g.ctrl);
                else
                    throw std::runtime_error("DM GATE NOT IMPLEMENTED");
            }
            else
            {
                throw std::runtime_error("Invalid gate operation");
            }
            return SVGate(OP::ID, 0, 0);
        }

        std::vector<SVGate> generateGates(const std::vector<Gate> &gates)
        {
            std::vector<SVGate> sim_gates;

            for (const auto &g : gates)
            {
                sim_gates.push_back(generateGate(g));
            }
            return sim_gates;
        }

    private:
        std::map<OP, GateMatrixFunction> gateMatrixFunctions;

        GateFactory() = default;
        GateFactory(const GateFactory &) = delete;
        GateFactory &operator=(const GateFactory &) = delete;
    };

    inline void registerGates()
    {

        GateFactory::getInstance().registerGate(OP::X, [](const Gate &g, SIM_TYPE)
                                                {
                                                  SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {0, 1, 1, 0};
                                                 ValType imag[4] = {0, 0, 0, 0};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::Y, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {0, 0, 0, 0};
                                                 ValType imag[4] = {0, -1, 1, 0};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::Z, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {1, 0, 0, -1};
                                                 ValType imag[4] = {0, 0, 0, 0};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::H, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {S2I, S2I, S2I, -S2I};
                                                 ValType imag[4] = {0, 0, 0, 0};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::S, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {1, 0, 0, 0};
                                                 ValType imag[4] = {0, 0, 0, 1};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });
        GateFactory::getInstance().registerGate(OP::SDG, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {1, 0, 0, 0};
                                                 ValType imag[4] = {0, 0, 0, -1};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::T, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {1, 0, 0, S2I};
                                                 ValType imag[4] = {0, 0, 0, S2I};
                                                  memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::TDG, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {1, 0, 0, S2I};
                                                 ValType imag[4] = {0, 0, 0, -S2I};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::RI, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {cos(g.theta), 0, 0, cos(g.theta)};
                                                 ValType imag[4] = {sin(g.theta), 0, 0, sin(g.theta)};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::RX, [](const Gate &g, SIM_TYPE)
                                                {
                                                    SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {cos(HALF * g.theta), 0, 0, cos(HALF * g.theta)};
                                                 ValType imag[4] = {0, -sin(HALF * g.theta), -sin(HALF * g.theta), 0};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::RY, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {cos(HALF * g.theta), -sin(HALF * g.theta), sin(HALF * g.theta), cos(HALF * g.theta)};
                                                 ValType imag[4] = {0};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::RZ, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {cos(HALF * g.theta), 0, 0, cos(HALF * g.theta)};
                                                 ValType imag[4] = {-sin(HALF * g.theta), 0, 0, sin(HALF * g.theta)};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::SX, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {HALF, HALF, HALF, HALF};
                                                 ValType imag[4] = {HALF, -HALF, -HALF, HALF};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::P, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {1, 0, 0, cos(g.theta)};
                                                 ValType imag[4] = {0, 0, 0, sin(g.theta)};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::U, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                 ValType real[4] = {cos(HALF * g.theta),
                                                                    -cos(g.lam) * sin(HALF * g.theta),
                                                                    cos(g.phi) * sin(HALF * g.theta),
                                                                    cos(g.phi + g.lam) * cos(HALF * g.theta)};
                                                 ValType imag[4] = {0,
                                                                    -sin(g.lam) * sin(HALF * g.theta),
                                                                    sin(g.phi) * sin(HALF * g.theta),
                                                                    sin(g.lam + g.phi) * cos(HALF * g.theta)};
                                                 memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CX, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 0, 1,
                                                                     0, 0, 1, 0};
                                                 ValType imag[16] = {0};
                                                 memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CY, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, 0};
                                                 ValType imag[16] = {0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, -1,
                                                                     0, 0, 1, 0};
                                                 memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CZ, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 1, 0,
                                                                     0, 0, 0, -1};
                                                 ValType imag[16] = {0};
                                                   memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CH, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, S2I, S2I,
                                                                     0, 0, S2I, -S2I};
                                                 ValType imag[16] = {0};
                                                 memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CS, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 1, 0,
                                                                     0, 0, 0, 0};
                                                 ValType imag[16] = {0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, 1};
                                                 memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CSDG, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 1, 0,
                                                                     0, 0, 0, 0};
                                                 ValType imag[16] = {0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, -1};
                                                 memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CT, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, 1, 0,
                                                                     0, 0, 0, S2I};
                                                 ValType imag[16] = {0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, S2I};
                                                 memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CTDG, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                ValType real[16] = {1, 0, 0, 0,
                                                0, 1, 0, 0,
                                                0, 0, 1, 0,
                                                0, 0, 0, S2I};
                                                ValType imag[16] = {0, 0, 0, 0,
                                                0, 0, 0, 0, 
                                                0, 0, 0, 0,
                                                0, 0, 0, -S2I};  
                                                memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CRX, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, cos(HALF * g.theta), 0,
                                                                     0, 0, 0, cos(HALF * g.theta)};
                                                 ValType imag[16] = {0, 0, 0, 0,
                                                                     0, 0, 0, 0,
                                                                     0, 0, 0, -sin(HALF * g.theta),
                                                                     0, 0, -sin(HALF * g.theta), 0};
                                                 memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CRY, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                 ValType real[16] = {1, 0, 0, 0,
                                                                     0, 1, 0, 0,
                                                                     0, 0, cos(HALF * g.theta), -sin(HALF * g.theta),
                                                                     0, 0, sin(HALF * g.theta), cos(HALF * g.theta)};
                                                 ValType imag[16] = {0};
                                                 memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                 memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CRZ, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                ValType real[16] = {1, 0, 0, 0,
                                                                    0, 1, 0, 0,
                                                                    0, 0, cos(HALF * g.theta), 0,
                                                                    0, 0, 0, cos(HALF * g.theta)};
                                                ValType imag[16] = {0, 0, 0, 0,
                                                                    0, 0, 0, 0,
                                                                    0, 0, -sin(HALF * g.theta), 0,
                                                                    0, 0, 0, sin(HALF * g.theta)};
                                                memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CSX, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                ValType real[16] = {1, 0, 0, 0,
                                                                    0, 1, 0, 0,
                                                                    0, 0, HALF, HALF,
                                                                    0, 0, HALF, HALF};
                                                ValType imag[16] = {0, 0, 0, 0,
                                                                    0, 0, 0, 0,
                                                                    0, 0, HALF, -HALF,
                                                                    0, 0, -HALF, HALF};
                                                memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CP, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                ValType real[16] = {1, 0, 0, 0,
                                                                    0, 1, 0, 0,
                                                                    0, 0, 1, 0,
                                                                    0, 0, 0, cos(g.theta)};
                                                ValType imag[16] = {0, 0, 0, 0,
                                                                    0, 0, 0, 0,
                                                                    0, 0, 0, 0,
                                                                    0, 0, 0, sin(g.theta)};
                                                memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::CU, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                ValType real[16] = {1, 0, 0, 0,
                                                                    0, 1, 0, 0,
                                                                    0, 0, cos(g.gamma) * cos(HALF * g.theta), -cos(g.gamma + g.lam) * sin(HALF * g.theta),
                                                                    0, 0, cos(g.gamma + g.phi) * sin(HALF * g.theta), cos(g.gamma + g.phi + g.lam) * cos(HALF * g.theta)};
                                                ValType imag[16] = {0, 0, 0, 0,
                                                                    0, 0, 0, 0,
                                                                    0, 0, sin(g.gamma) * cos(HALF * g.theta), -sin(g.gamma + g.lam) * sin(HALF * g.theta),
                                                                    0, 0, sin(g.gamma + g.phi) * sin(HALF * g.theta), sin(g.gamma + g.phi + g.lam) * cos(HALF * g.theta)};
                                                memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::RXX, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                ValType real[16] = {cos(HALF * g.theta), 0, 0, 0,
                                                                    0, cos(HALF * g.theta), 0, 0,
                                                                    0, 0, cos(HALF * g.theta), 0,
                                                                    0, 0, 0, cos(HALF * g.theta)};

                                                ValType imag[16] = {0, 0, 0, -sin(HALF * g.theta),
                                                                    0, 0, -sin(HALF * g.theta), 0,
                                                                    0, -sin(HALF * g.theta), 0, 0,
                                                                    -sin(HALF * g.theta), 0, 0, 0};
                                                memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::RYY, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                ValType real[16] = {cos(HALF * g.theta), 0, 0, 0,
                                                                    0, cos(HALF * g.theta), 0, 0,
                                                                    0, 0, cos(HALF * g.theta), 0,
                                                                    0, 0, 0, cos(HALF * g.theta)};

                                                ValType imag[16] = {0, 0, 0, sin(HALF * g.theta),
                                                                    0, 0, -sin(HALF * g.theta), 0,
                                                                    0, -sin(HALF * g.theta), 0, 0,
                                                                    sin(HALF * g.theta), 0, 0, 0};
                                                memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::RZZ, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                ValType real[16] = {cos(-HALF * g.theta), 0, 0, 0,
                                                                    0, cos(HALF * g.theta), 0, 0,
                                                                    0, 0, cos(HALF * g.theta), 0,
                                                                    0, 0, 0, cos(-HALF * g.theta)};

                                                ValType imag[16] = {sin(-HALF * g.theta), 0, 0, 0,
                                                                    0, sin(HALF * g.theta), 0, 0,
                                                                    0, 0, sin(HALF * g.theta), 0,
                                                                    0, 0, 0, sin(-HALF * g.theta)};
                                                memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::SX, [](const Gate &g, SIM_TYPE)
                                                {
                                                     SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                ValType real[4] = {HALF, HALF, HALF, HALF};
                                                ValType imag[4] = {HALF, -HALF, -HALF, HALF};
                                                memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::ID, [](const Gate &g, SIM_TYPE)
                                                {
                                                SVGate sim_gate(OP::C1, g.qubit, g.ctrl);  
                                                ValType real[4] = {1, 0, 0, 1};
                                                ValType imag[4] = {0, 0, 0, 0};
                                                memcpy(sim_gate.gm_real, real, 4 * sizeof(ValType));
                                                memcpy(sim_gate.gm_imag, imag, 4 * sizeof(ValType)); 
                                                 return sim_gate; });

        GateFactory::getInstance().registerGate(OP::SWAP, [](const Gate &g, SIM_TYPE)
                                                {
                                                SVGate sim_gate(OP::C2, g.qubit, g.ctrl);  
                                                ValType real[16] = {1, 0, 0, 0,
                                                                    0, 0, 1, 0,
                                                                    0, 1, 0, 0,
                                                                    0, 0, 0, 1};
                                                ValType imag[16] = {0, 0, 0, 0,
                                                                    0, 0, 0, 0,
                                                                    0, 0, 0, 0,
                                                                    0, 0, 0, 0};
                                                memcpy(sim_gate.gm_real, real, 16 * sizeof(ValType));
                                                memcpy(sim_gate.gm_imag, imag, 16 * sizeof(ValType)); 
                                                 return sim_gate; });
    }

}