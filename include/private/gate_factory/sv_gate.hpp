#pragma once

#include <functional>
#include <map>
#include <cmath>

#include "../../nwq_util.hpp"
#include "../../gate.hpp"

#include "../sim_gate.hpp"

namespace NWQSim
{
    class GateFactory
    {
    public:
        using GateMatrixFunction = std::function<OP(const Gate &, ValType *, ValType *)>;

        static GateFactory &getInstance()
        {
            static GateFactory getInstance;
            return getInstance;
        }

        void registerGate(OP op, GateMatrixFunction func)
        {
            gateMatrixFunctions[op] = func;
        }

        std::vector<SVGate> getSVGates(const std::vector<Gate> &gates)
        {
            std::vector<SVGate> sim_sv_gates;

            for (const auto &g : gates)
            {
                auto it = gateMatrixFunctions.find(g.op_name);

                SVGate sv_gate(g.op_name, g.qubit, g.ctrl, g.data);

                if (it != gateMatrixFunctions.end())
                {
                    sv_gate.op_name = it->second(g, sv_gate.gm_real, sv_gate.gm_imag);
                }
                else if (g.op_name == OP::MA)
                {
                    sv_gate.qubit = g.repetition;
                }
                else if (!(g.op_name == OP::M || g.op_name == OP::RESET || g.op_name == OP::EXPECT || g.op_name == OP::MOD_NOISE))
                {
                    safe_print("%s\n", g.op_name);
                    throw std::runtime_error("Invalid gate operation");
                }

                sim_sv_gates.push_back(sv_gate);
            }
            return sim_sv_gates;
        }

    private:
        std::map<OP, GateMatrixFunction> gateMatrixFunctions;

        GateFactory() = default;
        GateFactory(const GateFactory &) = delete;
        GateFactory &operator=(const GateFactory &) = delete;
    };

    inline void registerGates()
    {

        GateFactory::getInstance().registerGate(
            OP::X,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {0, 1, 1, 0};
                ValType imag[4] = {0, 0, 0, 0};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::Y,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {0, 0, 0, 0};
                ValType imag[4] = {0, -1, 1, 0};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::Z,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {1, 0, 0, -1};
                ValType imag[4] = {0, 0, 0, 0};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::H,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {S2I, S2I, S2I, -S2I};
                ValType imag[4] = {0, 0, 0, 0};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::S,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {1, 0, 0, 0};
                ValType imag[4] = {0, 0, 0, 1};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::SDG,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {1, 0, 0, 0};
                ValType imag[4] = {0, 0, 0, -1};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::T,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {1, 0, 0, S2I};
                ValType imag[4] = {0, 0, 0, S2I};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::TDG,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {1, 0, 0, S2I};
                ValType imag[4] = {0, 0, 0, -S2I};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::RI,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {cos(g.theta), 0, 0, cos(g.theta)};
                ValType imag[4] = {sin(g.theta), 0, 0, sin(g.theta)};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::RX,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {cos(HALF * g.theta), 0, 0, cos(HALF * g.theta)};
                ValType imag[4] = {0, -sin(HALF * g.theta), -sin(HALF * g.theta), 0};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::RY,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {cos(HALF * g.theta), -sin(HALF * g.theta), sin(HALF * g.theta), cos(HALF * g.theta)};
                ValType imag[4] = {0};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::RZ,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {cos(HALF * g.theta), 0, 0, cos(HALF * g.theta)};
                ValType imag[4] = {-sin(HALF * g.theta), 0, 0, sin(HALF * g.theta)};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::SX,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {HALF, HALF, HALF, HALF};
                ValType imag[4] = {HALF, -HALF, -HALF, HALF};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::P,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {1, 0, 0, cos(g.theta)};
                ValType imag[4] = {0, 0, 0, sin(g.theta)};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::U,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {cos(HALF * g.theta),
                                   -cos(g.lam) * sin(HALF * g.theta),
                                   cos(g.phi) * sin(HALF * g.theta),
                                   cos(g.phi + g.lam) * cos(HALF * g.theta)};
                ValType imag[4] = {0,
                                   -sin(g.lam) * sin(HALF * g.theta),
                                   sin(g.phi) * sin(HALF * g.theta),
                                   sin(g.lam + g.phi) * cos(HALF * g.theta)};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::CX,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, 0, 1,
                                    0, 0, 1, 0};
                ValType imag[16] = {0};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });
        GateFactory::getInstance().registerGate(
            OP::ECR,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {0, S2I, 0, 0,
                                    S2I, 0, 0, 0,
                                    0, 0, 0, S2I,
                                    0, 0, S2I, 0};
                ValType imag[16] = {0, 0, 0, S2I,
                                    0, 0, -S2I, 0,
                                    0, S2I, 0, 0,
                                    -S2I, 0, 0, 0};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });
        GateFactory::getInstance().registerGate(
            OP::CY,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, 0};
                ValType imag[16] = {0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, -1,
                                    0, 0, 1, 0};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::CZ,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, 1, 0,
                                    0, 0, 0, -1};
                ValType imag[16] = {0};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::CH,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, S2I, S2I,
                                    0, 0, S2I, -S2I};
                ValType imag[16] = {0};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::CS,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, 1, 0,
                                    0, 0, 0, 0};
                ValType imag[16] = {0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, 1};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::CSDG,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, 1, 0,
                                    0, 0, 0, 0};
                ValType imag[16] = {0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, -1};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::CT,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, 1, 0,
                                    0, 0, 0, S2I};
                ValType imag[16] = {0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, S2I};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::CTDG,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, 1, 0,
                                    0, 0, 0, S2I};
                ValType imag[16] = {0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, -S2I};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::CRX,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, cos(HALF * g.theta), 0,
                                    0, 0, 0, cos(HALF * g.theta)};
                ValType imag[16] = {0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, -sin(HALF * g.theta),
                                    0, 0, -sin(HALF * g.theta), 0};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::CRY,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, cos(HALF * g.theta), -sin(HALF * g.theta),
                                    0, 0, sin(HALF * g.theta), cos(HALF * g.theta)};
                ValType imag[16] = {0};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::CRZ,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, cos(HALF * g.theta), 0,
                                    0, 0, 0, cos(HALF * g.theta)};
                ValType imag[16] = {0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, -sin(HALF * g.theta), 0,
                                    0, 0, 0, sin(HALF * g.theta)};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::CSX,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, HALF, HALF,
                                    0, 0, HALF, HALF};
                ValType imag[16] = {0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, HALF, -HALF,
                                    0, 0, -HALF, HALF};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::CP,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, 1, 0,
                                    0, 0, 0, cos(g.theta)};
                ValType imag[16] = {0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, sin(g.theta)};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::CU,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, cos(g.gamma) * cos(HALF * g.theta), -cos(g.gamma + g.lam) * sin(HALF * g.theta),
                                    0, 0, cos(g.gamma + g.phi) * sin(HALF * g.theta), cos(g.gamma + g.phi + g.lam) * cos(HALF * g.theta)};
                ValType imag[16] = {0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, sin(g.gamma) * cos(HALF * g.theta), -sin(g.gamma + g.lam) * sin(HALF * g.theta),
                                    0, 0, sin(g.gamma + g.phi) * sin(HALF * g.theta), sin(g.gamma + g.phi + g.lam) * cos(HALF * g.theta)};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::RXX,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {cos(HALF * g.theta), 0, 0, 0,
                                    0, cos(HALF * g.theta), 0, 0,
                                    0, 0, cos(HALF * g.theta), 0,
                                    0, 0, 0, cos(HALF * g.theta)};

                ValType imag[16] = {0, 0, 0, -sin(HALF * g.theta),
                                    0, 0, -sin(HALF * g.theta), 0,
                                    0, -sin(HALF * g.theta), 0, 0,
                                    -sin(HALF * g.theta), 0, 0, 0};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::RYY,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {cos(HALF * g.theta), 0, 0, 0,
                                    0, cos(HALF * g.theta), 0, 0,
                                    0, 0, cos(HALF * g.theta), 0,
                                    0, 0, 0, cos(HALF * g.theta)};

                ValType imag[16] = {0, 0, 0, sin(HALF * g.theta),
                                    0, 0, -sin(HALF * g.theta), 0,
                                    0, -sin(HALF * g.theta), 0, 0,
                                    sin(HALF * g.theta), 0, 0, 0};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::RZZ,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {cos(-HALF * g.theta), 0, 0, 0,
                                    0, cos(HALF * g.theta), 0, 0,
                                    0, 0, cos(HALF * g.theta), 0,
                                    0, 0, 0, cos(-HALF * g.theta)};

                ValType imag[16] = {sin(-HALF * g.theta), 0, 0, 0,
                                    0, sin(HALF * g.theta), 0, 0,
                                    0, 0, sin(HALF * g.theta), 0,
                                    0, 0, 0, sin(-HALF * g.theta)};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });

        GateFactory::getInstance().registerGate(
            OP::SX,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {HALF, HALF, HALF, HALF};
                ValType imag[4] = {HALF, -HALF, -HALF, HALF};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::ID,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[4] = {1, 0, 0, 1};
                ValType imag[4] = {0, 0, 0, 0};
                memcpy(gm_real, real, 4 * sizeof(ValType));
                memcpy(gm_imag, imag, 4 * sizeof(ValType));
                return OP::C1;
            });

        GateFactory::getInstance().registerGate(
            OP::SWAP,
            [](const Gate &g, ValType *gm_real, ValType *gm_imag)
            {
                ValType real[16] = {1, 0, 0, 0,
                                    0, 0, 1, 0,
                                    0, 1, 0, 0,
                                    0, 0, 0, 1};
                ValType imag[16] = {0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, 0,
                                    0, 0, 0, 0};
                memcpy(gm_real, real, 16 * sizeof(ValType));
                memcpy(gm_imag, imag, 16 * sizeof(ValType));
                return OP::C2;
            });
    }

}