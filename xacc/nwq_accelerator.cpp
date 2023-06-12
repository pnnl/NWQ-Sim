
#include <cassert>

#include "nwq_accelerator.hpp"
#include "nwq_api.hpp"

#include "AllGateVisitor.hpp"

namespace xacc
{
    namespace quantum
    {
        class NWQCircuitVisitor : public AllGateVisitor
        {
        private:
            NWQSim::NWQBackend *m_backend;

        public:
            using AllGateVisitor::visit;
            NWQCircuitVisitor(NWQSim::NWQBackend *backend) : m_backend(backend) {}

            void visit(Hadamard &h) override
            {
                m_backend->add_gate("H", std::vector<int>{(int)h.bits()[0]});
            }

            void visit(X &x) override
            {
                m_backend->add_gate("X", std::vector<int>{(int)x.bits()[0]});
            }

            void visit(Y &y) override
            {
                m_backend->add_gate("Y", std::vector<int>{(int)y.bits()[0]});
            }

            void visit(Z &z) override
            {
                m_backend->add_gate("Z", std::vector<int>{(int)z.bits()[0]});
            }

            void visit(Rx &rx) override
            {
                m_backend->add_gate("RX", std::vector<int>{(int)rx.bits()[0]},
                                    {InstructionParameterToDouble(rx.getParameter(0))});
            }

            void visit(Ry &ry) override
            {
                m_backend->add_gate("RY", std::vector<int>{(int)ry.bits()[0]},
                                    {InstructionParameterToDouble(ry.getParameter(0))});
            }

            void visit(Rz &rz) override
            {
                m_backend->add_gate("RZ", std::vector<int>{(int)rz.bits()[0]},
                                    {InstructionParameterToDouble(rz.getParameter(0))});
            }

            void visit(S &s) override
            {
                m_backend->add_gate("S", std::vector<int>{(int)s.bits()[0]});
            }

            void visit(Sdg &sdg) override
            {
                m_backend->add_gate("SDG", std::vector<int>{(int)sdg.bits()[0]});
            }

            void visit(T &t) override
            {
                m_backend->add_gate("T", std::vector<int>{(int)t.bits()[0]});
            }

            void visit(Tdg &tdg) override
            {
                m_backend->add_gate("TDG", std::vector<int>{(int)tdg.bits()[0]});
            }

            void visit(CNOT &cnot) override
            {
                m_backend->add_gate("CX", std::vector<int>{(int)cnot.bits()[0],
                                                           (int)cnot.bits()[1]});
            }

            void visit(CY &cy) override
            {
                m_backend->add_gate("CY",
                                    std::vector<int>{(int)cy.bits()[0], (int)cy.bits()[1]});
            }

            void visit(CZ &cz) override
            {
                m_backend->add_gate("CZ",
                                    std::vector<int>{(int)cz.bits()[0], (int)cz.bits()[1]});
            }

            void visit(Swap &s) override
            {
                m_backend->add_gate("SWAP",
                                    std::vector<int>{(int)s.bits()[0], (int)s.bits()[1]});
            }

            void visit(CH &ch) override
            {
                m_backend->add_gate("CH",
                                    std::vector<int>{(int)ch.bits()[0], (int)ch.bits()[1]});
            }

            void visit(CPhase &cphase) override
            {
                m_backend->add_gate(
                    "CP",
                    std::vector<int>{(int)cphase.bits()[0], (int)cphase.bits()[1]},
                    {InstructionParameterToDouble(cphase.getParameter(0))});
            }

            void visit(CRZ &crz) override
            {
                m_backend->add_gate("CRZ",
                                    std::vector<int>{(int)crz.bits()[0], (int)crz.bits()[1]},
                                    {InstructionParameterToDouble(crz.getParameter(0))});
            }

            void visit(Identity &i) override {}

            void visit(U &u) override
            {
                const auto theta = InstructionParameterToDouble(u.getParameter(0));
                const auto phi = InstructionParameterToDouble(u.getParameter(1));
                const auto lambda = InstructionParameterToDouble(u.getParameter(2));

                m_backend->add_gate("U", std::vector<int>{(int)u.bits()[0]},
                                    {theta, phi, lambda});
            }

            void visit(Measure &measure) override
            {
                m_measureQubits.emplace_back(measure.bits()[0]);
            }

            // NOT SUPPORTED:
            void visit(IfStmt &ifStmt) override {}
            std::vector<size_t> getMeasureBits() const { return m_measureQubits; }

        private:
            std::vector<size_t> m_measureQubits;
        };

        void
        NWQAccelerator::initialize(const HeterogeneousMap &params)
        {
            if (params.stringExists("sim-type"))
            {
                simulation_type = params.getString("sim-type");
            }

            if (params.stringExists("backend"))
            {
                backend_name = params.getString("backend");
            }

            if (params.keyExists<int>("shots"))
            {
                m_shots = params.get<int>("shots");
            }
        }

        void NWQAccelerator::execute(
            std::shared_ptr<AcceleratorBuffer> buffer,
            const std::shared_ptr<CompositeInstruction> circuit)
        {
            auto nwq_sim = NWQSim::get_backend();
            if (!nwq_sim)
            {
                throw std::logic_error("NWQ-Sim was not installed or selected backend is not supported.");
            }
            nwq_sim->init(buffer->size());

            // Create a visitor that will map IR to NWQ-Sim
            NWQCircuitVisitor visitor(nwq_sim.get());

            // Walk the IR tree, and visit each node
            InstructionIterator it(circuit);
            while (it.hasNext())
            {
                auto nextInst = it.next();
                if (nextInst->isEnabled())
                {
                    nextInst->accept(&visitor);
                }
            }

            auto measured_bits = visitor.getMeasureBits();
            if (measured_bits.empty())
            {
                // Default is just measure alls:
                for (size_t i = 0; i < buffer->size(); ++i)
                {
                    measured_bits.emplace_back(i);
                }
            }
            std::sort(measured_bits.begin(), measured_bits.end());

            buffer->addExtraInfo("exp-val-z", nwq_sim->calcExpectationValueZ(measured_bits));

            //  const auto measured_results = nwq_sim->measure(m_shots);
            //     const auto nwqSimMeasureToBitString = [&measured_bits](const auto &val)
            //     {
            //         std::string bitString;
            //         for (const auto &bit : measured_bits)
            //         {
            //             if (val & (1ULL << bit))
            //             {
            //                 bitString.push_back('1');
            //             }
            //             else
            //             {
            //                 bitString.push_back('0');
            //             }
            //         }
            //         std::reverse(bitString.begin(), bitString.end()); // Reverse the bit string

            //         return bitString;
            //     };

            //     for (const auto &m : measured_results)
            //     {
            //         buffer->appendMeasurement(nwqSimMeasureToBitString(m));
            //     }
        }

        void NWQAccelerator::execute(
            std::shared_ptr<AcceleratorBuffer> buffer,
            const std::vector<std::shared_ptr<CompositeInstruction>>
                circuits)
        {
            for (auto &f : circuits)
            {
                auto tmpBuffer =
                    std::make_shared<xacc::AcceleratorBuffer>(f->name(), buffer->size());
                execute(tmpBuffer, f);
                buffer->appendChild(f->name(), tmpBuffer);
            }
        }

    } // namespace quantum
} // namespace xacc

