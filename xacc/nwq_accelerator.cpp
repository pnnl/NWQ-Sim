
#include <cassert>
#include <chrono>
#include <xacc.hpp>
#include "IRUtils.hpp"

#include "AllGateVisitor.hpp"

#include "nwq_accelerator.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/backendManager.hpp"
#include "../include/nwq_util.hpp"
namespace xacc
{
    namespace quantum
    {
        class NWQCircuitVisitor : public AllGateVisitor
        {

        private:
            std::shared_ptr<NWQSim::Circuit> m_circuit;

        public:
            using AllGateVisitor::visit;

            NWQCircuitVisitor(int n_qubits)
            {
                m_circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
            }

            void visit(Hadamard &h) override
            {
                m_circuit->H(h.bits()[0]);
            }

            void visit(X &x) override
            {
                m_circuit->X(x.bits()[0]);
            }

            void visit(Y &y) override
            {
                m_circuit->Y(y.bits()[0]);
            }

            void visit(Z &z) override
            {
                m_circuit->Z(z.bits()[0]);
            }

            void visit(Rx &rx) override
            {
                m_circuit->RX(InstructionParameterToDouble(rx.getParameter(0)), rx.bits()[0]);
            }

            void visit(Ry &ry) override
            {
                m_circuit->RY(InstructionParameterToDouble(ry.getParameter(0)), ry.bits()[0]);
            }

            void visit(Rz &rz) override
            {
                m_circuit->RZ(InstructionParameterToDouble(rz.getParameter(0)), rz.bits()[0]);
            }

            void visit(S &s) override
            {
                m_circuit->S(s.bits()[0]);
            }

            void visit(Sdg &sdg) override
            {
                m_circuit->SDG(sdg.bits()[0]);
            }

            void visit(T &t) override
            {
                m_circuit->T(t.bits()[0]);
            }

            void visit(Tdg &tdg) override
            {
                m_circuit->TDG(tdg.bits()[0]);
            }

            void visit(CNOT &cnot) override
            {
                m_circuit->CX(cnot.bits()[0], cnot.bits()[1]);
            }

            void visit(CY &cy) override
            {
                m_circuit->CY(cy.bits()[0], cy.bits()[1]);
            }

            void visit(CZ &cz) override
            {
                m_circuit->CZ(cz.bits()[0], cz.bits()[1]);
            }

            void visit(Swap &s) override
            {
                m_circuit->SWAP(s.bits()[0], s.bits()[1]);
            }

            void visit(CH &ch) override
            {
                m_circuit->CH(ch.bits()[0], ch.bits()[1]);
            }

            void visit(CPhase &cphase) override
            {
                m_circuit->CP(InstructionParameterToDouble(cphase.getParameter(0)), cphase.bits()[0], cphase.bits()[1]);
            }

            void visit(CRZ &crz) override
            {
                m_circuit->CRZ(InstructionParameterToDouble(crz.getParameter(0)), crz.bits()[0], crz.bits()[1]);
            }

            void visit(Identity &i) override {}

            void visit(U &u) override
            {
                const auto theta = InstructionParameterToDouble(u.getParameter(0));
                const auto phi = InstructionParameterToDouble(u.getParameter(1));
                const auto lambda = InstructionParameterToDouble(u.getParameter(2));

                m_circuit->U(theta, phi, lambda, u.bits()[0]);
            }

            void visit(Measure &measure) override
            {
                m_measureQubits.emplace_back(measure.bits()[0]);
            }

            // NOT SUPPORTED:
            void visit(IfStmt &ifStmt) override {}

            std::vector<size_t> getMeasureBits() const { return m_measureQubits; }

            std::shared_ptr<NWQSim::Circuit> getNWQCircuit() const { return m_circuit; }

        private:
            std::vector<size_t> m_measureQubits;
        };

        void NWQAccelerator::initialize(const HeterogeneousMap &params)
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

            if (params.keyExists<bool>("vqe_mode"))
            {
                vqe_mode = params.get<bool>("vqe_mode");
            }
        }

        std::unique_ptr<NWQCircuitVisitor> get_circuit_visitor(const int n_qubits, const std::shared_ptr<CompositeInstruction> circuit)
        {
            auto visitor = std::make_unique<NWQCircuitVisitor>(n_qubits);
            InstructionIterator it(circuit);
            while (it.hasNext())
            {
                auto nextInst = it.next();
                if (nextInst->isEnabled())
                {
                    nextInst->accept(visitor.get());
                }
            }

            return visitor;
        }

        void NWQAccelerator::execute(
            std::shared_ptr<AcceleratorBuffer> buffer,
            const std::shared_ptr<CompositeInstruction> circuit)
        {

            m_state = BackendManager::create_state(backend_name, buffer->size(), simulation_type);

            if (!m_state)
            {
                throw std::logic_error("NWQ-Sim was not installed or selected backend is not supported.");
            }

            std::unique_ptr<NWQCircuitVisitor> visitor = get_circuit_visitor(buffer->size(), circuit);

            m_state->sim(visitor->getNWQCircuit());

            auto measured_bits = visitor->getMeasureBits();

            if (m_shots == 0)
            {
                if (measured_bits.empty() || measured_bits.size() == buffer->size())
                {
                    // Default is just measure alls:
                    buffer->addExtraInfo("exp-val-z", m_state->get_exp_z());
                }
                else
                {
                    std::sort(measured_bits.begin(), measured_bits.end());
                    buffer->addExtraInfo("exp-val-z", m_state->get_exp_z(measured_bits));
                }
            }
            else
            {
                const auto measured_results = m_state->measure_all(m_shots);

                const auto nwqSimMeasureToBitString = [&measured_bits](const auto &val)
                {
                    std::string bitString;
                    for (const auto &bit : measured_bits)
                    {
                        if (val & (1ULL << bit))
                        {
                            bitString.push_back('1');
                        }
                        else
                        {
                            bitString.push_back('0');
                        }
                    }
                    std::reverse(bitString.begin(), bitString.end()); // Reverse the bit string

                    return bitString;
                };

                for (int i = 0; i < m_shots; i++)
                {
                    buffer->appendMeasurement(nwqSimMeasureToBitString(measured_results[i]));
                }
            }
        }

        void NWQAccelerator::execute(
            std::shared_ptr<AcceleratorBuffer> buffer,
            const std::vector<std::shared_ptr<CompositeInstruction>>
                circuits)
        {
            if (circuits.size() == 1 || !vqe_mode)
            {
                // Only one circuit, or vqe_mode is not enabled, just execute them one by one
                for (auto &f : circuits)
                {
                    auto tmpBuffer =
                        std::make_shared<xacc::AcceleratorBuffer>(f->name(), buffer->size());
                    execute(tmpBuffer, f);
                    buffer->appendChild(f->name(), tmpBuffer);
                }
            }
            else
            {
                // Multiple circuits, execute them in parallel
                auto kernelDecomposed = ObservedAnsatz::fromObservedComposites(circuits);
                // Always validate kernel decomposition in DEBUG
                assert(kernelDecomposed.validate(circuits));

                // Base kernel:
                auto baseKernel = kernelDecomposed.getBase();
                // Basis-change + measures
                auto obsCircuits = kernelDecomposed.getObservedSubCircuits();

                m_state = BackendManager::create_state(backend_name, buffer->size(), simulation_type);
                if (!m_state)
                {
                    throw std::logic_error("NWQ-Sim was not installed or selected backend is not supported.");
                }

                // parse ansatz circuit
                std::unique_ptr<NWQCircuitVisitor> visitor = get_circuit_visitor(buffer->size(), baseKernel);

                // simulate and store post-anstaz state
                m_state->sim(visitor->getNWQCircuit());
                m_state->save_state();

                for (auto &f : obsCircuits)
                {
                    // parse observable circuit
                    visitor = get_circuit_visitor(buffer->size(), f);

                    m_state->load_state();
                    m_state->sim(visitor->getNWQCircuit());
                    auto measured_bits = visitor->getMeasureBits();

                    auto tmpBuffer =
                        std::make_shared<xacc::AcceleratorBuffer>(f->name(), buffer->size());

                    if (measured_bits.empty() || measured_bits.size() == buffer->size())
                    {
                        // Default is just measure alls:
                        tmpBuffer->addExtraInfo("exp-val-z", m_state->get_exp_z());
                    }
                    else
                    {
                        std::sort(measured_bits.begin(), measured_bits.end());
                        tmpBuffer->addExtraInfo("exp-val-z", m_state->get_exp_z(measured_bits));
                    }

                    buffer->appendChild(f->name(), tmpBuffer);
                }
                m_state->clear_state();
            }

        }

    } // namespace quantum
} // namespace xacc
