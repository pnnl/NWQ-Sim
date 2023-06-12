#pragma once

#include "xacc.hpp"
#include <string>

namespace NWQSim
{
    // Forward declare
    class NWQBackend;
}

namespace xacc
{
    namespace quantum
    {

        class NWQAccelerator : public Accelerator
        {
        public:
            // Identifiable interface impls
            virtual const std::string name() const override { return "nwq-sim"; }
            virtual const std::string description() const override
            {
                return "XACC Simulation Accelerator based on NWQ-Sim library.";
            }

            // Accelerator interface impls
            virtual void initialize(const HeterogeneousMap &params = {}) override;
            virtual void updateConfiguration(const HeterogeneousMap &config) override
            {
                initialize(config);
            };
            virtual const std::vector<std::string> configurationKeys() override
            {
                return {};
            }
            virtual BitOrder getBitOrder() override { return BitOrder::LSB; }
            virtual void execute(std::shared_ptr<AcceleratorBuffer> buffer,
                                 const std::shared_ptr<CompositeInstruction>
                                     compositeInstruction) override;
            virtual void execute(std::shared_ptr<AcceleratorBuffer> buffer,
                                 const std::vector<std::shared_ptr<CompositeInstruction>>
                                     compositeInstructions) override;

        private:
            std::shared_ptr<NWQSim::NWQBackend> get_backend();

            int m_shots = 8192;
            std::string backend_name = "cpu";   //  cpu_omp
            std::string simulation_type = "sv"; // sv or dm
        };

    }
} // namespace xacc