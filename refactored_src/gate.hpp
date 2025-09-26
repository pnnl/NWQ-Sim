#pragma once

#include "nwqsim_core.hpp"

#include <cstddef>

namespace NWQSim
{
    enum class GateKind
    {
        RX,
        RZ,
        H,
        CX,
        Measure,
        MeasureAll,
        Reset
    };

    inline const char *gate_kind_name(GateKind kind)
    {
        switch (kind)
        {
        case GateKind::RX:
            return "RX";
        case GateKind::RZ:
            return "RZ";
        case GateKind::H:
            return "H";
        case GateKind::CX:
            return "CX";
        case GateKind::Measure:
            return "Measure";
        case GateKind::MeasureAll:
            return "MeasureAll";
        case GateKind::Reset:
            return "Reset";
        }
        return "Unknown";
    }

    struct Gate
    {
        GateKind kind = GateKind::H;
        IdxType target = -1;
        IdxType control = -1;
        IdxType repetition = 1;

        bool uses_param = false;
        size_t param_index = 0;
        ValType angle = 0.0;

        void set_angle(ValType value)
        {
            uses_param = false;
            angle = value;
        }

        void set_parameter_index(size_t index)
        {
            uses_param = true;
            param_index = index;
        }
    };
}

