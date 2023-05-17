#include <cassert>
#include <cmath>                  // For M_PI
#include "../include/circuit.hpp" // Assuming this is the file that contains your Circuit class

using namespace NWQSim;

bool compare_qasm(const std::string &generated, const std::string &expected)
{
    // This is a simple comparison. You might want to implement a more robust comparison that
    // ignores leading/trailing white space, differences in white space between tokens, etc.

    if (generated.compare(expected) == 0)
        return true;
    else
    {
        std::cout << generated << expected;
        return false;
    }
}

void run_tests()
{
    // X gate
    {
        Circuit circuit;
        circuit.X(0);

        assert(compare_qasm(circuit.circuitToString(), "X 0\n"));
    }

    // Y gate
    {
        Circuit circuit;
        circuit.Y(0);
        assert(compare_qasm(circuit.circuitToString(), "Y 0\n"));
    }

    // Z gate
    {
        Circuit circuit;
        circuit.Z(0);
        assert(compare_qasm(circuit.circuitToString(), "Z 0\n"));
    }

    // H gate
    {
        Circuit circuit;
        circuit.H(0);
        assert(compare_qasm(circuit.circuitToString(), "H 0\n"));
    }

    // S gate
    {
        Circuit circuit;
        circuit.S(0);
        assert(compare_qasm(circuit.circuitToString(), "S 0\n"));
    }

    // SDG gate
    {
        Circuit circuit;
        circuit.SDG(0);
        assert(compare_qasm(circuit.circuitToString(), "SDG 0\n"));
    }

    // T gate
    {
        Circuit circuit;
        circuit.T(0);
        assert(compare_qasm(circuit.circuitToString(), "T 0\n"));
    }

    // TDG gate
    {
        Circuit circuit;
        circuit.TDG(0);
        assert(compare_qasm(circuit.circuitToString(), "TDG 0\n"));
    }

    // RX gate
    {
        Circuit circuit;
        circuit.RX(M_PI / 2, 0);
        assert(compare_qasm(circuit.circuitToString(), "RX(1.5708) 0\n"));
    }

    // RY gate
    {
        Circuit circuit;
        circuit.RY(M_PI / 2, 0);
        assert(compare_qasm(circuit.circuitToString(), "RY(1.5708) 0\n"));
    }

    // RZ gate
    {
        Circuit circuit;
        circuit.RZ(M_PI / 2, 0);
        assert(compare_qasm(circuit.circuitToString(), "RZ(1.5708) 0\n"));
    }

    // U gate
    {
        Circuit circuit;
        circuit.U(M_PI / 2, M_PI / 2, M_PI / 2, 0);
        assert(compare_qasm(circuit.circuitToString(), "U(1.5708,1.5708,1.5708) 0\n"));
    }

    // CX gate
    {
        Circuit circuit;
        circuit.CX(0, 1);
        assert(compare_qasm(circuit.circuitToString(), "CX 0,1\n"));
    }

    // CY gate
    {
        Circuit circuit;
        circuit.CY(0, 1);
        assert(compare_qasm(circuit.circuitToString(), "CY 0,1\n"));
    }

    // CZ gate
    {
        Circuit circuit;
        circuit.CZ(0, 1);
        assert(compare_qasm(circuit.circuitToString(), "CZ 0,1\n"));
    }

    // CP gate
    {
        Circuit circuit;
        circuit.CP(M_PI / 2, 0, 1);
        assert(compare_qasm(circuit.circuitToString(), "CP(1.5708) 0,1\n"));
    }
    // CRX gate
    {
        Circuit circuit;
        circuit.CRX(M_PI / 2, 0, 1);
        assert(compare_qasm(circuit.circuitToString(), "CRX(1.5708) 0,1\n"));
    }

    // CRY gate
    {
        Circuit circuit;
        circuit.CRY(M_PI / 2, 0, 1);
        assert(compare_qasm(circuit.circuitToString(), "CRY(1.5708) 0,1\n"));
    }

    // CRZ gate
    {
        Circuit circuit;
        circuit.CRZ(M_PI / 2, 0, 1);
        assert(compare_qasm(circuit.circuitToString(), "CRZ(1.5708) 0,1\n"));
    }

    // SWAP gate
    {
        Circuit circuit;
        circuit.SWAP(0, 1);
        assert(compare_qasm(circuit.circuitToString(), "SWAP 0,1\n"));
    }

    // U3 gate
    {
        Circuit circuit;
        circuit.U3(M_PI / 2, M_PI / 2, M_PI / 2, 0);
        assert(compare_qasm(circuit.circuitToString(), "U(1.5708,1.5708,1.5708) 0\n"));
    }

    // U2 gate
    {
        Circuit circuit;
        circuit.U2(M_PI / 2, M_PI / 2, 0);
        assert(compare_qasm(circuit.circuitToString(), "U(1.5708,1.5708,1.5708) 0\n"));
    }

    // U1 gate
    {
        Circuit circuit;
        circuit.U1(M_PI / 2, 0);
        assert(compare_qasm(circuit.circuitToString(), "U(1.5708) 0\n"));
    }
}

int main()
{
    run_tests();
}