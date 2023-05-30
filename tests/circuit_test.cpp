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
        Circuit circuit(1);
        circuit.X(0);

        assert(compare_qasm(circuit.to_string(), "X 0\n"));
    }

    // Y gate
    {
        Circuit circuit(1);
        circuit.Y(0);
        assert(compare_qasm(circuit.to_string(), "Y 0\n"));
    }

    // Z gate
    {
        Circuit circuit(1);
        circuit.Z(0);
        assert(compare_qasm(circuit.to_string(), "Z 0\n"));
    }

    // H gate
    {
        Circuit circuit(1);
        circuit.H(0);
        assert(compare_qasm(circuit.to_string(), "H 0\n"));
    }

    // S gate
    {
        Circuit circuit(1);
        circuit.S(0);
        assert(compare_qasm(circuit.to_string(), "S 0\n"));
    }

    // SDG gate
    {
        Circuit circuit(1);
        circuit.SDG(0);
        assert(compare_qasm(circuit.to_string(), "SDG 0\n"));
    }

    // T gate
    {
        Circuit circuit(1);
        circuit.T(0);
        assert(compare_qasm(circuit.to_string(), "T 0\n"));
    }

    // TDG gate
    {
        Circuit circuit(1);
        circuit.TDG(0);
        assert(compare_qasm(circuit.to_string(), "TDG 0\n"));
    }

    // RX gate
    {
        Circuit circuit(1);
        circuit.RX(M_PI / 2, 0);
        assert(compare_qasm(circuit.to_string(), "RX(1.5708) 0\n"));
    }

    // RY gate
    {
        Circuit circuit(1);
        circuit.RY(M_PI / 2, 0);
        assert(compare_qasm(circuit.to_string(), "RY(1.5708) 0\n"));
    }

    // RZ gate
    {
        Circuit circuit(1);
        circuit.RZ(M_PI / 2, 0);
        assert(compare_qasm(circuit.to_string(), "RZ(1.5708) 0\n"));
    }

    // U gate
    {
        Circuit circuit(1);
        circuit.U(M_PI / 2, M_PI / 2, M_PI / 2, 0);
        assert(compare_qasm(circuit.to_string(), "U(1.5708,1.5708,1.5708) 0\n"));
    }

    // CX gate
    {
        Circuit circuit(2);
        circuit.CX(0, 1);
        assert(compare_qasm(circuit.to_string(), "CX 0,1\n"));
    }

    // CY gate
    {
        Circuit circuit(2);
        circuit.CY(0, 1);
        assert(compare_qasm(circuit.to_string(), "CY 0,1\n"));
    }

    // CZ gate
    {
        Circuit circuit(2);
        circuit.CZ(0, 1);
        assert(compare_qasm(circuit.to_string(), "CZ 0,1\n"));
    }

    // CP gate
    {
        Circuit circuit(2);
        circuit.CP(M_PI / 2, 0, 1);
        assert(compare_qasm(circuit.to_string(), "CP(1.5708) 0,1\n"));
    }
    // CRX gate
    {
        Circuit circuit(2);
        circuit.CRX(M_PI / 2, 0, 1);
        assert(compare_qasm(circuit.to_string(), "CRX(1.5708) 0,1\n"));
    }

    // CRY gate
    {
        Circuit circuit(2);
        circuit.CRY(M_PI / 2, 0, 1);
        assert(compare_qasm(circuit.to_string(), "CRY(1.5708) 0,1\n"));
    }

    // CRZ gate
    {
        Circuit circuit(2);
        circuit.CRZ(M_PI / 2, 0, 1);
        assert(compare_qasm(circuit.to_string(), "CRZ(1.5708) 0,1\n"));
    }

    // SWAP gate
    {
        Circuit circuit(2);
        circuit.SWAP(0, 1);
        assert(compare_qasm(circuit.to_string(), "SWAP 0,1\n"));
    }

    // U3 gate
    {
        Circuit circuit(1);
        circuit.U3(M_PI / 2, M_PI / 2, M_PI / 2, 0);
        assert(compare_qasm(circuit.to_string(), "U(1.5708,1.5708,1.5708) 0\n"));
    }

    // U2 gate
    {
        Circuit circuit(1);
        circuit.U2(M_PI / 2, M_PI / 2, 0);
        assert(compare_qasm(circuit.to_string(), "U(1.5708,1.5708,1.5708) 0\n"));
    }

    // U1 gate
    {
        Circuit circuit(1);
        circuit.U1(M_PI / 2, 0);
        assert(compare_qasm(circuit.to_string(), "U(1.5708) 0\n"));
    }
}

int main()
{
    run_tests();
}