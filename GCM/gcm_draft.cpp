#include "gcm_header.hpp"

int main() {

    std::vector<std::string> bases = {
        "R0",
        "+R1", "-R1",
        "+R2", "-R2",
        "+R3", "-R3",
        "+R4", "-R4",
        "+R3+R5", "+R4+R6",
        "+R2+R1", "+R2-R1", "-R2+R1", "-R2-R1"
    };
    std::cout << "Test";
    //Generate ts
    uint64_t seed = 1234; //Seed for the PCG64 generator
    Vector ts(7); //Generate random numbers w/ PCG64 and scale them to be between 0 and 100
    ts[0] = 1.0; //Inserting 1 for R0
    for (int i = 1; i < ts.size(); i++) {
        ts[i] = generateRandomNumber(seed);
    }

    Matrix Hamiltonian; // Hamiltonian matrix
    Matrix HFstate;     // Hartree-Fock state matrix
    Matrix Z;
  
    int numOrbitals = 8; // Define number of orbitals in the system
    int numParticles = 4;
   
    std::map<std::string, std::vector<Term>> basisDict = {
        {"R1", {{"+_2 -_1", 0.1}, {"+_6 -_5", -0.1}}}
    };

    


    // Step 1: Perform Jordan-Wigner transformation
    // Step 2: Generate unitaries
    // Step 3: Trotterize unitaries into Pauli strings
    // Step 4-8: Compute matrix elements classically
    int matLen = bases.size();
    int matSize = matLen * matLen;
    // Define matrices of complex numbers
    ComplexMatrix Sclass(matLen, matLen);
    ComplexMatrix Hclass(matLen, matLen);
    // Initialize matrices with zeros
    Sclass.setZero();
    Hclass.setZero();
    // Compute function fills in the S and H matrices
    ComputeSH(Sclass, Hclass, bases, ts, basisDict, numOrbitals, matLen, true);

    //Matrix SclassReal = Sclass.real();
    //Matrix HclassReal = Hclass.real();

    // Steps 9 & 10: Evaluate S and H quantumly
    //ComplexMatrix Squant(matLen, matLen);
    //ComplexMatrix Hquant(matLen, matLen);
    // Initialize matrices with zeros
    //Squant.setZero();
    //Hquant.setZero();
    //EvaluateSH(Squant, Hquant, bases, ts);

    // Step 13: Solve the general eigenvalue problem classically
    //Matrix overlapMatrix;
    //auto [eigenvalues, eigenvectors] = SolveEigenvalueProblem(evaluatedMatrixElements, overlapMatrix);

    // Output the results
    //std::cout << "Eigenvalues:\n" << eigenvalues << std::endl;
    //std::cout << "Eigenvectors:\n" << eigenvectors << std::endl;
    std::cout<<"Press Return to exit...";
    std::cin.sync();
    std::cin.ignore();
    return 0;
}