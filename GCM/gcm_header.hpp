#include <vector>
#include <string>
#include <map>
#include <complex>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <map>
#include <exception>
#include <utility>
#include <cmath>
#include "pcg_random.hpp"
#include "../../Eigen/Dense"

//Using dynamic Eigen MatrixXd to streamline calculations
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using ComplexMatrix2D = Eigen::Matrix2cd;
using ComplexMatrix = Eigen::MatrixXcd;
using ComplexVector = Eigen::VectorXcd;

typedef std::pair<std::string, double> Term; //A fermionic string <operator, coefficient>
typedef std::vector<Term> CircuitConfig; //Series of strings

int generateRandomNumber(uint64_t seed) {
    //Create a PCG random number generator initialized with the given seed
    pcg64 rng(seed);
    std::uniform_int_distribution<int> dist(0, 100); // Distribution that maps to the range [0,100]
    return dist(rng); // Generate and return the random number
}

//Decodes R bases into the a fermionic string
std::vector<Term> baseCircuitCodeSingle(const std::string& singleName, double t, 
                                        const std::map<std::string, std::vector<Term>>& basisDict, bool left) {
    std::vector<Term> totalConfig;
    char signStr = singleName[0]; // '+' or '-'
    double tAct = 0.0;

    // Determine the action based on the sign and whether it is on the left
    if (signStr == '+') {
        tAct = (left ? -t : t);
    } else if (signStr == '-') {
        tAct = (left ? t : -t);
    } else {
        throw std::invalid_argument("Unknown Sign");
    }

    // Extract the base string and retrieve the corresponding configuration
    std::string baseStr = singleName.substr(1); //'R<num>', e.g., 'R1'
    auto it = basisDict.find(baseStr);
    if (it == basisDict.end()) {
        throw std::invalid_argument("Basis configuration not found");
    }

    // Apply the time scaling factor to each term in the configuration
    for (const Term& term : it->second) {
        totalConfig.emplace_back(term.first, term.second * tAct);
    }

    return totalConfig;
}

//Decodes R bases into the a fermionic string
std::vector<CircuitConfig> baseCircuitCodeSeparate(const std::string& name, const Vector& ts,
                                                   const std::map<std::string, std::vector<Term>>& basisDict, bool left) {
    //Returns Fermionic Strings

    // Check for composites by counting the occurrences of 'R'
    size_t countR = std::count(name.begin(), name.end(), 'R');
    std::vector<std::string> nameList;
    std::vector<CircuitConfig> configs;

    if (countR > 1) { // Composite case
        size_t midPos = name.rfind('R') - 1;
        std::string name0 = name.substr(0, midPos);
        std::string name1 = name.substr(midPos);

        if (left) {
            nameList = {name1, name0};
        } else {
            nameList = {name0, name1};
        }

        // Process each component of the composite
        for (size_t i = 0; i < nameList.size(); ++i) {
            configs.push_back(baseCircuitCodeSingle(nameList[i], ts[i], basisDict, left));
        }
    } else { // Single operator case
        nameList = {name};
        configs.push_back(baseCircuitCodeSingle(name, ts[0], basisDict, left));
    }

    return configs;
}

//Function to compute the right R
std::vector<CircuitConfig> computeRightR(const std::string& nameRight, const Vector& tsRight,
                                         const std::map<std::string, std::vector<Term>>& basisDict, bool isLeft = false, bool debug = false) {
    if (nameRight.find("R0") != std::string::npos && nameRight.length() > 2) {
        throw std::invalid_argument("Do not include R0 in composite or with a sign. Not implemented.");
    }

    size_t countR = std::count(nameRight.begin(), nameRight.end(), 'R'); // Count R's to determine if composite
    std::vector<CircuitConfig> rightSubs;

    if (countR > 1) { // Composite case
        rightSubs = baseCircuitCodeSeparate(nameRight, tsRight, basisDict, isLeft);
        if (debug) {
            std::cout << "Composite case processed." << std::endl;
        }
    } else { // Single operator case
        rightSubs.push_back(baseCircuitCodeSeparate(nameRight, tsRight, basisDict, isLeft)[0]);
        if (debug) {
            std::cout << "Single case processed." << std::endl;
        }
    }

    return rightSubs;
}

//Function to compute the left R
std::vector<CircuitConfig> computeLeftR(const std::string& nameLeft, const Vector tsLeft,
                                        const std::map<std::string, std::vector<Term>>& basisDict, bool debug = false) {
    //Check if it's composite
    size_t countR = std::count(nameLeft.begin(), nameLeft.end(), 'R');
    std::vector<CircuitConfig> leftSubs;

    if (countR > 1) { //Composite case
        leftSubs = computeRightR(nameLeft, tsLeft, basisDict, true, debug);
        if (debug) {
            std::cout << "Composite case processed in computeLeftR." << std::endl;
        }
    } else { //Single operator case
        leftSubs.push_back(computeRightR(nameLeft, tsLeft, basisDict, true, debug)[0]);
        if (debug) {
            std::cout << "Single case processed in computeLeftR." << std::endl;
        }
    }
    return leftSubs;
}

//Generates the t values in the correct order depending if its the left or right operator
void Gen_tLR(Vector& tLR, const std::string& name, const Vector& ts, bool left) { 
    size_t count_R = std::count(name.begin(), name.end(), 'R');
    //Single operator case
    if (count_R == 1) {
        tLR[0] = ts[std::stoi(name.substr(name.find('R') + 1))];
    //Composite case
    } else if (count_R == 2) {
        std::string filtered = std::regex_replace(name, std::regex("[+-]"), "");
        std::vector<int> indices;
        size_t pos = 0;
        while ((pos = filtered.find('R', pos)) != std::string::npos) {
            indices.push_back(ts[std::stoi(filtered.substr(pos + 1, filtered.find_first_not_of("0123456789", pos + 1) - pos - 1))]);
            pos += 2; // Move past 'R' and the digit
        }
        if (indices.size() == 2) {
            //If its the left bases
            if(left)
                tLR = {indices[1], indices[0]};
            //If its the right bases
            else
                tLR = {indices[0], indices[1]};
        } else {
            throw std::runtime_error("Invalid basis name: " + name);
        }
    } else {
        throw std::runtime_error("Invalid basis name: " + name);
    }
}

class PauliOp {
public:
    std::map<std::string, double> terms;

    void addTerm(const std::string& term, double coefficient) {
        if (terms.find(term) == terms.end()) {
            terms[term] = coefficient;
        } else {
            terms[term] += coefficient;
        }
    }
    void display() const {
        for (const auto& term : terms) {
            std::cout << term.first << " * " << term.second << "\n";
        }
    }
};

//Functions to define individual Pauli matrices
ComplexMatrix2D pauliI() {
    ComplexMatrix2D I;
    I << 1, 0,
         0, 1;
    return I;
}

ComplexMatrix2D pauliX() {
    ComplexMatrix2D X;
    X << 0, 1,
         1, 0;
    return X;
}

ComplexMatrix2D pauliY() {
    ComplexMatrix2D Y;
    Y << 0, -std::complex<double>(0, 1),
         std::complex<double>(0, 1), 0;
    return Y;
}

ComplexMatrix2D pauliZ() {
    ComplexMatrix2D Z;
    Z << 1, 0,
         0, -1;
    return Z;
}

// Create a map to access the Pauli matrices
std::map<char, ComplexMatrix2D> createPauliMap() {
    std::map<char, ComplexMatrix2D> pauliMap;
    pauliMap['I'] = pauliI();
    pauliMap['X'] = pauliX();
    pauliMap['Y'] = pauliY();
    pauliMap['Z'] = pauliZ();
    return pauliMap;
}

// Function to convert a single Pauli string to a matrix
ComplexMatrix pauliStringToMatrix(const std::string& pauliString) {
    //Create a map for gates to their respective matrices
    auto pauliMap = createPauliMap();
    size_t numQubits = pauliString.length();
    size_t matrixDim = 1 << numQubits; //2^numQubits
    //Fill with identity to start
    ComplexMatrix result = ComplexMatrix::Identity(matrixDim, matrixDim);

    for (size_t i = 0; i < numQubits; ++i) {
        ComplexMatrix2D currentMatrix = pauliMap[pauliString[i]];
        ComplexMatrix fullMatrix = ComplexMatrix::Identity(1 << i, 1 << i).kron(currentMatrix);
        result = result.kron(fullMatrix);
    }
    return result;
}

//Converts a PauliOp to a full matrix
ComplexMatrix pauliOpToMatrix(const PauliOp& pauliOp, int numQubits) {
    size_t matrixDim = 1 << numQubits; //2^numQubits
    ComplexMatrix  result = ComplexMatrix ::Zero(matrixDim, matrixDim);

    for (const auto& term : pauliOp.terms) {
        const std::string& pauliString = term.first;
        double coefficient = term.second;
        result += coefficient * pauliStringToMatrix(pauliString);
    }
    return result;
}

//Maps fermionic creation and annihilation operators to Pauli strings
std::pair<std::string, std::string> fermionicToPauli(int index, int numQubits) {
    std::string creation(numQubits, 'I');
    std::string annihilation(numQubits, 'I');

    // Add the Z strings up to index
    for (int i = 0; i < index; ++i) {
        creation[i] = 'Z';
        annihilation[i] = 'Z';
    }

    // Add X and Y Pauli operators at the index position
    creation[index] = 'X';
    annihilation[index] = 'Y';

    return std::make_pair(creation, annihilation);
}

//Jordan-Wigner Transformation Function
PauliOp JWTransform(const CircuitConfig& circuitConfig, int numQubits) {
    PauliOp pauliOp;

    for (const auto& fop : circuitConfig) {
        std::string fermionStr = fop.first;
        double coefficient = fop.second;

        for (size_t i = 0; i < fermionStr.length(); ++i) {
            if (isdigit(fermionStr[i])) {
                int index = fermionStr[i] - '0';
                if (index < 0 || index >= numQubits) {
                    throw std::out_of_range("Index out of range for the number of orbitals");
                }

                // Determine if it's a creation or annihilation operator
                bool isCreation = fermionStr[i - 1] == '+';
                std::pair<std::string, std::string> pauliStrings = fermionicToPauli(index, numQubits);

                // Add the appropriate Pauli string to the operator
                pauliOp.addTerm(isCreation ? pauliStrings.first : pauliStrings.second, coefficient);
            }
        }
    }
    return pauliOp;
}
//Trotterization function
PauliOp trotterize(const PauliOp& op, int numSlices) {
    PauliOp trotterized;
    for (const auto& term : op.terms) {
        trotterized.addTerm(term.first, term.second / numSlices);
    }
    return trotterized;
}

//Placeholder function in case of a need for an evaluation step for hamiltonian matrices
void Evaluate(){

};
//Function to compute matrix elements classically
void ComputeSH(ComplexMatrix HFState, ComplexMatrix HamMat, Matrix& S, Matrix& H, const std::vector<std::string>& bases, const Vector& ts, 
               const std::map<std::string, CircuitConfig>& basisDict, int numOrbitals, int matLen, bool debug = false) {
    //Start with the left side of <R| |R> as the outer loop to fill in the matrices
    //Initialize some strings for the name of the right and left bases
    std::string nar, nal;
    //Initialize some intermediate matrices
    ComplexMatrix HMatMeasurement, HMatMeasurable, HMatExpectation, SMatMeasurement, SMatMeasurable, SMatExpectation;
    ComplexMatrix SObs, HObs;
    //Initialize the final complex matrices, the real part wil be place in S and H
    ComplexMatrix SComp, HComp;
    //Initialize circuitconfig objects for the left and right bases
    std::vector<CircuitConfig> leftR, rightR;
    //The t's that will correspond with the right and left bases respectively
    Vector tr, tl;
    //Left and right pauli operators and matrix representations of the Pauli operators
    PauliOp leftPauliOp, rightPauliOp;
    ComplexMatrix leftMatrix, rightMatrix;
    //Outer loop & left basis(bases) calculation
    for (int i = 0; i < matLen; i++) {
        nal = bases[i]; //Name left
        tl(-1, -1); //Default to -1 no value
        Gen_tLR(tl, nal, ts, true); //Get the pair of t's in order for the left side
        //Output the tuple
        //std::cout << "Tuple: (" << tl.first << ", " << (tl.second == -1 ? "" : std::to_string(tl.second)) << ")\n";
        //Generate fermionic strings for left side
        leftR = computeLeftR(nal, tl, basisDict, debug);
        leftPauliOp = JWTransform(leftR, numOrbitals);
        leftMatrix = pauliOpToMatrix(leftPauliOp);
        
        //Take the right side of <R| |R> as the inner loop to fill in the matrices
        //Inner loop & right basis(bases) calculation
        for (int j = 0; j < matLen; j++) {
            //Name of right basis(bases)
            nar = basis_set[j];
            //t values that we are interested in for this(these) particular basis(bases)
            tr(-1, -1);
            //Get the t or pair of t's in order for the right side
            Gen_tLR(tr, nar, ts, false);

            HMatMeasurement = HamMat.adjoint(); 
            //Diagonal entries in S and H matrices
            if(nal == nar){ 
                SComp[i,j] = 1;
                //Simplify a bit in R0 basis case
                if(nal == 'R0'){
                    HMatMeasurable = HMatMeasurement  * HFState;
                }
                else{
                    rightR = computeRightR(nar, tr, basisDict, debug);
                    rightPauliOp = JWTransform(rightR, numOrbitals);
                    rightMatrix = pauliOpToMatrix(rightPauliOp);
                    HMatMeasurable = HMatMeasurement * rightMatrix * HFState;
                }
                HMatExpectation = MatrixExpectation(HMatMeasurable);
                HComp[i,j] = Evaluate(HMatExpectation);
            }
            //Off diagonal
            else{ 
                //Simplify a bit in R0 basis case
                if(nar == 'R0'){
                    SObs = leftMatrix;
                    HObs = leftMatrix * HamMat;
                }
                //Simplify a bit in R0 basis case
                else if (nal == 'R0'){
                    rightR = computeRightR(nar, tr, basisDict, debug);
                    rightPauliOp = JWTransform(rightR, numOrbitals);
                    rightMatrix = pauliOpToMatrix(rightPauliOp);
                    SObs = rightMatrix;
                    HObs = HamMat * rightMatrix;
                }
                else{
                    rightR = computeRightR(nar, tl, basisDict, debug);
                    rightPauliOp = JWTransform(rightR, numOrbitals);
                    rightMatrix = pauliOpToMatrix(rightPauliOp);
                    SObs = leftMatrix * rightMatrix;
                    HObs = leftMatrix * HamMat * rightMatrix;
                }
                //S matrix
                SMatMeasurement = SObs.adjoint();
                SMatMeasurable = SMatMeasurement * HFState;
                SMatExpectation = MatrixExpectation(SMatMeasurable);
                SComp[i,j] = Evaluate(SMatExpectation);
                //H matrix
                HMatMeasurement = Hobs.adjoint();
                HMatMeasurable = HMatMeasurement * HFState;
                HMatExpectation = MatrixExpectation(HMatMeasurement);
                HComp[i,j] = Evaluate(HMatExpectation);
            }
        }
    }
    S = SComp.real();
    H = HComp.real();
}

//for each Vp;
//Step 7 & 8
//Step 9 & 10: Evaluate matrix elements on a quantum device
//Matrix evaluatedMatrixElements = EvaluateMatrixElements(matrixElements);
//Function to evaluate matrix elements on a quantum device
void QuantumSH(ComplexMatrix HFState, ComplexMatrix HamMat, Matrix& S, Matrix& H, const std::vector<std::string>& bases, const Vector& ts, 
               const std::map<std::string, CircuitConfig>& basisDict, int numOrbitals, int matLen, bool debug = false) {
    //Start with the left side of <R| |R> as the outer loop to fill in the matrices
    //Initialize some strings for the name of the right and left bases
    std::string nar, nal;
    //Initialize some intermediate matrices
    ComplexMatrix HMatMeasurement, HMatMeasurable, HMatExpectation, SMatMeasurement, SMatMeasurable, SMatExpectation;
    ComplexMatrix SObs, HObs;
    //Initialize the final complex matrices, the real part wil be place in S and H
    ComplexMatrix SComp, HComp;
    //Initialize circuitconfig objects for the left and right bases
    std::vector<CircuitConfig> leftR, rightR;
    //The t's that will correspond with the right and left bases respectively
    Vector tr, tl;
    //Left and right pauli operators and matrix representations of the Pauli operators
    PauliOp leftPauliOp, rightPauliOp;
    ComplexMatrix leftMatrix, rightMatrix;
    //Outer loop & left basis(bases) calculation
    for (int i = 0; i < matLen; i++) {
        nal = bases[i]; //Name left
        tl(-1, -1); //Default to -1 no value
        Gen_tLR(tl, nal, ts, true); //Get the pair of t's in order for the left side
        //Output the tuple
        //std::cout << "Tuple: (" << tl.first << ", " << (tl.second == -1 ? "" : std::to_string(tl.second)) << ")\n";
        //Generate fermionic strings for left side
        leftR = computeLeftR(nal, tl, basisDict, debug);
        leftPauliOp = JWTransform(leftR, numOrbitals);
        leftMatrix = pauliOpToMatrix(leftPauliOp);
        
        //Take the right side of <R| |R> as the inner loop to fill in the matrices
        //Inner loop & right basis(bases) calculation
        for (int j = 0; j < matLen; j++) {
        }
    }
}


