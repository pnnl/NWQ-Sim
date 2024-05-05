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
#include "Eigen/Dense"


//Using dynamic Eigen MatrixXd to streamline calculations
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using ComplexMatrix = Eigen::MatrixXcd;
using ComplexVector = Eigen::VectorXcd;

typedef std::pair<std::string, double> Term; //Terms in a fermionic string
typedef std::vector<Term> CircuitConfig;

int generateRandomNumber(uint64_t seed) {
    // Create a PCG random number generator initialized with the given seed
    pcg64 rng(seed);

    std::uniform_int_distribution<int> dist(0, 100); // Distribution that maps to the range [0,100]
    return dist(rng); // Generate and return the random number
}

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
    std::string baseStr = singleName.substr(1); // 'R<num>', e.g., 'R1'
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

// Function to compute the left R
std::vector<CircuitConfig> computeLeftR(const std::string& nameLeft, const Vector tsLeft,
                                        const std::map<std::string, std::vector<Term>>& basisDict, bool debug = false) {
    size_t countR = std::count(nameLeft.begin(), nameLeft.end(), 'R'); // Check if it's composite
    std::vector<CircuitConfig> leftSubs;

    if (countR > 1) { // Composite case
        leftSubs = computeRightR(nameLeft, tsLeft, basisDict, true, debug);
        if (debug) {
            std::cout << "Composite case processed in computeLeftR." << std::endl;
        }
    } else { // Single operator case
        leftSubs.push_back(computeRightR(nameLeft, tsLeft, basisDict, true, debug)[0]);
        if (debug) {
            std::cout << "Single case processed in computeLeftR." << std::endl;
        }
    }

    return leftSubs;
}

void Gen_tLeft(Vector& tl, const std::string& nal, const Vector& ts) {
    size_t count_R = std::count(nal.begin(), nal.end(), 'R');
    if (count_R == 1) {
        tl[0] = ts[std::stoi(nal.substr(nal.find('R') + 1))];
    } else if (count_R == 2) {
        std::string filtered = std::regex_replace(nal, std::regex("[+-]"), "");
        std::vector<int> indices;
        size_t pos = 0;
        while ((pos = filtered.find('R', pos)) != std::string::npos) {
            indices.push_back(ts[std::stoi(filtered.substr(pos + 1, filtered.find_first_not_of("0123456789", pos + 1) - pos - 1))]);
            pos += 2; // Move past 'R' and the digit
        }
        if (indices.size() == 2) {
            tl = {indices[1], indices[0]};
        } else {
            throw std::runtime_error("Invalid basis name: " + nal);
        }
    } else {
        throw std::runtime_error("Invalid basis name: " + nal);
    }
}

class PauliOperator {
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
            std::cout << "exp(" << term.second << " * " << term.first << ")\n";
        }
    }
};

//Jordan-Wigner transformation
PauliOperator jordanWignerTransform(const CircuitConfig& circuitConfig, int numOrbitals) {
    PauliOperator pauliOperator;
    for (const auto& fOp : circuitConfig) {
        std::string pauliString = "";
        int lastIndex = -1;

        for (size_t i = 0; i < fOp.first.length(); ++i) {
            if (isdigit(fOp.first[i])) {
                int index = fOp.first[i] - '0';
                for (int j = lastIndex + 1; j < index; ++j) {
                    pauliString += "I";
                }
                pauliString += (fOp.first[i-1] == '+' ? "X" : "Y");
                for (int k = index + 1; k < numOrbitals; ++k) {
                    pauliString += "Z";
                }
                lastIndex = index;
            }
        }

        pauliOperator.addTerm(pauliString, fop.second);
    }
    return pauliOperator;
}

// Trotterization function remains the same
PauliOperator trotterize(const PauliOperator& op, int numSlices) {
    PauliOperator trotterized;
    for (const auto& term : op.terms) {
        trotterized.addTerm(term.first, term.second / numSlices);
    }
    return trotterized;
}

//Function to compute matrix elements classically
void ComputeSH(ComplexMatrix& S, ComplexMatrix& H, const std::vector<std::string>& bases, const Vector& ts, 
               const std::map<std::string, std::vector<Term>>& basisDict, int numOrbitals, int matLen, bool debug = false) {
    for (int i = 0; i < matLen; i++) {
        std::string nal = bases[i];
        Vector tl(-1, -1); //Default to -1 no value
        Gen_tLeft(tl, nal, ts);
        //Output the tuple
        //std::cout << "Tuple: (" << tl.first << ", " << (tl.second == -1 ? "" : std::to_string(tl.second)) << ")\n";
        //Generate fermionic strings for LHS
        std::vector<CircuitConfig> leftR = computeLeftR(nal, tl, basisDict, debug);
        int trotterSlices = 2;  // Define number of Trotter slices
        for(int j = 0; j < (int)leftR.size(); j++){
            PauliOperator jwOp = jordanWignerTransform(leftR[j], numOrbitals);
            PauliOperator trotterOp = trotterize(jwOp, trotterSlices);
        }
    
        //for each Vp;
            // Step 7 & 8
            // Step 9 & 10: Evaluate matrix elements on a quantum device
            //Matrix evaluatedMatrixElements = EvaluateMatrixElements(matrixElements);
    }
}

// Function to evaluate matrix elements on a quantum device
void EvaluateSH(ComplexMatrix& S, ComplexMatrix& H, Vector& bases,const Vector& ts) {
    // placeholder for quantum computation
}


