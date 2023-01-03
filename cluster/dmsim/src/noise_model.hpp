// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Muqing Zheng and Ang Li
// Pacific Northwest National Laboratory(PNNL), U.S.
// GitHub repo: http://www.github.com/pnnl/SV-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
#ifndef NOISE_MODEL_HPP_
#define NOISE_MODEL_HPP_

#include <iostream> // std::cout, std::scientific, for prints only
#include <limits> // std::numeric_limits, for precision in prints only
#include <complex> // std::complex, std::conj, std::real, std::imag
#include <stdexcept> // std::invalid_argument
#include <cmath> // std::sqrt, std::exp, std:isinf, INFINITY, std::log2
#include <bitset> //  std::bitset, used for the tensor product of superoperators

namespace NWQSim {

// --------------------------------- References ---------------------------------
// 1. http://arxiv.org/abs/1111.6950 (which is also one of the main references used in Qiskit Aer)
// 2. Quantum Computation and Quantum Information (10th Anniversary Edition) by Michael A. Nielsen & Isaac L. Chuang
// 3. Qiskit Aer repository  https://github.com/Qiskit/qiskit-aer

//--------------------------------- Constants --------------------------------------
// Pauli matrices I, X, Y, Z
std::complex<double> pauliMats[4][2][2] = {
    {{{1,0},{0, 0}}, {{0,0},{ 1,0}}}, //Pauli-I
    {{{0,0},{1, 0}}, {{1,0},{ 0,0}}}, //Pauli-X
    {{{0,0},{0,-1}}, {{0,1},{ 0,0}}}, //Pauli-Y
    {{{1,0},{0, 0}}, {{0,0},{-1,0}}}  //Pauli-Z
};

// Thermal relaxation error operators
std::complex<double> trMats[3][2][2] = {
    {{{1,0},{0,0}}, {{0,0},{ 1,0}}}, //Pauli-I
    {{{1,0},{0,0}}, {{0,0},{-1,0}}}, //Pauli-Z
    {{{1,0},{0,0}}, {{0,0},{ 0,0}}}, //0-State
};
// Also define Amplitude and Phase Damping Kraus operators inside the function

// ---------------------------------------- Helper Functions ---------------------------------------------
/**
 *  Printing out gate according to scientific format.
*/
void printGate(std::complex<double>* mat, int dim) {
    std::cout.precision(std::numeric_limits<double>::max_digits10 - 1);
    for(int x = 0; x < dim; x++)  {
        for(int y = 0; y < dim; y++) {
            std::cout <<  std::scientific
                << "("
                << std::real(mat[x*dim+y])
                << (std::imag(mat[x*dim+y]) >= 0.0 ? "+" : "")
                << std::imag(mat[x*dim+y])
                << "j)   ";
        }
    std::cout << std::endl;
    }
}

/**
 * Kronecker product with scalar multiplication (if necessary) between two square matrices, A \otimes B = res, where dim is the dimension of A and B
 */
void kronProd(std::complex<double>*  A, std::complex<double>*  B, double c, std::complex<double>*  res, int dim) {
    int bigdim = dim*dim; // dimsion of res
    for (int r = 0; r < dim; r++) {
        for (int s = 0; s < dim; s++) {
            for (int v = 0; v < dim; v++) {
                for (int w = 0; w < dim; w++) {
                    res[(dim*r+v)*bigdim + (dim*s+w)] = c * A[r*dim+s] * B[v*dim+w];
                }
            }
        }
    }
}

/**
 * Dot product between square matrices A.B = res 
 */
void dotProd(std::complex<double>*  A, std::complex<double>*  B, std::complex<double>*  res, int dim) {
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            res[i*dim+j] = 0;
            for (int k = 0; k < dim; k++) {
                res[i*dim+j] += A[i*dim+k] * B[k*dim+j];
            }
        }
    }
}

/**
 * Conjugate transpose A^\dagger = res
 */
void dagger(std::complex<double>*  A, std::complex<double>*  res, int dim) {
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            res[j*dim+i] = std::conj(A[i*dim+j]);
        }
    }
}

/**
 * Conjugate \bar{A} = res
 */
void conjbar(std::complex<double>*  A, std::complex<double>*  res, int dim) {
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            res[i*dim+j] = std::conj(A[i*dim+j]);
        }
    }
}

/**
 * @brief Transfer the 1 Qubit Kraus representation of a map to a Liouville-Superoperator. 
 *        Ref: (4.5) in http://arxiv.org/abs/1111.6950
 * 
 * @param kr_chnnl_ops 3D array that stores 1-Qubit Kraus operators, std::complex<double> kr_chnnl_ops[n_ops][2][2]
 * @param sp_ops  Where to store the Liouville-Superoperator, std::complex<double> sp_ops[4][4]
 * @param n_ops number of Kraus operators
 */
void kr2sp1Q(std::complex<double> kr_chnnl_ops[][2][2], std::complex<double> sp_ops[4][4], int n_ops) {
    const int mat_dim = 2;
    const int big_dim = mat_dim*mat_dim;
    std::complex<double> kbar[mat_dim][mat_dim] = {}; // to save conjugate of K_\alpha
    std::complex<double> KK_kp[big_dim][big_dim] = {}; // to save \bar{K_\alpha} \otimes K_\alpha
    for (int a = 0; a < n_ops; a++) {// iterate over all Kraus operators
        // Compute conjugate of Kraus operator
        conjbar(kr_chnnl_ops[a][0], kbar[0], mat_dim);
        // \bar{K_\alpha} \otimes K_\alpha
        kronProd(kbar[0], kr_chnnl_ops[a][0], 1.0, KK_kp[0], mat_dim);
        // Added the tensor product to 
        for (int i = 0; i < big_dim; i++) {
            for (int j = 0; j < big_dim; j++) {
                if (a == 0) {
                    sp_ops[i][j] = KK_kp[i][j];
                } else {
                    sp_ops[i][j] += KK_kp[i][j];
                }
            }
        }
    }
}

// 2-Qubit version
// The only difference is different mat_dim
void kr2sp2Q(std::complex<double> kr_chnnl_ops[][4][4], std::complex<double> sp_ops[16][16], int n_ops) {
    const int mat_dim = 4;
    const int big_dim = mat_dim*mat_dim;
    std::complex<double> kbar[mat_dim][mat_dim] = {}; // to save conjugate of K_\alpha
    std::complex<double> KK_kp[big_dim][big_dim] = {}; // to save \bar{K_\alpha} \otimes K_\alpha
    for (int a = 0; a < n_ops; a++) {// iterate over all Kraus operators
        // Compute conjugate of Kraus operator
        conjbar(kr_chnnl_ops[a][0], kbar[0], mat_dim);
        // \bar{K_\alpha} \otimes K_\alpha
        kronProd(kbar[0], kr_chnnl_ops[a][0], 1.0, KK_kp[0], mat_dim);
        // Added the tensor product to 
        for (int i = 0; i < big_dim; i++) {
            for (int j = 0; j < big_dim; j++) {
                if (a == 0) {
                    sp_ops[i][j] = KK_kp[i][j];
                } else {
                    sp_ops[i][j] += KK_kp[i][j];
                }
            }
        }
    }
}

//--------------------------- For tensor product of Superoperators --------------------------------------
const int int_bits_size = 4;
/**
 * @brief Swith two bits in a bitset from two given position
 * 
 * @param origBitset a length-(int_bits_size) bitset
 * @param pos1 position
 * @param pos2 position
 * @return std::bitset<int_bits_size>, the bitset with swithed positions
 */
std::bitset<int_bits_size> swapBits(std::bitset<int_bits_size> origBitset, int pos1, int pos2) {
    std::bitset<int_bits_size> res (origBitset);
    res.set(pos1, origBitset[pos2]);
    res.set(pos2, origBitset[pos1]);
    return res;
}

/**
 * @brief Generate a 2-qubit superoperator from two 1-qubit superoperators. 
 *        A combination of Kronecker product and reshuffle. The general reshuffle logic is that,
 *        for example, Use bits to represent the coordinates of entries in a 4-by-4 matrix C, where C = A \otimes B
 *        Then, we have formula C_{ij,kl} = A_{i,k} \times B_{j,l} where i,j,k,l \in \{0,1\}.
 *        So In terms of superoperators, for example, A = \bar{K}_1 \otimes K_1, B = \bar{K}_2 \otimes K_2, are two superoperators applied on two different qubits,
 *        then the superoperator represents the combination of both operators are C = \bar{(K_1 \\otimes K_2)} \otimes (K_1 \\otimes K_2)
 *        while A \otimes B = (\bar{K}_1 \otimes K_1) \otimes (\bar{K}_2 \otimes K_2). 
 *        So if still using bits to represent coordinates of entries, shifting bits at position 1 and 2 (rightmost bit is position 0) of coordinates of all entries of A \otimes B will give C
 * 
 * @param A 4-by-4 array, std::complex<double>, an 1-qubit superoperator
 * @param B 4-by-4 array, std::complex<double>, an 1-qubit superoperator
 * @param res 16-by-16 array, std::complex<double>, an 2-qubit superoperator as the combined operator from A and B
 */
void SuOpskronProd1q(std::complex<double> A[4][4], std::complex<double> B[4][4], std::complex<double> res[16][16]) {
    // Compute regular tensor product first
    std::complex<double> dum_res[16][16] = {};
    kronProd(A[0], B[0], 1.0, dum_res[0], 4);
    // printGate(dum_res[0], 16);
    // Reshuffle
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            // transfer coordnates to bits
            std::bitset<int_bits_size> bit_i (i);
            std::bitset<int_bits_size> bit_j (j);
            // Swap pos 1 and pos 2 (starting from rightmost)
            std::bitset<int_bits_size> swapped_bit_i (swapBits(bit_i, 1, 2));
            std::bitset<int_bits_size> swapped_bit_j (swapBits(bit_j, 1, 2));
            // Convert to int
            int new_i = static_cast<int>(swapped_bit_i.to_ulong());
            int new_j = static_cast<int>(swapped_bit_j.to_ulong());
            // Swap the values
            res[new_i][new_j] = dum_res[i][j];
        }
    }
}

//------------------------------------ Depolarizing Errors ------------------------------------------------------
/**
 * @brief Apply depolarizing channel to a given idea unitary gate under rate lambda
 *        D_\lambda(\rho) = (1-\lambda)\rho + \lambda/d I
 *        where d is the dimension of the state space (2^num_of_qubits).
 *        Note that lambda = 4/3 p where p is the probability to have one of the X,Y,Z error
 * 
 *        So the ideal gate operation U \rho U^\dagger becomes
 *         $D(U \rho U^\dagger) = (1-0.75\lambda)U \rho U^\dagger + (\lambda/4)\sum_{E \in {X,Y,Z}}(EU)\rho(EU)^\dagger$
 * 
 * @param lam,  0<= lam <= 1/(d^2-1), where d is  the dimension of matrix (i.e., 2^n)
 * @param ideal_gate 2-by-2 matrix (1-qubit gate)
 * @param sp  where superoperator is saved, either the noise superoperator itself, or combining with the ideal gate
 * @param add_to_gate if true, then sp will be saved by a superoperator for the noisy gate; if false, then only noise superoperator will be saved to sp.
 */
void addDepErr1Q(double lam, 
                 std::complex<double> ideal_gate[2][2], 
                 std::complex<double> sp[4][4],
                 bool add_to_gate = true) {

    if ( lam < 0.0 ) { 
        throw std::invalid_argument( "Lambda should be non-negative." );
    }     

    // Compute super operator of ideal gate
    std::complex<double> ideal_gate_conj[2][2] = {}; // conjugate of ideal gate
    conjbar(ideal_gate[0], ideal_gate_conj[0], 2);
    std::complex<double> ideal_gate_sp[4][4] = {}; // super operator of ideal gate
    kronProd(ideal_gate_conj[0], ideal_gate[0], 1.0, ideal_gate_sp[0], 2);

    if (lam == 0.0 && add_to_gate) {
        std::copy(&ideal_gate_sp[0][0], &ideal_gate_sp[0][0]+16,&sp[0][0]); // could use more efficient method
        return;
    } else if (lam == 0.0) { // identity
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i == j) {
                    sp[i][j] = 1.0;
                } else {
                    sp[i][j] = 0.0;
                }
            }
        }
        return;
    }


    // Compute noise parameters
    double idea_coef = std::sqrt(1 - 0.75 * lam); // 0.75 is by absorbing the case when error operator is I
    double err_coef = 0.5 * std::sqrt(lam);
    // Computer Kraus operators
    std::complex<double> chnnl_ops[4][2][2] = {};
    std::complex<double> noise_sp[4][4] = {};
    //  $D(U \rho U^\dagger) = (1-0.75\lambda)U \rho U^\dagger + \sum_{E \in {X,Y,Z}}(\sqrt{\lambda/4}EU)\rho(\sqrt{\lambda/4}EU)^\dagger$
    // So 4 operators: sqrt(1-0.75lam)I, and  \sqrt(lam/4)EU for E \in {X,Y,Z}
    // sqrt(1-0.75lam)I
    chnnl_ops[0][0][0] = idea_coef;
    chnnl_ops[0][1][1] = idea_coef;
    // \sqrt(lam/4)X
    chnnl_ops[1][0][1] = err_coef;
    chnnl_ops[1][1][0] = err_coef;
    // \sqrt(lam/4)Y
    std::complex<double> sqno(0,1); // sqrt of -1
    chnnl_ops[2][0][1] = -1.0 * sqno * err_coef;
    chnnl_ops[2][1][0] = sqno * err_coef;
    // \sqrt(lam/4)Z
    chnnl_ops[3][0][0] = err_coef;
    chnnl_ops[3][1][1] = -1.0 * err_coef;

    // Transfer to a Superoperator
    kr2sp1Q(chnnl_ops, noise_sp, 4);
    // compute final superoperator
    if (add_to_gate) {
        dotProd(noise_sp[0], ideal_gate_sp[0], sp[0], 4); // dot product the noise superoperator and ideal-gate superoperator
    } else { // if just need superoperator of noise
        std::copy(&noise_sp[0][0], &noise_sp[0][0]+16, &sp[0][0]);
    }
}

/**
 * @brief Apply 2-qubit depolarizing channel to a given idea unitary gate under rate lambda
 *        D_\lambda(\rho) = (1-\lambda)\rho + \lambda/d I
 *        where d is the dimension of the state space (2^num_of_qubits).
 *        TODO: verify the relation between p and lambda
 *        TODO: Did not multiply the ideal gate, fix it! (need to write a dot product)
 * @param lam,  0<= lam <= 1/(d^2-1)
 * @param ideal_gate 4-by-4 matrix (2-qubit gate)
 * @param sp  where superoperator is saved, either the noise superoperator itself, or combining with the ideal gate
 * @param add_to_gate if true, then sp will be saved by a superoperator for the noisy gate; if false, then only noise superoperator will be saved to sp.
 */
void addDepErr2Q(double lam, 
                 std::complex<double> ideal_gate[4][4], 
                 std::complex<double> sp[16][16],
                 bool add_to_gate = true) {
    if ( lam < 0.0 ) { 
        throw std::invalid_argument( "Lambda should be non-negative." );
    }              
    // Compute super operator of ideal gate
    std::complex<double> ideal_gate_conj[4][4] = {}; // conjugate of ideal gate
    conjbar(ideal_gate[0], ideal_gate_conj[0], 4); 
    std::complex<double> ideal_gate_sp[16][16] = {}; // super operator of ideal gate
    kronProd(ideal_gate_conj[0], ideal_gate[0], 1.0, ideal_gate_sp[0], 4);

    if (lam == 0.0 && add_to_gate) {
        std::copy(&ideal_gate_sp[0][0], &ideal_gate_sp[0][0]+256,&sp[0][0]); // could use more efficient method
        return;
    } else if (lam == 0.0) { // identity
        for (int i = 0; i < 16; i++) {
            for (int j = 0; j < 16; j++) {
                if (i == j) {
                    sp[i][j] = 1.0;
                } else {
                    sp[i][j] = 0.0;
                }
            }
        }
        return;
    }

    // Compute noise parameters
    double idea_coef = std::sqrt(1 - 0.9375*lam);
    double err_coef = 0.25 * std::sqrt(lam);
    std::complex<double> chnnl_ops[16][4][4] = {};
    std::complex<double> noise_sp[16][16] = {};
    // D(U \rho U^\dagger) = (1-\lambda)U \rho U^\dagger + \sum_{E \in {I,X,Y,Z}^2}(\sqrt{p/16}EU)\rho(\sqrt{p/16}EU)^\dagger
    // So 16 operators: sqrt(1-0.9375\lambda)I, and  \sqrt(lam/16)EU for E \in {I,X,Y,Z}^2 and E \neq II
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i == 0 && j == 0) { // identity
                kronProd(pauliMats[i][0], pauliMats[j][0], idea_coef, chnnl_ops[i*4+j][0], 2);
            } else {
                kronProd(pauliMats[i][0], pauliMats[j][0], err_coef, chnnl_ops[i*4+j][0], 2);
            }
        }
    }
    // Transfer to a Superoperator
    kr2sp2Q(chnnl_ops, noise_sp, 16);
    // compute final superoperator
    if (add_to_gate) {
        dotProd(noise_sp[0], ideal_gate_sp[0], sp[0], 16); // dot product the noise superoperator and ideal-gate superoperator
    } else { // if just need superoperator of noise
        std::copy(&noise_sp[0][0], &noise_sp[0][0]+256, &sp[0][0]);
    }
}

//------------------------------------- Thermal Relaxation Errors -----------------------------------------------------------------
/**
 * @brief Provide Thermal Relaxation (TR) channel in Kraus representation.
 *        Two cases: T1 <= T2 <= 2T1, use a Choi matrix, then transferred into Kaus
 *                   T2 < T1, then Kraus operators are I, Z, and reset (i.e., |0><0|)
 *                   T2 > 2T1, then set T2 = 2T1 like in Qiskit
 *                   TODO: Support non-zero excited_state_population if necessary (i.e., when temperature is not absolute zero)
 * 
 * @param gate_time time duration of operating a gate
 * @param T1 
 * @param T2 
 * @param ideal_gate 2-by-2 matrix (1-qubit gate)
 * @param sp  where superoperator is saved, either the noise superoperator itself, or combining with the ideal gate
 * @param add_to_gate if true, then sp will be saved by a superoperator for the noisy gate; if false, then only noise superoperator will be saved to sp.
 */
void addTRErr1Q(double gate_time, double T1, double T2,
                 std::complex<double> ideal_gate[2][2], 
                 std::complex<double> sp[4][4],
                 bool add_to_gate = true) {

    if (T1 <= 0.0 || T2 <= 0.0 || gate_time < 0.0 ) { 
        throw std::invalid_argument( "T1 or T2 should be positive. Gate time should be non-negative." );
    }

    // Compute super operator of ideal gate
    std::complex<double> ideal_gate_conj[2][2] = {}; // conjugate of ideal gate
    std::complex<double> ideal_gate_sp[4][4] = {}; // super operator of ideal gate
    if (add_to_gate) {
        conjbar(ideal_gate[0], ideal_gate_conj[0], 2);
        kronProd(ideal_gate_conj[0], ideal_gate[0], 1.0, ideal_gate_sp[0], 2);
    }
    
    // If gate time is 0, then there is no noise
    if (gate_time == 0.0 && add_to_gate) {
        std::copy(&ideal_gate_sp[0][0], &ideal_gate_sp[0][0]+16,&sp[0][0]); // could use more efficient method
        return;
    } else if (gate_time == 0.0) { // identity
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i == j) {
                    sp[i][j] = 1.0;
                } else {
                    sp[i][j] = 0.0;
                }
            }
        }
        return;
    }
    // Check effectiveness of T1 and T2
    double eff_T2 = 0.0; // effective T1
    if (T2 > 2*T1) {
        eff_T2 = 2*T1;
    } else {
        eff_T2 = T2;
    }
    // Compute parameters
    double p_rate1;
    double p_rate2;
    double p_reset;
    if (std::isinf(T1)) {
        p_rate1 = 1;
        p_reset = 0;
    } else {
        p_rate1 = std::exp(-1 * gate_time / T1);
        p_reset = 1 - p_rate1;
    }
    if (std::isinf(eff_T2)) {
        p_rate2 = 1;
    } else {
        p_rate2 = std::exp(-1 * gate_time / eff_T2);  
    }
    double sq_preset = std::sqrt(p_reset);
    if (T1 <= eff_T2) {     // T1 <= T2
        // // Construct Choi matrix
        // std::complex<double> tr_choi[4][4] = {
        //     {{1,0}, {0,0}, {0,0}, {p_rate2, 0}},
        //     {{0,0}, {0,0}, {0,0}, {0, 0}},
        //     {{0,0}, {0,0}, {p_reset,0}, {0, 0}},
        //     {{p_rate2,0}, {0,0}, {0,0}, {1-p_reset, 0}},
        // };
        // Make it superoperator by column reshuffle
        std::complex<double> noise_sp[4][4] = {
            {{1,0}, {0,0}, {0,0}, {p_reset,0}},
            {{0,0}, {p_rate2,0}, {0,0}, {0, 0}},
            {{0,0}, {0,0}, {p_rate2, 0}, {0, 0}},
            {{0,0}, {0,0}, {0,0}, {1-p_reset, 0}},
        };
        // compute final superoperator
        if (add_to_gate) {
            dotProd(noise_sp[0], ideal_gate_sp[0], sp[0], 4); // dot product the noise superoperator and ideal-gate superoperator
        } else { // if just need superoperator of noise
            std::copy(&noise_sp[0][0], &noise_sp[0][0]+16,&sp[0][0]);
        }
    } else { // T1 > T2
        std::complex<double> chnnl_ops[3][2][2] = {};
        double p_z = 0.5 * (1 - p_reset) * (1 - p_rate2/p_rate1);
        double p_iden = 1 - p_z - p_reset;
        double sq_pz = std::sqrt(p_z);
        double sq_piden = std::sqrt(p_iden);

        chnnl_ops[0][0][0] = sq_piden;
        chnnl_ops[0][1][1] = sq_piden;
        // \sqrt(p_z)Z
        chnnl_ops[1][0][0] = sq_pz;
        chnnl_ops[1][1][1] = -1.0 * sq_pz;
        // \sqrt(p_reset) |0><0|
        chnnl_ops[2][0][0] = sq_preset;

        // Transfer to a Superoperator of noise
        std::complex<double> noise_sp[4][4] = {}; // pure Superoperator without ideal gate
        kr2sp1Q(chnnl_ops, noise_sp, 3);
        noise_sp[0][3] = p_reset;
        // compute final superoperator
        if (add_to_gate) {
            dotProd(noise_sp[0], ideal_gate_sp[0], sp[0], 4); // dot product the noise superoperator and ideal-gate superoperator
        } else { // if just need superoperator of noise
            std::copy(&noise_sp[0][0], &noise_sp[0][0]+16, &sp[0][0]);
        }
    }
}

// The 2-qubit version of previous function.
// Essentially the Kronecker product of two 1-qubit thermal relaxation noise superoperator.
void addTRErr2Q(double gate_time_1, double T1_1, double T2_1,
                double gate_time_2, double T1_2, double T2_2,
                 std::complex<double> ideal_gate[4][4], 
                 std::complex<double> sp[16][16],
                 bool add_to_gate = true) {
    // Compute 1-qubit thermal relaxation noise superoperators
    std::complex<double> dummy_ideal_gate[2][2] = {{{1,0},{0,0}},{{0,0},{1,0}}}; // whatever that is should not matter, not used
    std::complex<double> gate1_sp[4][4] = {};
    std::complex<double> gate2_sp[4][4] = {};
    addTRErr1Q(gate_time_1, T1_1, T2_1, dummy_ideal_gate, gate1_sp, false);
    addTRErr1Q(gate_time_2, T1_2, T2_2, dummy_ideal_gate, gate2_sp, false);
    // Compute 2-qubit noise operator
    std::complex<double> noise_sp[16][16] = {};
    SuOpskronProd1q(gate2_sp, gate1_sp, noise_sp); // "Tensor Product (with reshuffle)" of two super operators
    if (!add_to_gate) { // if just need superoperator of noise
        std::copy(&noise_sp[0][0], &noise_sp[0][0]+256, &sp[0][0]);
        return;
    }
    // Compute super operator of ideal gate
    std::complex<double> ideal_gate_conj[4][4] = {}; // conjugate of ideal gate
    conjbar(ideal_gate[0], ideal_gate_conj[0], 4); 
    std::complex<double> ideal_gate_sp[16][16] = {}; // super operator of ideal gate
    kronProd(ideal_gate_conj[0], ideal_gate[0], 1.0, ideal_gate_sp[0], 4);
    // compute final superoperator
    dotProd(noise_sp[0], ideal_gate_sp[0], sp[0], 16); // dot product the noise superoperator and ideal-gate superoperator
}


//-------------------------------------------- Amplitude and Phase Damping Errors ----------------------------------------------------------
/**
 * @brief Produce amplitude and phase damping noise operator and apply it to a give 1-qubit gate.
 *        By setting either damping parameters to 0, it can give only amplitude or phase damping.
 *        References:
 *          Amplitude Damping: Nielsen and Chuang, Page 380, Equation (8.108)
 *          Phase Damping: Nielsen and Chuang, Page 384, Equation (8.127) and (8.128)
 * 
 *        For now, do NOT consider effect of dissipation to an environment at finite temperature
 *        If consider, then refer Amplitude Damping in Nielsen and Chuang Page 382, Equation (8.116) to (8.120)
 * 
 * @param amp_decay_prob amplitude damping probability
 * @param phase_decay_prob phase damping probiability
 * @param ideal_gate the 1-qubit ideal gate in superoperator form
 * @param sp  where superoperator is saved, either the noise superoperator itself, or combining with the ideal gate
 * @param add_to_gate if true, then sp will be saved by a superoperator for the noisy gate; if false, then only noise superoperator will be saved to sp.
 */
void addAmpPhaseErr1Q(double amp_decay_prob, double phase_decay_prob,
                 std::complex<double> ideal_gate[2][2], 
                 std::complex<double> sp[4][4],
                 bool add_to_gate = true)  {
    if ( amp_decay_prob < 0.0 || amp_decay_prob > 1.0 || phase_decay_prob < 0.0 || phase_decay_prob > 1.0) { 
        throw std::invalid_argument( "0 <= Decay probability <= 1. Please check the parameter." );
    }     
    // Compute super operator of ideal gate
    std::complex<double> ideal_gate_conj[2][2] = {}; // conjugate of ideal gate
    conjbar(ideal_gate[0], ideal_gate_conj[0], 2);
    std::complex<double> ideal_gate_sp[4][4] = {}; // super operator of ideal gate
    kronProd(ideal_gate_conj[0], ideal_gate[0], 1.0, ideal_gate_sp[0], 2);
    // If no noise, skip the computation of noise superoperator
    if (amp_decay_prob == 0.0 && phase_decay_prob == 0.0 && add_to_gate) {
        std::copy(&ideal_gate_sp[0][0], &ideal_gate_sp[0][0]+16,&sp[0][0]); // could use more efficient method
        return;
    } else if (amp_decay_prob == 0.0 && phase_decay_prob == 0.0) { // identity
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i == j) {
                    sp[i][j] = 1.0;
                } else {
                    sp[i][j] = 0.0;
                }
            }
        }
        return;
    }
    // Otherwise
    // Amplitude damping, excited state population=0
    std::complex<double> ampPhMats[3][2][2] = {
        {{{1,0},{0,0}},
        {{0,0},{std::sqrt(1-amp_decay_prob-phase_decay_prob),0}}},
        {{{0,0},{std::sqrt(amp_decay_prob),0}},
        {{0,0},{0,0}}},
        {{{0,0},{0,0}},
        {{0,0},{std::sqrt(phase_decay_prob),0}}},
    };
    // Transfer to a Superoperator
    std::complex<double> noise_sp[4][4] = {};
    kr2sp1Q(ampPhMats, noise_sp, 3);
    // compute final superoperator
    if (add_to_gate) {
        dotProd(noise_sp[0], ideal_gate_sp[0], sp[0], 4); // dot product the noise superoperator and ideal-gate superoperator
    } else { // if just need superoperator of noise
        std::copy(&noise_sp[0][0], &noise_sp[0][0]+16, &sp[0][0]);
    }
}

// The 2-qubit version of previous function.
// Essentially the Kronecker product of two 1-qubit amplitude and phase damping noise superoperator.
void addAmpPhaseErr2Q(double amp_decay_prob_1, double phase_decay_prob_1,
                      double amp_decay_prob_2, double phase_decay_prob_2,
                      std::complex<double> ideal_gate[4][4], 
                      std::complex<double> sp[16][16],
                      bool add_to_gate = true) {
    
    // Compute 1-qubit thermal relaxation noise superoperators
    std::complex<double> dummy_ideal_gate[2][2] = {{{1,0},{0,0}},{{0,0},{1,0}}}; // whatever that is should not matter, not used
    std::complex<double> gate1_sp[4][4] = {};
    std::complex<double> gate2_sp[4][4] = {};
    addAmpPhaseErr1Q(amp_decay_prob_1, phase_decay_prob_1, dummy_ideal_gate, gate1_sp, false);
    addAmpPhaseErr1Q(amp_decay_prob_2, phase_decay_prob_2, dummy_ideal_gate, gate2_sp, false);
    // Compute 2-qubit noise operator
    std::complex<double> noise_sp[16][16] = {};
    SuOpskronProd1q(gate2_sp, gate1_sp, noise_sp); // "Tensor Product (with reshuffle)" of two super operators
    if (!add_to_gate) { // if just need superoperator of noise
        std::copy(&noise_sp[0][0], &noise_sp[0][0]+256, &sp[0][0]);
        return;
    }
    // Compute super operator of ideal gate
    std::complex<double> ideal_gate_conj[4][4] = {}; // conjugate of ideal gate
    conjbar(ideal_gate[0], ideal_gate_conj[0], 4); 
    std::complex<double> ideal_gate_sp[16][16] = {}; // super operator of ideal gate
    kronProd(ideal_gate_conj[0], ideal_gate[0], 1.0, ideal_gate_sp[0], 4);
    // compute final superoperator
    dotProd(noise_sp[0], ideal_gate_sp[0], sp[0], 16); // dot product the noise superoperator and ideal-gate superoperator
}

//---------------------------------------  Single qubit Measurement Noise -------------------------------------
/**
 * @brief  Krasus operators for measurement error given bit-flip probability p_flip (1 qubit)
 *          i.e., equivalent to add a new gate before the measurement
 *          If Pr(0|1) = Pr(1|0), it can be simplified to the following channel
 *              D(\rho) = (1 - p_flip)\rho + p_flip X \rho X^\dagger
 *          However, they are usually different. We just need to change the value of non-zero entries of above channel a little bit
 * 
 * @param p_1g0 Pr(1|0)
 * @param p_0g1 Pr(0|1)
 * @param sp  where superoperator is saved.
 * @param rd_len readout length, not used by default. If used, please also specify T1 and T2
 * @param T1 T1, not used by default
 * @param T2 T2, not used by default
 */
void addMeaErr1Q(double p_1g0, double p_0g1, std::complex<double> sp[4][4],
                 double rd_len = 0.0, double T1 = INFINITY, double T2 = INFINITY) {

    if ( p_0g1 < 0.0 || p_0g1 > 1.0 || p_1g0 < 0.0 || p_1g0 > 1.0) { 
        throw std::invalid_argument( "Bit-flip probability should between 0 to 1." );
    }     
    if ( rd_len < 0.0 || T1 < 0.0 || T2 < 0.0) { 
        throw std::invalid_argument( "Time parameter(s) should be non-negative." );
    }     
    // If no error
    if (p_0g1 == 0.0 && p_1g0 == 0.0) { // identity
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i == j) {
                    sp[i][j] = 1.0;
                } else {
                    sp[i][j] = 0.0;
                }
            }
        }
        return;
    }
    // Initilize to zero matrix
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            sp[i][j] = 0;
        }
    }
    // Fill in non-zeros
    double coef_K1 = std::sqrt(1 - p_1g0) * std::sqrt(1 - p_0g1);
    double coef_K2 = std::sqrt(p_1g0) * std::sqrt(p_0g1);
    sp[0][0] = 1-p_1g0;
    sp[1][1] = coef_K1;
    sp[2][2] = coef_K1;
    sp[3][3] = 1-p_0g1;

    sp[0][3] = p_0g1;
    sp[1][2] = coef_K2;
    sp[2][1] = coef_K2;
    sp[3][0] = p_1g0;
    // If thermal relaxation noise is not counted in p_0g1
    if (rd_len > 0.0) {
        std::complex<double> ideal_gate[2][2] = {}; // dummy parameter
        std::complex<double> tr_sp[4][4] = {};

        addTRErr1Q(rd_len, T1, T2, ideal_gate, tr_sp, false);

        // compute final superoperator
        dotProd(sp[0], tr_sp[0], sp[0], 4); // dot product the noise superoperator and ideal-gate superoperator
    }
}




} //end of namespace DM-Sim
#endif
