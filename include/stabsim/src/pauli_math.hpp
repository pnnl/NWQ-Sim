#include "../../../Eigen/Dense"
#include "../../../Eigen/Sparse"

//Using dynamic Eigen MatrixXd to streamline calculations
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Vector2D = Eigen::Vector2d;
using ComplexMatrix2D = Eigen::Matrix2cd;
using ComplexMatrix = Eigen::MatrixXcd;
using ComplexVector = Eigen::VectorXcd;
using Complex = std::complex<double>;
using ComplexSparseMatrix = Eigen::SparseMatrix<Complex>;
using Triplet = Eigen::Triplet<Complex>;

//Functions to define Pauli matrices (sigma X,Y,Z,+,-)
ComplexMatrix pauliI() {
    ComplexMatrix sigmaI = ComplexMatrix::Zero(2, 2);
    sigmaI(0, 0) = std::complex<double>(1, 0);
    sigmaI(1, 1) = std::complex<double>(1, 0);
    return sigmaI;
}

ComplexMatrix pauliX() {
    ComplexMatrix sigmaX = ComplexMatrix::Zero(2, 2);
    sigmaX(0, 1) = std::complex<double>(1, 0);
    sigmaX(1, 0) = std::complex<double>(1, 0);
    return sigmaX;
}

ComplexMatrix pauliY() {
    ComplexMatrix sigmaY = ComplexMatrix::Zero(2, 2);
    sigmaY(0, 1) = std::complex<double>(0, -1);
    sigmaY(1, 0) = std::complex<double>(0, 1);
    return sigmaY;
}

ComplexMatrix pauliZ() {
    ComplexMatrix sigmaZ = ComplexMatrix::Zero(2, 2);
    sigmaZ(0, 0) = std::complex<double>(1, 0);
    sigmaZ(1, 1) = std::complex<double>(-1, 0);
    return sigmaZ;
}

//Create a sparse identity matrix of whatever type is necessary
ComplexMatrix createIdentity(int n) {
    return ComplexMatrix::Identity(n, n);
}

ComplexMatrix kroneckerProduct(const ComplexMatrix& A, const ComplexMatrix& B) {
    int aRows = A.rows(), aCols = A.cols();
    int bRows = B.rows(), bCols = B.cols();

    ComplexMatrix result(aRows * bRows, aCols * bCols);

    for (int i = 0; i < aRows; i++) {
        for (int j = 0; j < aCols; j++) {
            result.block(i * bRows, j * bCols, bRows, bCols) = A(i, j) * B;
        }
    }

    return result;
}

ComplexSparseMatrix sparsePauliI() {
    ComplexSparseMatrix sigmaI(2, 2);
    sigmaI.insert(0, 0) = std::complex<double>(1, 0);
    sigmaI.insert(1, 1) = std::complex<double>(1, 0);
    return sigmaI;
}

ComplexSparseMatrix sparsePauliX() {
    ComplexSparseMatrix sigmaX(2, 2);
    sigmaX.insert(0, 1) = std::complex<double>(1, 0);
    sigmaX.insert(1, 0) = std::complex<double>(1, 0);
    return sigmaX;
}

ComplexSparseMatrix sparsePauliY() {
    ComplexSparseMatrix sigmaY(2, 2);
    sigmaY.insert(0, 1) = std::complex<double>(0, -1);
    sigmaY.insert(1, 0) = std::complex<double>(0, 1);
    return sigmaY;
}

ComplexSparseMatrix sparsePauliZ() {
    ComplexSparseMatrix sigmaZ(2, 2);
    sigmaZ.insert(0, 0) = std::complex<double>(1, 0);
    sigmaZ.insert(1, 1) = std::complex<double>(-1, 0);
    return sigmaZ;
}

//Create a sparse identity matrix of whatever type is necessary
ComplexSparseMatrix createSparseIdentity(int n) {
    ComplexSparseMatrix I(n, n);
    I.setIdentity();  //Sets the diagonal elements to 1 (implicitly (1,0) in complex)
    return I;
}

ComplexSparseMatrix kroneckerProduct(const ComplexSparseMatrix& A, const ComplexSparseMatrix& B) {
    int aRows = A.rows(), aCols = A.cols();
    int bRows = B.rows(), bCols = B.cols();
    ComplexSparseMatrix result(aRows * bRows, aCols * bCols);

    std::vector<Triplet> triplets;
    triplets.reserve(A.nonZeros() * B.nonZeros());  //Pre-allocate memory for triplets to avoid reallocation

    //Iterate over all elements in A
    for (int j = 0; j < aCols; ++j) {
        for (ComplexSparseMatrix::InnerIterator itA(A, j); itA; ++itA) {
            //Iterate over all elements in B
            for (int k = 0; k < bRows; ++k) {
                for (int l = 0; l < bCols; ++l) {
                    if (B.coeff(k, l) != Complex(0.0, 0.0)) { //Check if B's element is non-zero to avoid unnecessary entries
                        Complex val = itA.value() * B.coeff(k, l);
                        int row_index = itA.row() * bRows + k;
                        int col_index = itA.col() * bCols + l;
                        triplets.emplace_back(row_index, col_index, val);
                    }
                }
            }
        }
    }

    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

//Computes the matrix exponential using a Taylor expansion e^{M}
//This feature is available in unsupported Eigen using a more accurate Padé approximant, it is explicitly implemented here for reliability
ComplexMatrix matrixExponential(const ComplexMatrix& A, int numTerms = 10) {
    ComplexMatrix result(A.rows(), A.cols());
    result.setIdentity(); //Initialize result to identity matrix

    ComplexMatrix term(A.rows(), A.cols());
    term.setIdentity(); //Initialize term to identity matrix

    for (int i = 1; i < numTerms; i++) {
        term = (term * A) / static_cast<double>(i);

        result += term;
        
    }
    return result;
}