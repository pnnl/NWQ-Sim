import numpy as np
from scipy.optimize import minimize

def single_qubit_state_tomography(measurement_counts):
    """
    Reconstruct a single-qubit density matrix from measurement outcomes along X, Y, and Z bases.
    Basic state tomography function.

    Args:
        measurement_counts: List of 3 lists, each containing two integers.
            Format: [[count_x_0, count_x_1], [count_y_0, count_y_1], [count_z_0, count_z_1]]

    Returns:
        rho: A 2x2 complex numpy array representing the reconstructed density matrix.
    """
    def compute_expectation(count_0, count_1):
        total = count_0 + count_1
        if total == 0:
            return 0.0
        return (count_0 - count_1) / total

    # Extract measurement counts
    (cx0, cx1), (cy0, cy1), (cz0, cz1) = measurement_counts

    # Estimate expectation values
    rx = compute_expectation(cx0, cx1)
    ry = compute_expectation(cy0, cy1)
    rz = compute_expectation(cz0, cz1)

    # Form Bloch vector
    bloch_vector = np.array([rx, ry, rz])
    norm = np.linalg.norm(bloch_vector)

    # Project back into Bloch sphere if necessary
    if norm > 1.0:
        bloch_vector = bloch_vector / norm  # Normalize to unit vector
        bloch_vector *= 0.999999  # Slightly shrink to stay strictly within Bloch ball

    # Construct density matrix: ρ = (I + r⋅σ) / 2
    I = np.eye(2)
    sigma_x = np.array([[0, 1], [1, 0]])
    sigma_y = np.array([[0, -1j], [1j, 0]])
    sigma_z = np.array([[1, 0], [0, -1]])
    rho = 0.5 * (I +
                 bloch_vector[0] * sigma_x +
                 bloch_vector[1] * sigma_y +
                 bloch_vector[2] * sigma_z)

    return rho

########################################
#
# MLE-based state tomography function
#
########################################

class MLE_SQStateTomography:
    BMO = np.array([
        np.array([[1, 1],[1, 1]])/2, # +X
        np.array([[1, -1],[1, -1]])/2, # -X
        np.array([[1, -1j],[1j, 1]])/2, # +Y
        np.array([[1, 1j],[-1j, 1]])/2, # -Y
        np.array([[1, 0],[0, 0]]), # +Z
        np.array([[0, 0],[0, 1]]) # -Z
        ])


    def __init__(self, measurements, init_params=None):
        """
        Initialize the MLE-based single-qubit state tomography with measurement data.

        Args:
        measurements (list of list of int): A list of three [n_0, n_1] pairs
            representing measurement counts in the X, Y, and Z Pauli bases.
            Format:
                measurements = [[n_x0, n_x1], [n_y0, n_y1], [n_z0, n_z1]]
            - [n_x0, n_x1] are counts of outcomes in the X basis (|+⟩, |−⟩)
            - [n_y0, n_y1] are counts in the Y basis (|+i⟩, |−i⟩)
            - [n_z0, n_z1] are counts in the Z basis (|0⟩, |1⟩)

            Each entry n_bk must be a non-negative integer corresponding to
            how many times outcome k was observed in basis b.
        """
        self.measurements = measurements
        if init_params is not None:
            if len(init_params) != 4:
                raise ValueError("Initial parameters must be a list of 4 real numbers.")
            self.params = init_params
        else:
            self.params = [1, 1, 0, 0]  # This gives rho = I / 2
        
        self.density_matrix = self.param_to_rho(self.params)

        return
        
    @staticmethod
    def param_to_rho(params):
        """Map 4 real parameters to 2x2 physical density matrix via Cholesky decomposition"""
        a, b, c, d = params
        T = np.array([[a, 0], [c + 1j * d, b]])
        rho = T.conj().T @ T
        return rho / np.trace(rho)

    def _negative_log_likelihood(self, params, measurements):
        """
        Compute the negative log-likelihood function for a single-qubit quantum state,
        used in Maximum Likelihood Estimation (MLE) tomography.

        This function assumes the quantum state (density matrix) is parameterized via
        a Cholesky decomposition:
            ρ = T† T / Tr(T† T)
        where T is a lower-triangular complex matrix with real parameters.

        Args:
            params (list or array-like of float): A list of 4 real numbers [a, b, c, d]
                which define the lower-triangular matrix T as follows:
                    T = [[a,       0      ],
                         [c + i*d, b      ]]
                The corresponding density matrix is then constructed as:
                    ρ = T† T / Tr(T† T)
                This ensures that ρ is Hermitian, positive semidefinite, and trace-1.

            measurements (list of list of int): A list of three [n_0, n_1] pairs
                representing measurement counts in the X, Y, and Z Pauli bases.
                Format:
                    measurements = [[n_x0, n_x1], [n_y0, n_y1], [n_z0, n_z1]]
                - [n_x0, n_x1] are counts of outcomes in the X basis (|+⟩, |−⟩)
                - [n_y0, n_y1] are counts in the Y basis (|+i⟩, |−i⟩)
                - [n_z0, n_z1] are counts in the Z basis (|0⟩, |1⟩)

                Each entry n_bk must be a non-negative integer corresponding to
                how many times outcome k was observed in basis b.

        Returns:
            float: The negative log-likelihood value corresponding to the input parameters.
                This function is meant to be minimized by a numerical optimizer such as
                scipy.optimize.minimize. The resulting optimal parameters produce the
                physical density matrix that best fits the measurement data under the
                maximum likelihood principle.

        Notes:
            - The computed probabilities are clipped from below to avoid log(0),
              which ensures numerical stability.
            - The output density matrix ρ (via param_to_rho(params)) is always physical:
              Hermitian, positive semidefinite, and trace-1.
        """
        rho = self.param_to_rho(params)

        prob_array = np.einsum('ij, kji -> k', rho, self.BMO)
        measurements = np.array(measurements).flatten()

        prob_array = np.clip(prob_array, 1e-10, None)
        log_likelihood = np.dot(measurements, np.log(prob_array))
        
        return -log_likelihood

    def update_measurements(self, new_measurements):
        self.measurements = new_measurements
        return

    def state_tomography_mle(self, measurements=None):
        """
        Perform single-qubit state tomography using MLE.
        Input: measurements = [[nx0, nx1], [ny0, ny1], [nz0, nz1]]
        Output: 2x2 density matrix
        """
        # Initial guess: identity state
        if measurements is None:
            measurements = self.measurements

        result = minimize(self._negative_log_likelihood, self.params,
                          args=(measurements,), method='l-BFGS-B', 
                          options={"maxiter":1000, "gtol": 1e-6, "disp": True},
                          )

        if not result.success:
            raise RuntimeError("MLE optimization failed: " + result.message)

        self.density_matrix = self.param_to_rho(result.x)
        return self.density_matrix


class MLE_SQProcessTomography:
    GAMMA_MTX = np.array([
        [1, 0, 0, 1],
        [0, 1, 1, 0],
        [0, 1, -1, 0],
        [1, 0, 0, -1]
    ]) / 2
    def __init__(self, measurements):
        '''
        Note:
        the measurement should be aligned in the roder of:
        [counts for |0> state, 
        counts for |1> state,
        counts for |+> state,
        counts for |+y> state]
        '''
        self.measurements = np.array(measurements)
        if self.measurements.shape != (4, 6):
            try:
                self.measurements = np.array(measurements).reshape(4, 6)
            except ValueError:
                raise ValueError("Measurements must be a 4x6 array or list of lists.")
        
        self.chi_matrix = np.zeros((4, 4), dtype=np.complex128)
        return
    
    def process_tomography(self):
        """
        This process tomography function is based on Nielsen and Chuang's book, Page 393, Box 8.5
        """
        dm_list = []
        for m in self.measurements:
            temp = MLE_SQStateTomography(m)
            dm_list.append(temp.state_tomography_mle())

        dm2 = dm_list[2] - 1j * dm_list[3] - (1 - 1j)/2 * (dm_list[0] + dm_list[1])
        dm3 = dm_list[2] + 1j * dm_list[3] - (1 + 1j)/2 * (dm_list[0] + dm_list[1])
        
        dms = np.bmat([[dm_list[0], dm2],
                       [dm3, dm_list[1]]])
        
        self.chi_matrix = self.GAMMA_MTX @ dms @ self.GAMMA_MTX
        return self.chi_matrix



    