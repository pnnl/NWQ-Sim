import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
 
def bloch_measurement_operators():
    """Return Pauli measurement POVMs for 0/1 outcome along X/Y/Z"""
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])
 
    def meas_ops(pauli):
        eigvals, eigvecs = np.linalg.eigh(pauli)
        return [np.outer(eigvecs[:,i], eigvecs[:,i].conj()) for i in range(2)]
 
    return {
        'X': meas_ops(X),
        'Y': meas_ops(Y),
        'Z': meas_ops(Z),
    }
 
def param_to_rho(params):
    """Map 4 real parameters to 2x2 physical density matrix via Cholesky decomposition"""
    a, b, c, d = params
    T = np.array([[a, 0], [c + 1j * d, b]])
    rho = T.conj().T @ T
    trace = np.trace(rho)
    if np.abs(trace) < 1e-8 or np.isnan(trace):
        raise ValueError("Optimizer failed, trace is approaching 0.")
    
    return rho / trace
 
def negative_log_likelihood(params, measurements):
    """Compute negative log likelihood for given parameters.
    Essentially, the routine to minimize the density matrix with
    feedback from the tomography routine."""

    rho = param_to_rho(params)
    povms = bloch_measurement_operators()
    log_likelihood = 0
    for axis, (n0, n1) in measurements.items():
        E0, E1 = povms[axis.upper()]
        epsilon = 1e-8
        p0 = np.clip(np.real(np.trace(rho @ E0)), epsilon, 1 - epsilon)
        p1 = np.clip(np.real(np.trace(rho @ E1)), epsilon, 1 - epsilon)
        log_likelihood += n0 * np.log(p0) + n1 * np.log(p1)

    if np.isnan(p0) or np.isnan(p1):
        raise ValueError("NaN detected in probs", params)

    return -log_likelihood

def plot_density_matrix(rho):
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    n = rho.shape[0]
    real_part = np.real(rho)
    imag_part = np.imag(rho)

    # Plot real part
    im1 = axs[0].imshow(real_part, cmap='coolwarm', vmin=-1, vmax=1)
    axs[0].set_title(r'$\Re(\rho)$')
    if(n < 3):
        axs[0].set_xticks([0, 1])
        axs[0].set_yticks([0, 1])
    elif(n == 4):
        axs[0].set_xticks([0, 1, 2, 3])
        axs[0].set_yticks([0, 1, 2, 3])
        # axs[0].set_xticklabels(['0', '1', '+', 'i'])
        # axs[0].set_yticklabels(['0', '1', '+', 'i'])
    for i in range(n):
        for j in range(n):
            axs[0].text(j, i, f"{real_part[i, j]:.2f}", ha='center', va='center', color='black')

    # Plot imaginary part
    im2 = axs[1].imshow(imag_part, cmap='coolwarm', vmin=-1, vmax=1)
    axs[1].set_title(r'$\Im(\rho)$')
    if(n < 3):
        axs[1].set_xticks([0, 1])
        axs[1].set_yticks([0, 1])
    elif(n == 4):
        axs[1].set_xticks([0, 1, 2, 3])
        axs[1].set_yticks([0, 1, 2, 3])
        # axs[1].set_xticklabels(['0', '1', '+', 'i'])
        # axs[1].set_yticklabels(['0', '1', '+', 'i'])
    for i in range(n):
        for j in range(n):
            axs[1].text(j, i, f"{imag_part[i, j]:.2f}", ha='center', va='center', color='black')

    fig.colorbar(im1, ax=axs[0])
    fig.colorbar(im2, ax=axs[1])
    plt.tight_layout()
    plt.show()

 
def state_tomography_mle(measurements):
    init_params = [1, 0.5, 0, 0]
    result = minimize(
        negative_log_likelihood,
        init_params,
        args=(measurements,),
        method='L-BFGS-B',  # more stable than BFGS
        options={'maxiter': 1000}
    )

    if not result.success:
        #Accept near-optimal result if it's close to converging
        if result.status == 2 and result.fun < 1e6:
            print("Warning: Optimization did not fully converge, but result deemed close enough.")
        else:
            raise RuntimeError("MLE optimization failed: " + result.message)

    return param_to_rho(result.x)

GAMMA_MTX = np.array([
        [1, 0, 0, 1],
        [0, 1, 1, 0],
        [0, 1, -1, 0],
        [1, 0, 0, -1]
    ]) / 2

def process_tomography(measurements):
    """
    This process tomography function is based on Nielsen and Chuang's book, Page 393, Box 8.5
    """
    '''
    Note:
    the measurement should be aligned in the roder of:
    [counts for |0> state, 
    counts for |1> state,
    counts for |+> state,
    counts for |+y> state]
    '''
    # if measurements.shape != (4, 6):
    #     try:
    #         measurements = np.array(measurements).reshape(4, 6)
    #     except ValueError:
    #         raise ValueError("Measurements must be a 4x6 array or list of lists.")
        
    dm_list = []
    chi_matrix = np.zeros((4, 4), dtype=np.complex128)
    for m in measurements:
        dm_list.append(state_tomography_mle(m))

    dm2 = dm_list[2] - 1j * dm_list[3] - (1 - 1j)/2 * (dm_list[0] + dm_list[1])
    dm3 = dm_list[2] + 1j * dm_list[3] - (1 + 1j)/2 * (dm_list[0] + dm_list[1])
    
    dms = np.bmat([[dm_list[0], dm3],
                    [dm2, dm_list[1]]])
    
    chi_matrix = GAMMA_MTX @ dms @ GAMMA_MTX
    return chi_matrix

