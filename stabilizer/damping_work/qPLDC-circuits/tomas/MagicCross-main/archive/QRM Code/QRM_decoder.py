import numpy as np

def repetition_ml_decoder(y: np.ndarray) -> np.ndarray:
    """
    RM(0, m): all-0 or all-1 codeword; choose the better one
    """
    factor = 0.99 # Avoid extreme values of y (+1 or -1)
    score0 = np.prod((1 + factor * y) / 2)   # likelihood of all 0 codeword (+1)
    score1 = np.prod((1 - factor * y) / 2)   # likelihood of all 1 codeword (-1)

    c0 = np.ones_like(y)
    c1 = -np.ones_like(y)

    # print(f"Repetition ML scores: {score0:.4f} (0), {score1:.4f} (1)")
    if score0 > score1:
        return c0
    elif score1 > score0:
        return c1
    else:
        # If scores are equal, return with equal probability, unbiased
        return c0 if np.random.rand() < 0.5 else c1

def rm_list_decode(y, r, m, L=4) -> np.ndarray:
    """
    Recursive list decoder Ψ_m^r(L)
    y: input vector in [-1, 1]^n
    r, m: RM(r, m)
    L: list size
    Returns: list of (codeword, path_cost)
    """
    # print("r = ", r, ", m = ", m)
    factor = 0.99  # Avoid extreme values of y (+1 or -1)
    N = len(y)

    # Base cases
    if r == 0:
        return repetition_ml_decoder(y)
    if r == m:
        # RM(m,m) is full space — just use hard-decision
        return np.sign(y + (1-factor) * np.random.uniform(-1, 1, size=N))  # Add noise to avoid extreme values

    yL, yR = y[:N//2], y[N//2:]

    # Step 1: estimate v from componentwise product y_v = yL * yR
    y_v = yL * yR
    v_hat = rm_list_decode(y_v, r-1, m-1, L)

    # Step 2: for each v_hat, estimate u_hat
    y_hat = yR * v_hat  # Compute ŷ = yR * v_hat
    # Compute y_u = (yL + ŷ) / (1 + yL * ŷ)
    denom = 1 + factor * yL * y_hat
    y_u = (yL + y_hat) / denom
    u_hat = rm_list_decode(y_u, r, m-1, L)

    cw = np.concatenate([u_hat, u_hat * v_hat])

    return cw

def RM_decoder(r, m, noisy_cws): # punctured RM code decoder
    decoded_cws = []
    # Insert the 0-th bit (punctured)
    for noisy_cw in noisy_cws:
        y = (-1) ** noisy_cw
        y0 = np.insert(y, 0, 1)  # Insert a 1 at the 0-th bit
        y1 = np.insert(y, 0, -1)  # Insert a -1 at the 0-th bit
        decoded_cw0 = rm_list_decode(y0, r, m, L=4)
        decoded_cw1 = rm_list_decode(y1, r, m, L=4)
        overlap0 = np.sum(y0 * decoded_cw0)
        overlap1 = np.sum(y1 * decoded_cw1)
        if overlap0 > overlap1:
            decoded_cw = decoded_cw0
        else:
            decoded_cw = decoded_cw1
        decoded_cw = (1 - decoded_cw[1:]) // 2  # Convert from ±1 to {0, 1}
        decoded_cws.append(decoded_cw)
    return np.array(decoded_cws, dtype=int)