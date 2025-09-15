

import numpy as np
from typing import List, Tuple

def gf2_rank(A: np.ndarray) -> int:
    A = (A % 2).astype(np.uint8).copy()
    m, n = A.shape
    r = 0
    for c in range(n):
        pivot = None
        for i in range(r, m):
            if A[i, c]:
                pivot = i
                break
        if pivot is None:
            continue
        if pivot != r:
            A[[r, pivot]] = A[[pivot, r]]
        for i in range(m):
            if i != r and A[i, c]:
                A[i, :] ^= A[r, :]
        r += 1
    return r

def gf2_inv(A: np.ndarray) -> np.ndarray:
    A = (A % 2).astype(np.uint8)
    n = A.shape[0]
    aug = np.concatenate([A.copy(), np.eye(n, dtype=np.uint8)], axis=1)
    r = 0
    for c in range(n):
        pivot = None
        for i in range(r, n):
            if aug[i, c]:
                pivot = i
                break
        if pivot is None:
            raise ValueError("Matrix is singular over GF(2).")
        if pivot != r:
            aug[[r, pivot]] = aug[[pivot, r]]
        for i in range(n):
            if i != r and aug[i, c]:
                aug[i, :] ^= aug[r, :]
        r += 1
    return aug[:, n:] % 2

def synthesize_cnot_from_A(A: np.ndarray) -> List[Tuple[int,int]]:
    """
    Given an invertible n x n binary matrix A over GF(2) representing the images of X_i
    (rows = transformed X_i as products of original X_j), produce a list of CNOTs (control, target)
    in Stim format (CX control target) such that applying them (in returned order) maps the standard
    X basis to A.
    
    Strategy: perform Gaussian elimination on a working copy of A to reduce it to the identity, using
    only row-xor operations (R_i ^= R_j), and record each such operation as a CNOT (i j).
    These operations map A -> I. Reversing the list gives I -> A (the desired forward circuit).
    Row swaps are implemented as three XORs (i ^= j; j ^= i; i ^= j) — i.e., three CNOTs.
    """
    A = (A % 2).astype(np.uint8).copy()
    n = A.shape[0]
    ops = []  # operations that reduce A to I (each is a tuple (i, j) meaning R_i ^= R_j)

    for col in range(n):
        # Find a pivot in or below row=col
        pivot = None
        for r in range(col, n):
            if A[r, col]:
                pivot = r
                break
        if pivot is None:
            raise ValueError(f"No pivot found in column {col}; A must be invertible.")
        # If pivot is not current row, swap rows using 3 XORs
        if pivot != col:
            i, j = col, pivot
            # swap rows i <-> j using 3 XORs
            A[i, :] ^= A[j, :]; ops.append((i, j))
            A[j, :] ^= A[i, :]; ops.append((j, i))
            A[i, :] ^= A[j, :]; ops.append((i, j))
        # Now eliminate all other 1s in this column
        for r in range(n):
            if r != col and A[r, col]:
                A[r, :] ^= A[col, :]
                ops.append((r, col))

    # At this point, A is identity.
    # The ops we recorded map A -> I. Reverse them to get I -> A.
    ops_forward = list(reversed(ops))
    return ops_forward

def verify_AX_B_relation(A: np.ndarray, B: np.ndarray) -> bool:
    invA_T = gf2_inv(A).T % 2
    return np.array_equal(invA_T, B % 2)

# (reuse gf2_rank, gf2_inv, synthesize_cnot_from_A, stim_from_ops from earlier)

def map_ops_to_physical(
    ops: List[Tuple[int,int]],
    inject_indices: List[int],
) -> List[Tuple[int,int]]:
    """
    Map operations (control, target) in matrix-index space 0..n-1
    into physical indices using inject_indices[i] -> physical qubit id.

    Assumes the ordering of inject_indices aligns with the rows/cols of A.
    """
    n = len(inject_indices)
    # basic sanity checks
    assert len(set(inject_indices)) == n, "inject_indices must be unique."
    assert all(0 <= c < n and 0 <= t < n for c, t in ops), "op index out of range."
    return [(inject_indices[c], inject_indices[t]) for (c, t) in ops]


def synthesize_and_map(
    X_support: np.ndarray,
    Z_support: np.ndarray = None,
    inject_indices: List[int] = None,
    verify_relation: bool = True,
):
    """
    Full pipeline:
      - optional verify: Z == (X^{-1}).T  (no X/Z mixing, CNOT-only)
      - synthesize forward CNOTs for X_support
      - map to physical indices using inject_indices
      - also produce inverse (reverse order)

    Returns:
      (ops_forward, ops_inverse, stim_forward_phys, stim_inverse_phys)
    """
    X = (X_support % 2).astype(np.uint8)
    n = X.shape[0]

    if X.shape[0] != X.shape[1]:
        raise ValueError("X_support must be square.")
    if inject_indices is None or len(inject_indices) != n:
        raise ValueError("inject_indices must be provided and match X dimension.")

    if verify_relation and Z_support is not None:
        invXT = gf2_inv(X).T % 2
        if not np.array_equal(invXT, (Z_support % 2)):
            raise ValueError("Z_support != (X_support^{-1}).T; not CNOT-only (would need X/Z-mixing).")

    # Synthesize in matrix index space
    ops_forward = synthesize_cnot_from_A(X)
    ops_inverse = list(reversed(ops_forward))  # CX is self-inverse

    # Map to physical indices
    ops_forward_phys = map_ops_to_physical(ops_forward, inject_indices)
    ops_inverse_phys = map_ops_to_physical(ops_inverse, inject_indices)

    return ops_inverse_phys

# --- Example usage (pseudo-context) ---
# inject_data = ...  # from your code
# coord_to_index = ...  # from bb_code
# inject_indices = [coord_to_index[q] for q in inject_data]  # your line

# X_support = ...  # n x n binary numpy array
# Z_support = ...  # (optional) n x n binary numpy array

# ops_f_phys, ops_i_phys, stim_f_phys, stim_i_phys = synthesize_and_map(
#     X_support, Z_support, inject_indices, verify_relation=True
# )

# with open("forward_physical.stim", "w") as f:
#     f.write(stim_f_phys + "\n")
# with open("inverse_physical.stim", "w") as f:
#     f.write(stim_i_phys + "\n")