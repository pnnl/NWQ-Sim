import numpy as np

def sparsify_binary_matrix(A, max_iter=10_000):
    """
    Greedily sparsify a binary matrix A (entries in {0,1}) over F2.
    By repeatedly replacing row i with row_i XOR row_j whenever that
    strictly reduces the number of 1's in row i.

    Parametersx
    ----------
    A : ndarray of shape (m, n), dtype {0,1}
        Input binary matrix.
    max_iter : int
        Maximum number of row-add operations to try before giving up.

    Returns
    -------
    A_sparse : ndarray of shape (m, n), dtype {0,1}
        A permuted copy of A with (heuristically) maximized number of 0's.
    ops : list of tuples
        List of operations performed, each of the form (i, j),
        meaning row i := row_i XOR row_j was applied.
    """
    A = A.copy().astype(np.uint8)
    m, n = A.shape

    def ones_count(row):
        return row.sum()

    ops = []
    for it in range(max_iter):
        best_reduction = 0
        best_i = best_j = None

        # search for the best row-pair to reduce ones
        for i in range(m):
            cnt_i = ones_count(A[i])
            # skip rows already zero
            if cnt_i == 0:
                continue
            for j in range(m):
                if i == j:
                    continue
                # XOR i with j
                new_row = A[i] ^ A[j]
                reduction = cnt_i - ones_count(new_row)
                if reduction > best_reduction:
                    best_reduction = reduction
                    best_i, best_j = i, j

        # if no improvement, stop
        if best_reduction <= 0:
            break

        # apply the best operation
        A[best_i] ^= A[best_j]
        ops.append((best_i, best_j))

    # finally, permute rows in non-decreasing order of their weight
    weights = A.sum(axis=1)
    perm = np.argsort(weights)
    A_sparse = A[perm]

    return A_sparse, ops