import numpy as np
import itertools

def solve_matrix_partition_p2_strict(A: np.ndarray, k: int):
    """
    Finds an optimal permutation of a binary matrix A to fit a specific structure
    where the lower half has no "U" columns.

    The goal is to permute the rows and columns of A to form A' such that:
    A' = | U       |
         |---------|
         | 0   L12 |
         | L21 0   |
    
    This is achieved by partitioning rows into (U_rows, L1_rows, L2_rows) and ALL
    columns into (L1_cols, L2_cols). This is stricter than the previous version.
    A row partition is only valid if all columns can be partitioned this way.
    The cost function to minimize is:
    cost = len(U_rows) + max(len(L1_rows), len(L2_rows)).

    A constraint requires that the first k rows of the original matrix A can only end up
    in the U block or the L1 block.

    Args:
        A (np.ndarray): The input binary matrix (m x n).
        k (int): The number of initial rows with placement constraints.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: The optimally permuted matrix A'.
            - dict: A dictionary describing the row/column partitions of the optimal solution.
            - float: The minimum cost found.
        Returns (None, None, "Error message") if no solution is found or if inputs are invalid.
    """
    if not isinstance(A, np.ndarray) or A.ndim != 2:
        return None, None, "Error: Input matrix A must be a 2D numpy array."
    if not (0 <= k <= A.shape[0]):
        return None, None, f"Error: k must be between 0 and the number of rows ({A.shape[0]})."

    m, n = A.shape
    min_cost = float('inf')
    best_config = None

    # Iterate through all 3^m possible ways to partition the m rows
    # into three sets: U (group 0), L1 (group 1), and L2 (group 2).
    for partition in itertools.product(range(3), repeat=m):
        # 1. Assign rows to U, L1, and L2 based on the current partition
        U_rows = [r for r, group in enumerate(partition) if group == 0]
        L1_rows = [r for r, group in enumerate(partition) if group == 1]
        L2_rows = [r for r, group in enumerate(partition) if group == 2]

        # 2. Check the k-row constraint.
        if any(r < k for r in L2_rows):
            continue

        # 3. Determine column partitions.
        # Find all columns that are candidates for L1_cols (all zeros in L1_rows)
        L1_potential_cols = set(range(n))
        if L1_rows:
            is_zero_for_L1 = np.all(A[L1_rows, :] == 0, axis=0)
            L1_potential_cols = {j for j, is_zero in enumerate(is_zero_for_L1) if is_zero}

        # Find all columns that are candidates for L2_cols (all zeros in L2_rows)
        L2_potential_cols = set(range(n))
        if L2_rows:
            is_zero_for_L2 = np.all(A[L2_rows, :] == 0, axis=0)
            L2_potential_cols = {j for j, is_zero in enumerate(is_zero_for_L2) if is_zero}

        # 4. NEW (Stricter) Constraint Check:
        # The union of potential L1 and L2 columns MUST cover all columns.
        # If not, this row partition is invalid because some columns would have
        # non-zero entries in both L1 and L2 sections, which is not allowed.
        if L1_potential_cols.union(L2_potential_cols) != set(range(n)):
            continue

        # If the check passes, we can successfully partition all columns.
        # We assign columns that are candidates for L1 to L1_cols.
        # All remaining columns must be candidates for L2_cols.
        L1_cols = sorted(list(L1_potential_cols))
        L2_cols = sorted(list(set(range(n)) - L1_potential_cols))

        # 5. Calculate the cost for this specific partition.
        cost = len(U_rows) + max(len(L1_rows), len(L2_rows))

        # 6. If this partition is better than the best one found so far, save it.
        if cost < min_cost:
            min_cost = cost
            best_config = {
                'U_rows': sorted(U_rows), 'L1_rows': sorted(L1_rows), 'L2_rows': sorted(L2_rows),
                'L1_cols': L1_cols, 'L2_cols': L2_cols,
                'cost': cost
            }

    # 7. After checking all partitions, reconstruct the matrix from the best configuration.
    if best_config is None:
        return None, None, "No valid partition could be found that satisfies the stricter constraints."

    row_permutation = best_config['U_rows'] + best_config['L1_rows'] + best_config['L2_rows']
    col_permutation = best_config['L1_cols'] + best_config['L2_cols']

    A_prime = A[np.ix_(row_permutation, col_permutation)]

    return A_prime, best_config, min_cost

# --- Example Usage ---
if __name__ == '__main__':
    # Example 1: A simple matrix that can be perfectly partitioned.
    print("--- Example 1 ---")
    A1 = np.array([
    # even-weight rows (rows 1 – 13) – this is G₀ in App. B of Bravyi-Haah
    list(map(int, "1111111111111110101010101010101010101010101010101")),
    list(map(int, "0000000000000000000111100110011000011001100110011")),
    list(map(int, "0000000000000001100000011001100110000000000000000")),
    list(map(int, "0000000000000000000000000000000001111000000001111")),
    list(map(int, "0000000000000000011110000000000000000111100000000")),
    list(map(int, "0000000000000000000001111000011110000000000000000")),
    list(map(int, "0000000000000000000000000111111110000000000000000")),
    list(map(int, "0000000000000000000000000000000001111111100000000")),
    list(map(int, "0000000000000000000000000000000000000000011111111")),
    list(map(int, "1010101010101010000000000000000000000000000000000")),
    list(map(int, "0110011001100110000000000000000000000000000000000")),
    list(map(int, "0001111000011110000000000000000000000000000000000")),
    list(map(int, "0000000111111110000000000000000000000000000000000"))
 ], dtype=int)
    k1 = 1 # First 2 rows must be in U or L1
    
    A_prime1, config1, cost1 = solve_matrix_partition_p2_strict(A1, k1)

    if A_prime1 is not None:
        print("Original Matrix A1:")
        print(A1)
        print("\nOptimal Permuted Matrix A1':")
        print(A_prime1)
        print(f"\nAchieved Minimum Cost: {cost1}")
        print("\nPartition Configuration:")
        print(config1)

