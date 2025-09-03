import numpy as np
import pulp
import matrixsparse

def solve_block_structure(A, k, max_p=None, time_limit=None):
    """
    Solve for the optimal row/column permutation into a p x p anti-diagonal block structure,
    respecting that the first k rows of A may only go to U or the first row-partition of L.

    Parameters:
    - A: binary 2D numpy array of shape (m, n)
    - k: integer count of rows from original matrix that are restricted to U or row-partition 1
    - max_p: maximum p to try (defaults to min(m, n))
    - time_limit: optional solver time limit in seconds

    Returns:
    - best: dict with keys 'p', 'obj', 'row_group', 'col_group'
    """
    m, n = A.shape
    if max_p is None:
        max_p = min(m, n)

    best = {'p': None, 'obj': float('inf'), 'row_group': None, 'col_group': None}

    # Loop over possible block counts p
    for p in range(1, max_p + 1):
        prob = pulp.LpProblem(f"block_struct_p{p}", pulp.LpMinimize)

        # Variables
        x = pulp.LpVariable.dicts('x', (range(m), range(p+1)), cat='Binary')  # row assignments
        y = pulp.LpVariable.dicts('y', (range(n), range(1, p+1)), cat='Binary')  # column assignments
        M = pulp.LpVariable('M', lowBound=0, cat='Integer')  # max anti-diagonal block height

        # Objective: |U| + M
        prob += pulp.lpSum(x[a][0] for a in range(m)) + M

        # Each row in exactly one group
        for a in range(m):
            prob += pulp.lpSum(x[a][r] for r in range(p+1)) == 1
        # Restrict first k rows
        for a in range(min(k, m)):
            for r in range(2, p+1):
                prob += x[a][r] == 0
        # Each column in exactly one partition
        for b in range(n):
            prob += pulp.lpSum(y[b][c] for c in range(1, p+1)) == 1
        # M bounds
        for r in range(1, p+1):
            prob += M >= pulp.lpSum(x[a][r] for a in range(m))
        # Anti-diagonal constraints
        for a in range(m):
            for b in range(n):
                if A[a, b] == 1:
                    for r in range(1, p+1):
                        for c in range(1, p+1):
                            if r + c != p + 1:
                                prob += x[a][r] + y[b][c] <= x[a][0] + 1
        # Solve
        if time_limit:
            prob.solve(pulp.PULP_CBC_CMD(timeLimit=time_limit))
        else:
            prob.solve()

        # Check and update
        if prob.status == pulp.LpStatusOptimal:
            obj = pulp.value(prob.objective)
            if obj < best['obj']:
                best['p'] = p
                best['obj'] = obj
                # extract groups
                row_group = np.zeros(m, dtype=int)
                for a in range(m):
                    for r in range(p+1):
                        if pulp.value(x[a][r]) > 0.5:
                            row_group[a] = r
                col_group = np.zeros(n, dtype=int)
                for b in range(n):
                    for c in range(1, p+1):
                        if pulp.value(y[b][c]) > 0.5:
                            col_group[b] = c
                best['row_group'] = row_group
                best['col_group'] = col_group
    return best


def get_permuted_matrix(A, row_group, col_group):
    """
    Permute matrix A according to row_group and col_group.

    Returns the permuted matrix P and the row/column split indices.

    - Rows: U (group 0) first, then groups 1..p.
    - Columns: groups 1..p in order.
    """
    m, n = A.shape
    p = max(col_group)
    ordered_rows = np.hstack([np.where(row_group == r)[0] for r in range(p+1)])
    ordered_cols = np.hstack([np.where(col_group == c)[0] for c in range(1, p+1)])
    P = A[ordered_rows][:, ordered_cols]

    row_counts = [np.sum(row_group == r) for r in range(p+1)]
    col_counts = [np.sum(col_group == c) for c in range(1, p+1)]
    row_splits = np.cumsum(row_counts)  # indices after each group
    col_splits = np.cumsum(col_counts)
    return P, row_splits, col_splits


if __name__ == "__main__":
    # Example usage
    tri49 = np.array([
    # even-weight rows (rows 1 – 13) – this is G0 in App. B of Bravyi-Haah
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

    a15to1 = np.array([
        [0, 0, 0, 0, 0, 1, 1, 1, 1, 1 ,1 ,1 ,0, 0, 0],
        [1, 0, 0, 0, 1, 1, 1 ,0 ,0 ,1, 0, 1, 0, 1, 1],
        [0, 1, 0, 0, 1, 1, 0 ,1 ,0 ,0, 1, 1, 1, 0, 1],
        [0, 0, 1, 0, 1, 0, 1 ,1 ,1 ,0, 0, 1, 1, 1, 0],
        [0, 0, 0, 1, 0, 0, 0 ,0 ,1 ,1, 1, 1, 1, 1, 1]
    ], dtype=int)

    a14to2 = np.array([
        [0, 0, 0, 0, 1, 1, 1, 1, 1 ,1 ,1 ,0, 0, 0],
        [0, 0, 0, 1, 1, 1 ,0 ,0 ,1, 0, 1, 0, 1, 1],
        [1, 0, 0, 1, 1, 0 ,1 ,0 ,0, 1, 1, 1, 0, 1],
        [0, 1, 0, 1, 0, 1 ,1 ,1 ,0, 0, 1, 1, 1, 0],
        [0, 0, 1, 0, 0, 0 ,0 ,1 ,1, 1, 1, 1, 1, 1]
    ], dtype=int)

    a20to4 = np.array([
        [0, 0, 0, 0, 1, 1, 1, 0, 1 ,0 ,1 ,0, 1, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0 ,0 ,1 ,1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0 ,0 ,0 ,1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0 ,0 ,0 ,1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1],
        [0, 0, 1, 1, 1, 0 ,1 ,1 ,0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0],
        [0, 1, 0, 1, 1, 1 ,0 ,1 ,1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1],
        [1, 0, 0, 1, 0, 1 ,1 ,0 ,0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1]
    ], dtype=int)

    hex_values = [
    "000000000000ffff",
    "000f000f000f000f",
    "1111111111111111",
    "000000ff000000ff",
    "0303030303030303",
    "0000000055555555",
    "ffffffffffffffff",
    "00000000ffffffff",
    "0000ffff0000ffff",
    "00ff00ff00ff00ff",
    "0f0f0f0f0f0f0f0f",
    "3333333333333333",
    "5555555555555555",
    "0000000000ff00ff",
    "0033003300330033",
    "0000000033333333",
    "0000000000330033",
    ]


    a64to2ccz = np.array([
    list(map(int, bin(int(h, 16))[2:].zfill(64)))
    for h in hex_values
    ], dtype=np.uint8)


    A = a64to2ccz
    k = 0

    A_sparse, permutation = matrixsparse.sparsify_binary_matrix(A, k)
   

    A=A_sparse

    result = solve_block_structure(A, k=k, max_p=3, time_limit=1000)
    P, row_splits, col_splits = get_permuted_matrix(A, result['row_group'], result['col_group'])

    output_path = 'compress/permuted_visualization_a15.txt'
    with open(output_path, 'w') as f:
        f.write(f"Best p: {result['p']}\n")
        f.write(f"Objective: {result['obj']}\n\n")
        f.write("Row groups (0=U,1..p): " + ' '.join(map(str, result['row_group'])) + "\n")
        f.write("Column partitions (1..p): " + ' '.join(map(str, result['col_group'])) + "\n\n")
        f.write(f"Permutation applied: {permutation}\n")
        f.write(f"Sparsified matrix: {A_sparse}\n")
        f.write("Permuted matrix P with splitters:\n")
        # Write matrix with column splitters '|' and row splitters '---'
        for i, row in enumerate(P):
            elements = []
            for j, val in enumerate(row):
                elements.append(str(int(val)))
                if j+1 in col_splits[:-1]:
                    elements.append('|')
            row_str = ' '.join(elements)
            f.write(row_str + "\n")
            if i+1 in row_splits[:-1]:
                f.write('-' * len(row_str) + "\n")
        f.write("\nRow split indices: " + ' '.join(map(str, row_splits)) + "\n")
        f.write("Column split indices: " + ' '.join(map(str, col_splits)) + "\n")
    print(f"Results written to {output_path}")

