import numpy as np
import pulp
from matrixsparse import sparsify_matrix


def gaussian_eliminate_partial(A, k, max_pivots):
    """
    Perform up to `max_pivots` Gaussian-elimination steps over F2 using only even-parity rows
    (rows k..m-1) as pivots.  Pivot rows are chosen in descending order of odd-parity-row elimination potential
    (i.e., densest even rows first) to zero out as many 1s as possible.
    """
    A2 = A.copy().astype(np.uint8) % 2
    m, n = A2.shape

    # Identify even-parity rows (rows k..m-1) and compute their weights
    parities = A2.sum(axis=1) % 2
    even_rows = [i for i in range(k, m) if parities[i] == 0]
    # Rank them by descending weight (sum of 1s)
    even_rows.sort(key=lambda i: A2[i].sum(), reverse=True)

    pivot_cols = set()
    pivots_done = 0
    for i in even_rows:
        if pivots_done >= max_pivots:
            break
        # Find a fresh pivot column in this row
        for j in range(n):
            if j not in pivot_cols and A2[i, j] == 1:
                # Use column j as a pivot
                pivot_cols.add(j)
                # Eliminate that column from all other rows
                for r in range(m):
                    if r != i and A2[r, j] == 1:
                        A2[r, :] ^= A2[i, :]
                pivots_done += 1
                break
    return A2


def solve_block_structure(A, k, max_p=None, time_limit=None):
    """
    (Unchanged) Solve for optimal p×p anti-diagonal block permutation on A with parity restrictions:
    - First k rows must be odd-parity (and only go to U or block 1).
    - Even rows (k..m-1) cannot go into U.
    """
    m, n = A.shape
    parities = A.sum(axis=1) % 2
    if not np.all(parities[:k] == 1):
        raise ValueError("First k rows must be odd-parity")
    if not np.all(parities[k:] == 0):
        raise ValueError("Rows k..m-1 must be even-parity")

    if max_p is None:
        max_p = min(m, n)

    best = {'p': None, 'obj': float('inf'), 'row_group': None, 'col_group': None}

    for p in range(1, max_p + 1):
        prob = pulp.LpProblem(f"block_struct_p{p}", pulp.LpMinimize)
        x = pulp.LpVariable.dicts('x', (range(m), range(p+1)), cat='Binary')
        y = pulp.LpVariable.dicts('y', (range(n), range(1, p+1)), cat='Binary')
        M = pulp.LpVariable('M', lowBound=0, cat='Integer')

        prob += pulp.lpSum(x[a][0] for a in range(m)) + M

        # Row assignments
        for a in range(m):
            prob += pulp.lpSum(x[a][r] for r in range(p+1)) == 1
        # Odd rows → U or block 1 only
        for a in range(min(k, m)):
            for r in range(2, p+1):
                prob += x[a][r] == 0
        # Even rows not to U
        for a in range(k, m):
            prob += x[a][0] == 0

        # Column assignments
        for b in range(n):
            prob += pulp.lpSum(y[b][c] for c in range(1, p+1)) == 1

        # Bound M
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

        solver = pulp.PULP_CBC_CMD(timeLimit=time_limit) if time_limit else None
        prob.solve(solver)

        if prob.status == pulp.LpStatusOptimal:
            obj = pulp.value(prob.objective)
            if obj < best['obj']:
                best['p'], best['obj'] = p, obj
                row_group = np.zeros(m, int)
                for a in range(m):
                    for r in range(p+1):
                        if pulp.value(x[a][r]) > 0.5:
                            row_group[a] = r
                col_group = np.zeros(n, int)
                for b in range(n):
                    for c in range(1, p+1):
                        if pulp.value(y[b][c]) > 0.5:
                            col_group[b] = c
                best['row_group'], best['col_group'] = row_group, col_group

    return best


def solve_with_elimination(A, k, max_p=None, time_limit=None, elimination_levels=None):
    """
    Try multiple partial-elimination levels before solving the block-structure MIP,
    and pick the outcome with smallest objective.
    ``elimination_levels`` is a list of integers (number of pivot steps).
    """
    m = A.shape[0]
    # Default levels: no elimination, half, full
    if elimination_levels is None:
        parities = A.sum(axis=1) % 2
        total_even = int(np.sum(parities[k:]))
        elimination_levels = [0, total_even // 2, total_even]

    best_overall = {'obj': float('inf')}
    best_data = {}

    for lvl in elimination_levels:
        # Apply partial elimination
        A2 = A.copy() if lvl == 0 else gaussian_eliminate_partial(A, k, lvl)
        # Solve block structure
        res = solve_block_structure(A2, k, max_p=max_p, time_limit=time_limit)
        if res['obj'] < best_overall['obj']:
            best_overall = res.copy()
            best_data = {
                'elimination_level': lvl,
                'A2': A2
            }

    best_overall.update(best_data)
    return best_overall


def get_permuted_matrix(A, row_group, col_group):
    """(Same as before)"""
    m, n = A.shape
    p = max(col_group)
    ordered_rows = np.hstack([np.where(row_group == r)[0] for r in range(p+1)])
    ordered_cols = np.hstack([np.where(col_group == c)[0] for c in range(1, p+1)])
    P = A[ordered_rows][:, ordered_cols]
    row_counts = [np.sum(row_group == r) for r in range(p+1)]
    col_counts = [np.sum(col_group == c) for c in range(1, p+1)]
    return P, np.cumsum(row_counts), np.cumsum(col_counts)


if __name__ == '__main__':
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

    A, k = tri49, 0

    # Run heuristic + optimization
    result = solve_with_elimination(A, k, max_p=3, time_limit=1000)
    print(f"Best elimination steps: {result['elimination_level']}")

    P, row_splits, col_splits = get_permuted_matrix(
        result['A2'], result['row_group'], result['col_group']
    )
    # (Then write to file as before)

    output_path = 'compress/permuted_visualization_per.txt'
    with open(output_path, 'w') as f:
        f.write(f"Best p: {result['p']}\n")
        f.write(f"Objective: {result['obj']}\n\n")
        f.write("Row groups (0=U,1..p): " + ' '.join(map(str, result['row_group'])) + "\n")
        f.write("Column partitions (1..p): " + ' '.join(map(str, result['col_group'])) + "\n\n")
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
