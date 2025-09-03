import numpy as np
import pulp

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

        # Variables:
        # x[a][r]: row a in group r (0 = U, 1..p = partitions)
        x = pulp.LpVariable.dicts('x', (range(m), range(p+1)), cat='Binary')
        # y[b][c]: column b in partition c (1..p)
        y = pulp.LpVariable.dicts('y', (range(n), range(1, p+1)), cat='Binary')
        # M: max size of any anti-diagonal row partition
        M = pulp.LpVariable('M', lowBound=0, cat='Integer')

        # Objective: |U| + M
        prob += pulp.lpSum(x[a][0] for a in range(m)) + M

        # Each row assigned to exactly one group
        for a in range(m):
            prob += pulp.lpSum(x[a][r] for r in range(p+1)) == 1

        # Restrict first k rows: they can only go to U (r=0) or first partition (r=1)
        for a in range(min(k, m)):
            for r in range(2, p+1):
                prob += x[a][r] == 0

        # Each column assigned to exactly one partition
        for b in range(n):
            prob += pulp.lpSum(y[b][c] for c in range(1, p+1)) == 1

        # M bounds: M >= number of rows in each anti-diagonal partition
        for r in range(1, p+1):
            prob += M >= pulp.lpSum(x[a][r] for a in range(m))

        # Anti-diagonal zero constraints
        for a in range(m):
            for b in range(n):
                if A[a, b] == 1:
                    for r in range(1, p+1):
                        for c in range(1, p+1):
                            if r + c != p + 1:
                                # row a in partition r and col b in partition c cannot both be 1
                                prob += x[a][r] + y[b][c] <= x[a][0] + 1

        # Solve with optional time limit
        if time_limit:
            prob.solve(pulp.PULP_CBC_CMD(timeLimit=time_limit))
        else:
            prob.solve()

        # Check for optimal and update best
        if prob.status == pulp.LpStatusOptimal:
            obj = pulp.value(prob.objective)
            if obj < best['obj']:
                best.update({ 'obj': obj, 'p': p })

                # Extract assignments
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


if __name__ == "__main__":
    # Example usage
    A49to1 = np.array([
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
    A = A49to1
    k = 1
    result = solve_block_structure(A, k=k, max_p=5, time_limit=60)
    print(f"Best p: {result['p']}, objective: {result['obj']}")
    print("Row assignments (0=U,1..p):", result['row_group'])
    print("Column partitions (1..p):", result['col_group'])
