import numpy as np
import pulp
import random

def solve_block_structure(A, k, max_p=None, time_limit=None):
    """
    Exact MILP solver for the block-structure problem:
      min |U| + M, subject to anti-diagonal zero-blocks
      and the first-k-rows restriction.
    Returns a dict with keys 'p','obj','row_group','col_group'.
    """
    m, n = A.shape
    if max_p is None:
        max_p = min(m, n)

    best = {'p': None, 'obj': float('inf'),
            'row_group': None, 'col_group': None}

    for p in range(1, max_p + 1):
        prob = pulp.LpProblem(f"block_struct_p{p}", pulp.LpMinimize)

        # Decision variables
        x = pulp.LpVariable.dicts('x',
            (range(m), range(p+1)), cat='Binary')    # row a → group r
        y = pulp.LpVariable.dicts('y',
            (range(n), range(1, p+1)), cat='Binary')  # col b → part c
        M = pulp.LpVariable('M', lowBound=0, cat='Integer')

        # Objective: |U| + M
        prob += pulp.lpSum(x[a][0] for a in range(m)) + M

        # Exactly-one assignment for rows & columns
        for a in range(m):
            prob += pulp.lpSum(x[a][r] for r in range(p+1)) == 1
        for b in range(n):
            prob += pulp.lpSum(y[b][c] for c in range(1, p+1)) == 1

        # First-k rows restriction: rows 0..k-1 only in U (r=0) or first block (r=1)
        for a in range(min(k, m)):
            for r in range(2, p+1):
                prob += x[a][r] == 0

        # M bounds the size of each anti-diagonal row block
        for r in range(1, p+1):
            prob += M >= pulp.lpSum(x[a][r] for a in range(m))

        # Anti-diagonal zero constraints
        for a in range(m):
            for b in range(n):
                if A[a, b] == 1:
                    for r in range(1, p+1):
                        for c in range(1, p+1):
                            if r + c != p + 1:
                                prob += x[a][r] + y[b][c] <= x[a][0] + 1

        # Solve
        if time_limit is not None:
            prob.solve(pulp.PULP_CBC_CMD(timeLimit=time_limit))
        else:
            prob.solve()

        if prob.status == pulp.LpStatusOptimal:
            obj = pulp.value(prob.objective)
            if obj < best['obj']:
                best['p'] = p
                best['obj'] = obj
                # extract assignments
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
    Permute A by row_group & col_group, and compute split indices.
    Returns (P, row_splits, col_splits).
    """
    m, n = A.shape
    p = int(col_group.max())

    # new row order: group 0 first, then 1..p
    ordered_rows = np.hstack([np.where(row_group == r)[0]
                              for r in range(p+1)])
    # new col order: 1..p
    ordered_cols = np.hstack([np.where(col_group == c)[0]
                              for c in range(1, p+1)])
    P = A[ordered_rows][:, ordered_cols]

    row_counts = [np.sum(row_group == r) for r in range(p+1)]
    col_counts = [np.sum(col_group == c) for c in range(1, p+1)]
    row_splits = np.cumsum(row_counts)  # boundaries after each row-part
    col_splits = np.cumsum(col_counts)

    return P, row_splits, col_splits


def generate_random_basis(m, num_ops, k):
    """
    Build an invertible T ∈ GL(m,2) by random swaps/XORs,
    but never move a row across the first-k boundary.
    Returns (T, ops) where ops is list of ('swap',i,j) or ('xor',i,j).
    """
    T = np.eye(m, dtype=int)
    ops = []
    for _ in range(num_ops):
        # pick a valid pair (i,j)
        while True:
            i, j = random.sample(range(m), 2)
            # allow swap only if both in top-k or both in bottom
            if random.random() < 0.5:
                if (i < k and j < k) or (i >= k and j >= k):
                    break
            else:
                # allow XOR only if source j is in same region as target i
                if i < k:
                    if j < k:
                        break
                else:
                    break

        if random.random() < 0.5:
            T[[i, j]] = T[[j, i]]
            ops.append(('swap', i, j))
        else:
            T[i] ^= T[j]
            ops.append(('xor', i, j))

    return T, ops


def solve_extended(A, k, max_p=None, time_limit=None,
                   iterations=30, ops_per_iter=10):
    """
    Heuristic search over random bases respecting the first-k constraint.
    Returns (best_obj, best_T, best_ops, best_B, best_res).
    """
    m, n = A.shape
    best = None

    for _ in range(iterations):
        T, ops = generate_random_basis(m, ops_per_iter, k)
        B = (T @ A) % 2

        res = solve_block_structure(B, k, max_p, time_limit)
        obj = res['obj']

        if best is None or obj < best[0]:
            best = (obj, T, ops, B, res)

    return best


if __name__ == "__main__":
    # -------------------------
    # Example usage
    # -------------------------
    # Define your input matrix A and k here:
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

    a14to2 = np.array([
        [0, 0, 0, 0, 1, 1, 1, 1, 1 ,1 ,1 ,0, 0, 0],
        [0, 0, 0, 1, 1, 1 ,0 ,0 ,1, 0, 1, 0, 1, 1],
        [1, 0, 0, 1, 1, 0 ,1 ,0 ,0, 1, 1, 1, 0, 1],
        [0, 1, 0, 1, 0, 1 ,1 ,1 ,0, 0, 1, 1, 1, 0],
        [0, 0, 1, 0, 0, 0 ,0 ,1 ,1, 1, 1, 1, 1, 1]
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
    k = 6
    
    max_p = 3
    time_limit = 40
    iterations = 100
    ops_per_iter = 15

    # Run the extended solver
    best_obj, best_T, best_ops, best_B, best_res = solve_extended(
        A, k, max_p, time_limit,
        iterations, ops_per_iter
    )

    # Build the final permuted matrix + splits
    P, row_splits, col_splits = get_permuted_matrix(
        best_B, best_res['row_group'], best_res['col_group']
    )

    # Write everything to a text file
    out = "extended_permuted.txt"
    with open(out, 'w') as f:
        f.write(f"Best objective: {best_obj}\n\n")

        f.write("Row operations performed:\n")
        for op, i, j in best_ops:
            if op == 'swap':
                f.write(f"  swap rows {i} ↔ {j}\n")
            else:
                f.write(f"  row {i} ^= row {j}\n")
        f.write("\n")

        f.write("Final row grouping (0=U,1..p):\n")
        f.write(" ".join(map(str, best_res['row_group'])) + "\n\n")

        f.write("Final column grouping (1..p):\n")
        f.write(" ".join(map(str, best_res['col_group'])) + "\n\n")

        f.write("Permuted matrix P with splitters:\n")
        for i, row in enumerate(P):
            elems = []
            for j, v in enumerate(row):
                elems.append(str(int(v)))
                if j+1 in col_splits[:-1]:
                    elems.append("|")
            line = " ".join(elems)
            f.write(line + "\n")
            if i+1 in row_splits[:-1]:
                f.write("-" * len(line) + "\n")

    print(f"Done. Results written to '{out}'.")
