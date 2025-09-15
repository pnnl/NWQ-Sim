from ortools.sat.python import cp_model
import numpy as np

def solve_block_cp_sat(A, k, pmax, time_limit_s=30):
    m, n = A.shape
    best = None

    for p in range(1, pmax+1):
      for r in range(k, m+1):
        model = cp_model.CpModel()
        model.parameters.max_time_in_seconds = time_limit_s

        # 1) Permutation booleans
        X = {}
        for i in range(m):
          for u in range(m):
            # i<k → u<k ; i>=k→u>=k
            if (i<k and u<k) or (i>=k and u>=k):
              X[(i,u)] = model.NewBoolVar(f"X_{i}_{u}")
            else:
              X[(i,u)] = None

        Y = { (j,v): model.NewBoolVar(f"Y_{j}_{v}")
              for j in range(n) for v in range(n) }

        # each original row i maps to exactly one u
        for i in range(m):
          valid_us = [u for u in range(m) if X[(i,u)] is not None]
          model.Add(sum(X[(i,u)] for u in valid_us) == 1)
        # each new position u receives exactly one original row
        for u in range(m):
          valid_is = [i for i in range(m) if X[(i,u)] is not None]
          model.Add(sum(X[(i,u)] for i in valid_is) == 1)

        # columns
        for j in range(n):
          model.Add(sum(Y[(j,v)] for v in range(n)) == 1)
        for v in range(n):
          model.Add(sum(Y[(j,v)] for j in range(n)) == 1)

        # 2) Block assignments for rows>r and all columns
        H = {}
        for u in range(m):
          for b in range(p):
            h = model.NewBoolVar(f"H_{u}_{b}")
            H[(u,b)] = h
            if u < r:
              model.Add(h == 0)   # row in U, not in L
          if u >= r:
            model.Add(sum(H[(u,b)] for b in range(p)) == 1)

        W = { (v,b): model.NewBoolVar(f"W_{v}_{b}")
              for v in range(n) for b in range(p) }
        for v in range(n):
          model.Add(sum(W[(v,b)] for b in range(p)) == 1)

        # 3) Anti-diagonal AND‑vars E[u,v,b] = H[u,b] ∧ W[v,p-1-b]
        E = {}
        for u in range(m):
          for v in range(n):
            for b in range(p):
              e = model.NewBoolVar(f"E_{u}_{v}_{b}")
              E[(u,v,b)] = e
              # e → H[u,b] and e → W[v,p-1-b]
              model.AddImplication(e, H[(u,b)])
              model.AddImplication(e, W[(v, p-1-b)])
              # H[u,b] ∧ W[v,p-1-b] → e
              model.AddBoolAnd([H[(u,b)], W[(v,p-1-b)]]).OnlyEnforceIf(e)

        # 4) Build A' entries via a two‑step AND for X & Y
        P = {}
        for u in range(m):
          for v in range(n):
            for i in range(m):
              for j in range(n):
                if X[(i,u)] is None: 
                  continue
                pvar = model.NewBoolVar(f"P_{i}_{j}_{u}_{v}")
                P[(i,j,u,v)] = pvar
                model.AddImplication(pvar, X[(i,u)])
                model.AddImplication(pvar, Y[(j,v)])
                model.AddBoolAnd([X[(i,u)], Y[(j,v)]]).OnlyEnforceIf(pvar)

        # 5) Z_b = sum_{u>=r,v,i,j} A[i,j]*P[i,j,u,v]*E[u,v,b]
        #    We linearize the product P⋀E via a third AND Q = P ∧ E
        Q = {}
        Zmax = model.NewIntVar(0, m*n, "Zmax")
        for b in range(p):
          # accumulate all terms for block b
          terms = []
          for u in range(r, m):
            for v in range(n):
              for i in range(m):
                for j in range(n):
                  if (i,j,u,v) not in P: 
                    continue
                  q = model.NewBoolVar(f"Q_{i}_{j}_{u}_{v}_{b}")
                  Q[(i,j,u,v,b)] = q
                  # q → P[i,j,u,v]
                  model.AddImplication(q, P[(i,j,u,v)])
                  # q → E[u,v,b]
                  model.AddImplication(q, E[(u,v,b)])
                  # P ∧ E → q
                  model.AddBoolAnd([
                    P[(i,j,u,v)],
                    E[(u,v,b)]
                  ]).OnlyEnforceIf(q)
                  if A[i,j] == 1:
                    terms.append(q)
          # Z_b ≤ Zmax
          zb = model.NewIntVar(0, m*n, f"Z_{b}")
          model.Add(zb == sum(terms))
          model.Add(Zmax >= zb)

        # 6) Objective: minimize r + Zmax
        #    r is constant in this subproblem
        model.Minimize(r + Zmax)

        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = time_limit_s
        result = solver.Solve(model)
        if result == cp_model.OPTIMAL or result == cp_model.FEASIBLE:
          obj = solver.ObjectiveValue()
          if best is None or obj < best['obj']:
            # extract perms
            row_perm = [ next(i for i in range(m)
                              if X[(i,u)] and solver.Value(X[(i,u)]) )
                         for u in range(m) ]
            col_perm = [ next(j for j in range(n)
                              if solver.Value(Y[(j,v)]))
                         for v in range(n) ]
            best = {
              'obj': obj,
              'p': p, 'r': r,
              'row_perm': row_perm,
              'col_perm': col_perm
            }

    # Write out final best
    with open('solution_output.txt','w') as f:
      f.write(f"Best p: {best['p']}\n")
      f.write(f"Best |U| (r): {best['r']}\n")
      f.write(f"Row permutation: {best['row_perm']}\n")
      f.write(f"Col permutation: {best['col_perm']}\n")
      # you can also dump the permuted matrix here if you like

    return best




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

    A = a14to2

    res = solve_block_cp_sat(A, k=2, pmax=2, time_limit_s=20)
    print("Done; see solution_output.txt for details.")
    
