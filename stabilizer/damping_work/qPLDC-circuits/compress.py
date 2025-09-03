from ortools.sat.python import cp_model

def solve_block_cp_sat(A, k, pmax, time_limit_s=20.0):
    m, n = A.shape
    best = None

    for p in range(1, pmax + 1):
      for r in range(k, m + 1):
        # 1) Build model
        model = cp_model.CpModel()

        # ... define variables and constraints here ...

        # 2) Create and configure solver
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = time_limit_s  # ← correct placement

        # 3) Solve
        result = solver.Solve(model)
        if result in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            obj = solver.ObjectiveValue()
            # ... process solution ...
            best = {'obj': obj, 'p': p, 'r': r}
    return best

if __name__ == '__main__':
    import numpy as np
    A = np.random.randint(0, 2, (8, 10))
    sol = solve_block_cp_sat(A, k=2, pmax=2, time_limit_s=20.0)
    print("Best solution:", sol)
