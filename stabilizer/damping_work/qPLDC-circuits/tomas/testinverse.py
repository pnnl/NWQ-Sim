import numpy as np

def convert_index(indices):
  new_indices = []
  for val in indices:
    if val >= 72:
      # Isolate the part of the value that corresponds to 6*x + y
      temp_val = val - 72
      # Extract x and y from the 6*x + y format
      x = temp_val // 6
      y = temp_val % 6
      # Calculate the new value in the 12*y + x format and add the offset back
      new_val = 72 + 12 * y + x
    else:
      # Direct conversion for values less than 72
      x = val // 6
      y = val % 6
      new_val = 12 * y + x
    new_indices.append(new_val)
  return new_indices

def computepolyfromTomas(val):
    x = val // 6
    y = val % 6
    return x, y

import numpy as np

def find_reversal_circuit(M):
    """
    Reduce an invertible binary matrix M (GF(2)) to a permutation matrix (one '1' per row/col)
    using only row additions (CNOTs) and no row swaps, to cut down gate count.

    Returns:
        gates: list[(control_row, target_row)]  # CNOT(control -> target)
        pivots: dict[int -> int]               # column j -> pivot row r
        P: np.ndarray                          # resulting permutation matrix
    Raises:
        ValueError if M is singular (no perfect matching).
    """
    m = M.copy() % 2
    n = m.shape[0]
    assert m.shape[0] == m.shape[1], "M must be square"

    # --- Phase 0: choose one pivot row per column via bipartite matching (rows <-> cols on 1-entries) ---
    # Build adjacency: for each column j, rows that have 1 in that column
    col_to_rows = [np.flatnonzero(m[:, j]).tolist() for j in range(n)]

    match_row_for_col = [-1] * n   # for column j, which row is assigned
    match_col_for_row = [-1] * n   # for row i, which column is assigned

    def try_assign(col, seen_rows):
        for r in col_to_rows[col]:
            if seen_rows[r]:
                continue
            seen_rows[r] = True
            # if row r is free or we can reassign its current column elsewhere
            if match_col_for_row[r] == -1 or try_assign(match_col_for_row[r], seen_rows):
                match_row_for_col[col] = r
                match_col_for_row[r] = col
                return True
        return False

    for j in range(n):
        if not try_assign(j, [False] * n):
            raise ValueError("Matrix is singular (no perfect row↔col matching).")

    pivots = {j: match_row_for_col[j] for j in range(n)}

    # We’ll define an order over pivot rows to mimic forward/backward elimination
    # without physically swapping rows: order index of row r is the column it pivots.
    row_order = {r: j for j, r in pivots.items()}  # row -> "position"

    gates = []

    # --- Phase 1: clear BELOW each pivot (forward elimination in pivot-order) ---
    # Process columns in increasing pivot-order (j = 0..n-1)
    for j in range(n):
        r = pivots[j]  # pivot row for column j
        # Make sure pivot bit is 1 (it is by construction)
        # Clear all rows with "larger order" (i.e., considered 'below') that have a 1 in column j
        for i in range(n):
            if i == r:
                continue
            # Only clear those considered "below" in our logical order
            if row_order[i] > row_order[r] and m[i, j] == 1:
                m[i, :] = (m[i, :] + m[r, :]) % 2
                gates.append((r, i))  # CNOT(r -> i)

    # --- Phase 2: clear ABOVE each pivot (backward elimination in reverse pivot-order) ---
    for j in range(n - 1, -1, -1):
        r = pivots[j]
        for i in range(n):
            if i == r:
                continue
            # Only clear those considered "above" in our logical order
            if row_order[i] < row_order[r] and m[i, j] == 1:
                m[i, :] = (m[i, :] + m[r, :]) % 2
                gates.append((r, i))  # CNOT(r -> i)

    # At this point, each column j has a single '1' at row r=pivots[j], and every row
    # appears as a pivot exactly once -> m is a permutation matrix.
    P = m
    return gates

def find_synthesis_circuit(matrix_to_generate):
    """
    Finds the CNOT sequence that generates the transformation matrix.
    This is the reverse of the reversal circuit.
    """
    # First, find the circuit that reverses the transformation
    reversal_gates = find_reversal_circuit(matrix_to_generate)
    
    # The synthesis circuit is simply this sequence in reverse order
    synthesis_gates = reversal_gates[::-1]
    
    return synthesis_gates