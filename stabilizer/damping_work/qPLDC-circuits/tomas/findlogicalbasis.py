# -------- robust input normalization (lists or numpy arrays) --------
import numpy as _np



def _to01_vector(v):
    # Accept list/tuple/range or numpy array; return list of 0/1 ints
    if _np is not None and isinstance(v, _np.ndarray):
        return [int(x) & 1 for x in v.astype(int).tolist()]
    return [int(x) & 1 for x in list(v)]

def _to01_matrix(M):
    # Accept list-of-lists or numpy 2D array; return list of list of 0/1 ints
    if _np is not None and isinstance(M, _np.ndarray):
        if M.ndim != 2:
            raise ValueError("Matrix must be 2D")
        return [[int(x) & 1 for x in row] for row in M.astype(int).tolist()]
    # assume iterable of iterables
    Ml = list(M)
    return [[int(x) & 1 for x in row] for row in Ml]

def _normalize_L(L):
    # Accept 1D or 2D (list or numpy). Return list-of-rows (each row list of 0/1).
    if _np is not None and isinstance(L, _np.ndarray):
        if L.ndim == 1:
            return [_to01_vector(L)]
        elif L.ndim == 2:
            return _to01_matrix(L)
        else:
            raise ValueError("L must be 1D or 2D")
    # Non-numpy path
    L_list = list(L)
    if len(L_list) == 0:
        raise ValueError("L must be non-empty")
    # If first element is scalar-like, treat as single row
    if not hasattr(L_list[0], "__iter__") or isinstance(L_list[0], (str, bytes)):
        return [_to01_vector(L_list)]
    # Otherwise assume rows
    return _to01_matrix(L_list)

def _transpose(M):
    return [list(col) for col in zip(*M)] if M else []

def _lincomb_rows_mod2(rows, coeffs):
    if not rows:
        return []
    n = len(rows[0])
    out = [0]*n
    for i, c in enumerate(coeffs):
        if (c & 1) == 1:
            ri = rows[i]
            for j in range(n):
                out[j] ^= (ri[j] & 1)
    return out

# Solve A x = b over GF(2). Returns one solution (list of 0/1) or None if infeasible.
def gf2_solve(A, b):
    m = len(A)
    n = len(A[0]) if m else 0
    aug = [A[i][:] + [b[i] & 1] for i in range(m)]
    row = 0
    pivcol = [-1]*n
    for col in range(n):
        sel = None
        for r in range(row, m):
            if aug[r][col] & 1:
                sel = r; break
        if sel is None: 
            continue
        if sel != row:
            aug[row], aug[sel] = aug[sel], aug[row]
        pivcol[col] = row
        for r in range(m):
            if r != row and (aug[r][col] & 1):
                for c in range(col, n+1):
                    aug[r][c] ^= aug[row][c]
        row += 1
        if row == m:
            break
    for r in range(m):
        if all((aug[r][c] & 1) == 0 for c in range(n)) and (aug[r][n] & 1):
            return None
    x = [0]*n
    for col in range(n-1, -1, -1):
        pr = pivcol[col]
        if pr == -1:
            x[col] = 0
        else:
            s = 0
            for j in range(col+1, n):
                s ^= (aug[pr][j] & 1) & (x[j] & 1)
            x[col] = (aug[pr][n] ^ s) & 1
    return x

def _resolve_allowed(n, *, basis=None, allowed=None, disallowed=None, allowed_symbol='0'):
    if disallowed is not None:
        D = sorted({int(i) for i in (disallowed.tolist() if _np is not None and isinstance(disallowed, _np.ndarray) else disallowed)})
        A = [i for i in range(n) if i not in D]
        return A
    if allowed is not None:
        return sorted({int(i) for i in (allowed.tolist() if _np is not None and isinstance(allowed, _np.ndarray) else allowed)})
    if basis is not None:
        # basis can be list/array of '0'/'+' (or 0/1). Keep it simple: treat exactly 'allowed_symbol' as allowed.
        b = list(basis) if not (_np is not None and isinstance(basis, _np.ndarray)) else basis.tolist()
        return [i for i,val in enumerate(b) if str(val) == allowed_symbol]
    raise ValueError("Provide one of: basis, allowed, or disallowed.")

# ---------------------- main CSS function ----------------------
def css_restrict_support(L, H, *, basis=None, allowed=None, disallowed=None, allowed_symbol='0'):
    """
    Find stabilizer-equivalent reps of logicals L using stabilizers H (same type),
    with support only on an allowed set of qubits.

    Inputs:
      - L: (k x n) or (n,) logical(s) of the same type as H (Z with Hz, or X with Hx)
      - H: (r x n) stabilizers of that type
      - basis/allowed/disallowed: choose one. If using basis, allowed_symbol is '0' (Z) or '+' (X).

    Returns: list of dicts (one per logical)
    """
    # Normalize inputs safely (no truthiness on arrays)
    L_rows = _normalize_L(L)
    H_ll = _to01_matrix(H)

    if len(H_ll) == 0:
        raise ValueError("H must be non-empty.")
    n = len(H_ll[0])
    if any(len(row) != n for row in H_ll):
        raise ValueError("All rows of H must have the same length n.")
    if any(len(row) != n for row in L_rows):
        raise ValueError("Each row of L must have length n.")

    A_allowed = _resolve_allowed(n, basis=basis, allowed=allowed, disallowed=disallowed, allowed_symbol=allowed_symbol)
    D = [i for i in range(n) if i not in A_allowed]

    # Build (H[:,D])^T
    H_cols_D = [[row[j] & 1 for j in D] for row in H_ll]  # r x |D|
    A_sys = _transpose(H_cols_D)                           # |D| x r

    results = []
    for l in L_rows:
        b_sys = [l[j] & 1 for j in D]
        beta = gf2_solve(A_sys, b_sys)                     # (H[:,D])^T beta = l_D
        if beta is None:
            results.append({
                'feasible': False, 'beta': None, 'rep': None,
                'support_indices': None, 'used_stabilizer_indices': None,
                'reason': "No combination of same-type stabilizers clears all disallowed qubits."
            })
            continue
        delta = _lincomb_rows_mod2(H_ll, beta)
        rep = [(l[j] ^ delta[j]) & 1 for j in range(n)]
        if any(rep[j] for j in D):
            results.append({
                'feasible': False, 
                'beta': None, 'rep': None,
                'support_indices': None, 'used_stabilizer_indices': None,
                'reason': "Internal check failed; please report."
            })
            continue 
        results.append({
            'feasible': True,
            'beta': beta,
            'rep': rep,
            'support_indices': [i for i,v in enumerate(rep) if v],
            'used_stabilizer_indices': [i for i,v in enumerate(beta) if v],
            'reason': ""
        })
    return results

# Convenience wrappers
def restrict_Z(Lz, Hz, *, basis=None, allowed=None, disallowed=None):
    return css_restrict_support(Lz, Hz, basis=basis, allowed=allowed, disallowed=disallowed, allowed_symbol='0')

def restrict_X(Lx, Hx, *, basis=None, allowed=None, disallowed=None):
    return css_restrict_support(Lx, Hx, basis=basis, allowed=allowed, disallowed=disallowed, allowed_symbol='+')

def findnewL(L,H,measure_raw,type):
    basis_all_zero = ['0'] * H.shape[1]
    basis_all_plus = ['+'] * H.shape[1]
    if type == 'Z':
        for i in measure_raw:
            basis_all_zero[i] = '+'
        newLz=[]
        out = restrict_Z(L, H, basis=basis_all_zero)  # allowed_symbol='0' by default
        for i,res in enumerate(out):
            newLz.append(res['rep'])
        return newLz
    elif type == 'X':
        for i in measure_raw:
            basis_all_plus[i] = '0'
        newLx=[]
        out = restrict_X(L, H, basis=basis_all_plus)
        for i,res in enumerate(out):
            newLx.append(res['rep'])
        return newLx
