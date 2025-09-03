import numpy as np
import ldpc
from ldpc import bposd_decoder
from bposd.css import css_code

def GBB_code(ell, m, p1, p2):
    # returns the parity check matrices of a bivariate bicycle code arXiv:2308.07915
    # unrestricted weight and no mixed terms
    # define cyclic shift matrices
    I_ell = np.identity(ell, dtype=int)
    I_m = np.identity(m, dtype=int)
    x = np.kron(np.roll(I_ell, 1, axis=1), I_m)
    y = np.kron(I_ell, np.roll(I_m, 1, axis=1))

    A = np.zeros([ell*m, ell*m])
    B = np.zeros_like(A)
    # define parity check matrices
    for monomial in p1:
        A += np.linalg.matrix_power(x, monomial[0]) @ np.linalg.matrix_power(y, monomial[1])
    for monomial in p2:
        B += np.linalg.matrix_power(x, monomial[0]) @ np.linalg.matrix_power(y, monomial[1])
    A = A % 2
    B = B % 2
    AT = np.transpose(A)
    BT = np.transpose(B)

    hx = np.hstack((A, B))
    hz = np.hstack((BT, AT))

    return hx, hz


def CBB_code(ell, m, p1, p2):
    # l,m are two primes
    # p1 = (a,b,c), p2 = (d,e,f)
    # A = \sum xy^p1[i]
    # B = \sum xy^p2[i]
    # define cyclic shift matrices
    I_ell = np.identity(ell, dtype=int)
    I_m = np.identity(m, dtype=int)
    x = np.kron(np.roll(I_ell, 1, axis=1), I_m)
    y = np.kron(I_ell, np.roll(I_m, 1, axis=1))

    A = np.zeros([ell*m, ell*m])
    B = np.zeros_like(A)
    # define parity check matrices
    for i in p1:
        A += np.linalg.matrix_power(x@y, i)
    for i in p2:
        B += np.linalg.matrix_power(x@y, i)
    A = A % 2
    B = B % 2
    AT = np.transpose(A)
    BT = np.transpose(B)

    hx = np.hstack((A, B))
    hz = np.hstack((BT, AT))

    return hx, hz

class qLDPC:
    def __init__(self, type, l, m, p1, p2, d):
        self.type = type
        self.l = l
        self.m = m
        self.p1 = p1
        self.p2 = p2
        if type == 'BB' or type == 'bb':
            self.hx, self.hz = GBB_code(l, m, p1, p2)
        elif type == 'CBB' or type == 'cbb':
            self.hx, self.hz = CBB_code(l, m, p1, p2)
        else:
            raise ValueError('Unsupported type')
        self.supp = self.get_supports()
        code = css_code(self.hx, self.hz)
        self.lx = code.lx
        self.lz = code.lz
        self.k = code.K
        self.n = code.N
        self.d = d

    def get_supports(self):
        support = dict()
        for i in range(self.l*self.m):
            support['X' + str(i)] = []
            support['Z' + str(i)] = []
            for j in range(2*self.l*self.m):
                if j < self.l*self.m:
                    name = 'L' + str(j)
                else:
                    name = 'R' + str(j - self.l*self.m)
                if self.hx[i][j] == 1:
                    support['X' + str(i)].append(name)
                if self.hz[i][j] == 1:
                    support['Z' + str(i)].append(name)
        return support