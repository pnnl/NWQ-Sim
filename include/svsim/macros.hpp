#pragma once

/***********************************************
 * Key Macros
 ***********************************************/
#define PUT(arr, i, val) (arr[(i)] = (val))
#define GET(arr, i) (arr[(i)])
#define BARR  \
    while (0) \
    {         \
    };

// For C2 and C4 gates
#define DIV2E(x, y) ((x) >> (y))
#define MOD2E(x, y) ((x) & (((IdxType)1 << (y)) - (IdxType)1))
#define EXP2E(x) ((IdxType)1 << (x))
#define SV16IDX(x) (((x >> 3) & 1) * EXP2E(qubit0) + ((x >> 2) & 1) * EXP2E(qubit1) + ((x >> 1) & 1) * EXP2E(qubit2) + ((x & 1) * EXP2E(qubit3)))
