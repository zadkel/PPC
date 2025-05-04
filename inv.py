import numpy as np
import random as rand
import time
import os

def div(x, mod):
    a = mod
    b = x
    i = 0
    j = 1
    while (b > 0):
        q = a // b
        a, b = b, a - q * b
        i, j = j, i - q * j
    return i % mod

def inv(mat, mod):
    m = np.shape(mat)[0]
    matrix = np.zeros((m, m), dtype="object")
    for i in range(m):
        matrix[i] = mat[i]
    inverse = np.eye(m, dtype="object")
    for i in range(m):
        if (matrix[i][i] == 0):
            for j in range(i + 1, m):
                if (matrix[i][j] != 0):
                    matrix[i] += matrix[j]
                    inverse[i] += inverse[j]
                    break
                assert (j == m - 1)
        n = div(matrix[i][i], mod)
        matrix[i] = (matrix[i] * n) % mod
        inverse[i] = (inverse[i] * n) % mod
        for j in range(i):
            inverse[j] = (inverse[j] - inverse[i] * matrix[j][i]) % mod
            matrix[j] = (matrix[j] - matrix[i] * matrix[j][i]) % mod
        for j in range(i + 1, m):
            inverse[j] = (inverse[j] - inverse[i] * matrix[j][i]) % mod
            matrix[j] = (matrix[j] - matrix[i] * matrix[j][i]) % mod
    return inverse

p = 11
m = 3
n0 = np.empty((m, m), dtype="object")
for i in range(m):
    for j in range(m):
        n0[i][j] = rand.randint(0, p - 1)
print("A")
print(n0)
i0 = inv(n0, p)
print("A^-1")
print(i0)
print("I")
print((n0 @ i0) % p)
