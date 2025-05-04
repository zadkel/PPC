import numpy as np
import random as rand
import time
import os

def matexp(matrix, m, e, mod):
    res = np.eye(m, dtype="object")
    pointer = 1
    while (pointer < e):
        pointer <<= 1
    while (pointer):
        res = (res @ res) % mod
        if (e & pointer):
            res = (res @ matrix) % mod
        pointer >>= 1
    return res

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


p = int(input("Enter the prime p. "))
m = int(input("Enter the degree of the recurrence. "))

gp = p ** m - 1

basic = np.zeros((m, m), dtype="object")
for i in range(m - 1):
    basic[i][i + 1] = 1

print("Enter the primitive polynomial.")
for i in range(m):
    basic[m - 1][i] = int(input())

a0 = np.zeros(m, dtype="object")
ax = np.zeros(m, dtype="object")
ak = np.zeros(m, dtype="object")
akxms = np.zeros(m, dtype="object")

print("Enter A0.")
for i in range(m):
    a0[i] = int(input())
print("Enter Ax.")
for i in range(m):
    ax[i] = int(input())
print("Enter Ak.")
for i in range(m):
    ak[i] = int(input())
print("Enter Akx + mv.")
for i in range(m):
    akxms[i] = int(input())

x = int(input("Enter the private key. "))

start = time.time()
#Akx (모방행렬 이용)
n0 = np.zeros((m, m), dtype="object")
nk = np.zeros((m, m), dtype="object")
for i in range(m):
    n0[0][i] = a0[i]
    nk[0][i] = ak[i]
for j in range(1, m):
    n0[j] = (n0[j - 1] @ np.transpose(basic)) % p
    nk[j] = (nk[j - 1] @ np.transpose(basic)) % p
ik = (nk @ inv(n0, p)) % p
akx = (a0 @ matexp(np.transpose(ik), m, x, p)) % p

#메시지 벡터
ms = (akxms - akx) % p

M = 0
for i in range(m):
    M = M * p + ms[i]

print("The message is", M)
end = time.time()
print("Total time :", (int)((end - start) * 1000), "ms")
os.system("pause")
