import numpy as np
import random as rand
import time
import os

def ex(base, e, mod):
    res = 1
    pointer = 1
    while (pointer < e):
        pointer <<= 1
    while (pointer):
        res = (res * res) % mod
        if (e & pointer):
            res = (res * base) % mod
        pointer >>= 1
    return res

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

def MillerRabin(n):
    if (n == 1):
        return 2
    s = 1
    pt = 2
    while (n & pt == 0):
        pt <<= 1
        s += 1
    d = n >> s
    k = s
    while (n > pt):
        pt <<= 1
        k += 1
    for t in range(k):
        a = rand.randint(0, n - 1)
        if (ex(a, n - 1, n) != 1):
            return 0
        e = ex(a, d, n)
        if (e == 1 or e == n - 1):
            continue
        c = 1
        for i in range(s - 1):
            e = (e * e) % n
            if (e == n - 1):
                c = 0
                break
        if (c):
            return 0
    return 1

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
                if (j == m - 1):
                    print(matrix)
                    print(inverse)
                    os.system("pause")
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

def v2n(vector, m, p):
    number = 0
    for i in range(m):
        number *= p
        number += vector[i]
    return number

def n2v(number, m, p):
    vector = np.zeros(m, dtype="object")
    for i in range(m):
        vector[m - 1 - i] = number % p
        number //= p
    return vector

################################################################
# setup

bit = 3
m = 4

print("Enter the bits of the prime p.", bit)
print("Enter the degree of the polynomial.", m)

start = time.time()
# 소수 판정
print("Finding a " + str(bit) + "-bit prime number...")
if (bit == 1):
    p = 2
else:
    p = rand.randint(2 ** (bit - 2), 2 ** (bit - 1) - 1) * 2 + 1
    while (MillerRabin(p) == 0):
        p += 2
        if (p > 2 ** bit):
            p = 2 ** (bit - 1) + 1
print("p =", p)
ptime = time.time()
print()

# 최대주기 소인수분해
gp = p ** m - 1
print("Factoring " + str(gp) + "...")
composite = []
factor = []
if (p != 2):
    factor.append(2)
    print(2)
    g = p - 1
    while (g & 1 == 0):
        g >>= 1
    if (g != 1):
        if (MillerRabin(g)):
            factor.append(g)
            print(g)
        else:
            composite.append(g)
pt = 1
while (pt < m):
    g = p ** pt + 1
    while (g & 1 == 0):
        g >>= 1
    if (g != 1):
        if (MillerRabin(g)):
            factor.append(g)
            print(g)
        else:
            composite.append(g)
    pt <<= 1
    assert (pt <= m)
for g in composite:
    i = 3
    while (g >= i * i):
        if (g % i == 0):
            factor.append(i)
            print(i)
            g //= i
            while (g % i == 0):
                g //= i
            if (MillerRabin(g)):
                break
        i += 2
    if (g != 1):
        factor.append(g)
        print(g)
gp = p ** m - 1
ftime = time.time()
print()

# 점화식 생성 (p가 합성수일 경우 무조건 무한 루프에 빠짐)
print("Finding a primitive polynomial...")
basic = np.zeros((m, m), dtype="object")
one = np.eye(m, dtype="object")
for i in range(m - 1):
    basic[i][i + 1] = 1
    basic[m - 1][i] = rand.randint(0, 1)
trial = 0
while (1):
    trial += 1
    if (trial % bit == 0):
        print("trial #" + str(trial))
    r = rand.randint(0, m - 1)
    basic[m - 1][r] = (basic[m - 1][r] + 1) % p
    if (r == m - 1 and basic[m - 1][m - 1] == 0):
        basic[m - 1][m - 1] == 1
    if (np.array_equal(matexp(basic, m, gp, p), one) == False):
        continue
    primitive = 1
    for f in factor:
        if (np.array_equal(matexp(basic, m, gp // f, p), one) == True):
            primitive = 0
            break
    if (primitive == 1):
        break
print(basic[m - 1])
otime = time.time()
print()

#공개키 완성
print("Completing keys...")
while (1):
    s = rand.randint(1, gp - 1)
    rp = 1
    for f in factor:
        if (s % f == 0):
            rp = 0
            break
    if (rp):
        break
while (1):
    x = rand.randint(1, gp - 1)
    rp = 1
    for f in factor:
        if (x % f == 0):
            rp = 0
            break
    if (rp):
        break
shift = matexp(basic, m, s, p) % p
a0 = np.ones(m, dtype="object")
while (1):
    rp = 0
    for i in range(m):
        a0[i] = rand.randint(0, p - 1)
        if (a0[i] != 0):
            rp = 1
    if (rp):
        break
ax = (a0 @ matexp(np.transpose(shift), m, x, p)) % p
print("p = " + str(p))
print("m = " + str(m))
print("S = ")
print(shift)
print("A0 = " + str(a0))
print("Ax = " + str(ax))
print("and the secret key is " + str(x))
end = time.time()

################################################################
# encryption

print()
print()

print("Enter the prime p.", p)
print("Enter the degree of the polynomial.", m)

#gp = p ** m - 1

#shift = np.zeros((m, m), dtype="object")

print("Enter the shift matrix.")
print(shift)

#a0 = np.zeros(m, dtype="object")
#ax = np.zeros(m, dtype="object")

print("Enter A0.")
print(a0)
print("Enter Ax.")
print(ax)

M = rand.randint(0, gp - 1)
print("Enter the message.", M)
#M = int(input("Enter the message. "))

encstart = time.time()
#Ak
k = rand.randint(1, gp - 1)
ak = (a0 @ matexp(np.transpose(shift), m, k, p)) % p

#Akx (모방행렬 이용)
n0 = np.zeros((m, m), dtype="object")
nx = np.zeros((m, m), dtype="object")
for i in range(m):
    n0[0][i] = a0[i]
    nx[0][i] = ax[i]
for j in range(1, m):
    n0[j] = (n0[j - 1] @ np.transpose(shift)) % p
    nx[j] = (nx[j - 1] @ np.transpose(shift)) % p
ix = (nx @ inv(n0, p)) % p
akx = (a0 @ matexp(np.transpose(ix), m, k, p)) % p

c1 = v2n(ak, m, p)
c2 = (v2n(akx, m, p) + M) % gp

print("c1 =", c1)
print("c2 =", c2)
encend = time.time()

################################################################
# decryption

print()

print("Enter the prime p.", p)
print("Enter the degree of the recurrence.", m)

#gp = p ** m - 1

#shift = np.zeros((m, m), dtype="object")

print("Enter the shift matrix.")
print(shift)

#a0 = np.zeros(m, dtype="object")
#ax = np.zeros(m, dtype="object")

print("Enter A0.")
print(a0)
print("Enter Ax.")
print(ax)
print("Enter c1.")
print(c1)
print("Enter c2.")
print(c2)

print("Enter the private key.", x)

decstart = time.time()
ak = n2v(c1, m, p)

#Akx (모방행렬 이용)
n0 = np.zeros((m, m), dtype="object")
nk = np.zeros((m, m), dtype="object")
for i in range(m):
    n0[0][i] = a0[i]
    nk[0][i] = ak[i]
for j in range(1, m):
    n0[j] = (n0[j - 1] @ np.transpose(shift)) % p
    nk[j] = (nk[j - 1] @ np.transpose(shift)) % p
ik = (nk @ inv(n0, p)) % p
akx = (a0 @ matexp(np.transpose(ik), m, x, p)) % p

#해독
M = (c2 - v2n(akx, m, p)) % gp
print("The message is", M)

decend = time.time()

################################################################
print()
print()

print(int((ptime - start) * 1000), "ms for a prime")
print(int((ftime - ptime) * 1000), "ms for factorization")
print(int((otime - ftime) * 1000), "ms for a primitive polynomial")
print(int((end - otime) * 1000), "ms for a key")
print(int((end - start) * 1000), "ms for setup")
print(int((encend - encstart) * 1000), "ms for enc")
print(int((decend - decstart) * 1000), "ms for dec")
os.system("pause")
