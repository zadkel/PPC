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

def gcd(x, mod):
    a = mod
    b = x
    while (b > 0):
        a, b = b, a % b
    return a

bit = int(input("Enter the bits of the message. "))
m = int(input("Enter the degree of the recurrence. "))

# 간단한 소수 판정
print("Finding a " + str(bit) + "-bit prime number...")
if (bit == 1):
    p = 2
else:
    p = rand.randint(2 ** (bit - 2), 2 ** (bit - 1) - 1) * 2 + 1
    while (ex(2, p - 1, p) != 1):
        p += 2
        if (p > 2 ** (bit + 1)):
            p = 2 ** bit + 1
print("p =", p)
print()

# 최대주기 소인수분해
start = time.time()
gp = p ** m - 1
print("Factoring " + str(gp) + "...")
cofactor = p - 1
while (cofactor & 1 == 0):
    cofactor //= 2
composite = []
if (p != 2):
    composite.append(cofactor)
pointer = 1
while (pointer != m):
    cofactor = p ** pointer + 1
    while (cofactor & 1 == 0):
        cofactor //= 2
    composite.append(cofactor)
    pointer <<= 1
    assert (pointer <= m)
factor = []
if (p != 2):
    factor.append(2)
    print(2)
for g in composite:
    i = 3
    while (g > i * i):
        if (g % i == 0):
            factor.append(i)
            print(i)
            g //= i
            while (g % i == 0):
                g //= i
            if (ex(2, g - 1, g) == 1):
                break
        i += 2
    if (g != 1):
        factor.append(g)
        print(g)
factortime = time.time()
print("Factoring time :", (int)((factortime - start) * 1000), "ms")
print()


# 점화식 생성 (p가 합성수일 경우 무조건 무한 루프에 빠짐)
print("Finding a primitive recurrence...")
basic = np.zeros((m, m), dtype="object")
one = np.eye(m, dtype="object")
for i in range(m - 1):
    basic[i][i + 1] = 1
while (1):
    for i in range(m - 1):
        basic[m - 1][i] = rand.randint(0, p - 1)
    basic[m - 1][m - 1] = rand.randint(1, p - 1)
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
polytime = time.time()
print("Polynomial finding time :", (int)((polytime - factortime) * 1000), "ms")
print()

#공개키 완성
print("Completing keys...")
while (1):
    x = rand.randint(1, gp - 1)
    rp = 1
    for f in factor:
        if (x % f == 0):
            rp = 0
            break
    if (rp):
        break
a0 = np.zeros(m, dtype="object")
for i in range(m):
    a0[i] = rand.randint(1, p - 1)
ax = (a0 @ matexp(np.transpose(basic), m, x, p)) % p
print("A0 = " + str(a0))
print("Ax = " + str(ax))
print("and the secret key is " + str(x))
end = time.time()
print("Key completing time :", (int)((end - polytime) * 1000), "ms")
print("Total time :", (int)((end - start) * 1000), "ms")
os.system("pause")
