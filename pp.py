import numpy as np
import math
import random
import time

from math import floor, sqrt
try:
    long
except NameError:
    long = int

def mod(a, p):
    _, r = divmod(a, p)
    return r.astype(int)

def extEuclid(a, b):
    s0, s1, t0, t1 = 1, 0, 0, 1
    while b > 0:
        q, r = divmod(a, b)
        a, b = b, r
        s0, s1, t0, t1 = s1, s0 - q * s1, t1, t0 - q * t1
        pass
    return s0, t0, a

def mulinv(A, q):
    if type(A) != int:
        assert(A.ndim <= 2), 'The dimension of A must be less than 2'
    
    if type(A) == int:
        return extEuclid(A, q)[0] % q

    elif A.ndim == 1:
        res = np.zeros(A.shape)
        for i in range(len(A)):
            res[i] = extEuclid(A[i], q)[0] % q
            
   
            
    else:
        res = np.zeros(A.shape)
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                res[i][j] = extEuclid(A[i][j], q)[0] % q
                
    return res.astype(int)
 
def matpow(mat, k, p):
    result = np.eye(mat.shape[0])
    while k>0:
        
        if k&1 == 1:
            result = mod(np.matmul(result, mat), p)
        else:
            pass
        
        mat = mod(np.matmul(mat, mat), p)
        k = k >> 1
        
    return result.astype(int)

class PrimPoly(object):
    def __init__(self, p):
        assert(self.is_prime(p) == 1)
        self.p = p
    
    def is_prime(self, p):
        test = [2]
        result = 1
        
        for i in test:
            if pow(i, p-1, p) != 1:
                result = 0
                break
            
        if p == 2:
            result = 1
            
        if result == 1:
            pass
            # print('{} is prime'.format(p))
        else:
            print('{} is composite'.format(p))
            raise ValueError("p must be prime.")
    
        return result
    
    def genPoly(self, deg):
        coefarr = [1]
        
        for i in range(deg):
            coef = random.randint(0, self.p-1)
            coefarr.append(coef)
            
        return coefarr
        
    def genIV(self, poly):
        deg = len(poly) - 1
        iv = np.random.randint(0, self.p, size = (deg, 1))
        
        return iv
    
    '''
    def gcd(self, poly1, poly2):
        while len(poly1) >= len(poly2) or len(poly1) != 0:
            _, r = np.polydiv(poly1, poly2)
            poly1, poly2 = poly2, np.polydiv(poly1, poly2)
            
        return poly1
    '''
    
    def chrmat(self, poly):
        deg = len(poly) - 1
        poly = poly[1:]
        mat = np.zeros((deg, deg), dtype = np.int32)
        
        for i in range(deg-1):
            mat[i][i+1] = 1
        
        for j, x in enumerate(reversed(poly)):
            mat[deg-1][j] = -x
            
        return mat.astype(int) 
        
    def MMA(self,poly, A0, Ax, k):
        deg = len(poly) - 1
        mat = self.chrmat(poly)

        assert(A0.shape == (deg, 1) == Ax.shape),"The input shape should be {}.".format((deg, 1))
        assert(mat.shape == (deg, deg)),"The matrix shape should be {}.".format((deg, deg))
        
        saveA0 = A0[:]
        
        V = A0[:]
        b = Ax[:]
        P = b[:]
        
        for i in range(deg - 1):
            b = mod(np.matmul(mat, b), self.p)[:]
            A0 = mod(np.matmul(mat, A0), self.p)[:]
            
            P = np.hstack((P, b))
            V = np.hstack((V, A0))
        
        '''
        print("IV is : ", IV, end ='\n\n')
        print("P is : ", P, end = '\n\n')
        print("Q is : ", Q, end = '\n\n')
        '''
        
        Vinv = np.linalg.inv(V)
        M = mod(np.matmul(P, Vinv), p)
        
        Akx = matpow(M, k, self.p)
        Akx = mod(np.matmul(Akx, saveA0), p)
        
        return Akx
        
    def next_seq(self, poly, A_prev):
        tic = time.time()
        mat = self.chrmat(poly)
        toc = time.time()
        

        #print("Time consumed : ", toc - tic)
        
        return mod(np.matmul(mat, A_prev), self.p)
    
    def given_seq(self, poly, A0, x):
        mat = self.chrmat(poly)
        res = matpow(mat, x, self.p)
        
        return mod(np.matmul(res, A0), self.p)
    
    def complete_seq(self, poly, A0):
        
        print("Compute the whole sequence.")
        
        tic = time.time()
        
        res = A0[:]
        
        res = self.next_seq(poly, res)
        #print(res.T)
        
        while np.array_equal(res, A0) is not True:
            res = self.next_seq(poly, res)
            #print(res.T)
            
        toc = time.time()
        
        
        print("Time consumed : ", toc - tic)
    

class PPC(object):
    def __init__(self, poly, p, x, A0):
        self.poly = poly
        self.p = p
        self.x = x
        self.A0 = A0     
        
        self.deg = len(poly) - 1
        self.fsize = self.p ** self.deg
        
        pp = PrimPoly(self.p)
        
        self.mat = pp.chrmat(self.poly)
        self.Ax = np.matmul(matpow(self.mat, self.x, self.p), A0)
        
        assert(self.deg == A0.size)
        

    def padic(self, x):
        res = []

        while x>0:
            tmp = x%self.p
            res.append(tmp)
            x = (x - tmp) // self.p
        
        return res

    def toint(self, x):
        res = 0
        p = 1
        for i in reversed(x):
            res = res + int(i) * p
            p *= self.p

        return res
        
    def enc(self, m):
        k = random.randint(1, p-1)
        Ms = np.zeros((self.A0.shape))
        for i, x in enumerate(reversed(self.padic(m))):
            Ms[0, Ms.shape[1]-i-1] = int(x)

        c1 = np.matmul(matpow(self.mat, k, self.p), self.A0)
        Akx = pp.MMA(self.poly, self.A0, self.Ax, k)
        
        
        c2 = mod(Akx + Ms, self.p)
        
        return (c1, c2)

    def dec(self, x, ct):
        c1, c2 = ct
        Akx = pp.MMA(self.poly, self.A0, c1, x)
       
        pt = mod(c2 - Akx, self.p)

        ptstr = []
        for i in range(self.A0.shape[1]):
            ptstr.append(pt[i])

        pt = self.toint(ptstr)
            
        
        return pt
        #return perm, perm_sq
    
def setup(bit, deg, rand = True):
    
    tic = time.time()
    # Find p
    
    if rand:
        tic1 = time.time()
        print("Generating a " + str(bit + 1) + "-bit prime number")
        p = random.randint(2 ** (bit - 1), 2 ** bit - 1) * 2 + 1
        while pow(2, p-1, p) != 1 or pow(2, int((p-3)/2), int((p-1)/2)) != 1:
            p += 2
            if p > 2 ** (bit + 1):
                p == 2 ** bit + 1
        
        toc1 = time.time()        
        
        print("p = ",p)
        print("Time consumed for finding prime: ",toc1-tic1, end='\n\n')
        
    else:
        p = int(input("Input {}-bit prime number : ".format(bit+1)))
        
        while len(format(p, 'b')) <bit + 1:
            print("Your prime is too small. Choose again.")
            p = int(input())
            
        while pow(2, p-1, p) != 1 and p != 2:
            print("Input value is not prime. Choose again." )
            p = int(input())
    
    
    # Find gp
    tic2 = time.time()
    gp = p ** deg - 1
    ogp = gp

    
    print("Factoring " + str(gp))
    
    factor = ['2']
    while gp & 1 == 0:
        gp = gp >> 1
        
    i = 3
    while gp > i*i:
        if gp % i == 0:
            factor.append(str(i))
            gp = gp // i
            
            while gp % i == 0:
                gp = gp // i
                
            if pow(2, gp-1, gp) == 1:
                break
        i += 2
    
    if gp!= 1:
        factor.append(str(gp))
    
    toc2 = time.time()    
    print('Prime factor of gp :  ' + ' , '.join(factor))
    print("Time consumed for factorization : ", toc2-tic2, end = "\n\n")
    pp = PrimPoly(p)

    # Find primitive polynomial
    tic3 = time.time()
    
    one = np.eye(deg)
    
    while True:
        poly = pp.genPoly(deg)
        mat = pp.chrmat(poly).astype(int)
        if np.array_equal(matpow(mat, ogp,p), one):
            continue
        
        primitive = 1
        for f in factor:
            if np.array_equal(matpow(mat, int(ogp/int(f)), p), one):
                primitive = 0
                break
            
        if primitive:
            break
            
        
    toc3 = time.time()
    print("Primitive polynomial : ", poly)
    print("Time consumed for finding polynomial : ", toc3-tic3, end = "\n\n")

     
    
    
    toc = time.time()
    
    print('Time consumed for entire process : ', toc - tic)
    
    return p, factor, poly

def initvec(poly):
    length = len(poly) - 1
    result = np.zeros((length, 1))
    result[length-1, 0] = 1
    
    return result

def str2byte(msg):
    res = ''
    for i in msg:
        blk = format(ord(i), 'x')
        res += blk
        
    return res

def byte2str(byte):
    res = ''
    assert(len(byte) % 2 == 0), 'Byte length must be even.'
    size = len(byte) // 2
    
    for i in range(size):
        blk = int('0x' + byte[2 * i : 2 * (i+1)], 16)
        blk = chr(blk)
        res += blk
        
    return res
        


if __name__ == '__main__':
    bit = 1
    deg = 256
    
    setparam = input("Need setup? [Y/N] : ")
    
    if setparam == 'Y' or NameError:
        print("Setting up parameters.")
        p, factor, poly = setup(bit, deg, rand = False)
    
    else:
        pass

    pp = PrimPoly(p)
    
    #msg = input("Input plaintext : ")
    #msg = format(str2byte(msg), 'i')
    
    msg = 1732
    
    x = int(input("Input secret key : "))
    
    A0 = initvec(poly)
    
    crypt = PPC(poly, p, x, A0)
    (c1, c2) = crypt.enc(msg)
    
    pt = crypt.dec(x, (c1,c2))
    
    print('Initial message : ', msg)
    print('Decrypted message : ', pt)
    
    
    


'''
def fac(n):
    tic = time.time()
    step = lambda x: 1 + (x<<2) - ((x>>1)<<1)
    maxq = long(floor(sqrt(n)))
    d = 1
    q = n % 2 == 0 and 2 or 3 
    while q <= maxq and n % q != 0:
        q = step(d)
        d += 1
    toc = time.time()
    
    print("Time Consumed : ", toc-tic)
    
    return q <= maxq and [q] + fac(n//q) or [n]
'''

#cpt = PPC(poly, p, x, A0)


'''
# 점화식 생성
print("Finding a primitive recurrence...")
basic = np.zeros((m, m), dtype="object")
one = np.eye(m, dtype="object")
for i in range(m - 1):
    basic[i][i + 1] = 1
while (1):
    basic[m - 1][rand.randint(0, m - 1)] += 1
    if (np.array_equal(matexp(basic, m, gp, p), one) == False):
        continue
    primitive = 1
    for f in factor:
        if (np.array_equal(matexp(basic, m, gp, p / f), one)):
            primitive = 0
            break
    if (primitive):
        break
print(basic[m - 1])
print()
'''
