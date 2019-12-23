from math import log
from collections import Counter
from functools import reduce
import itertools

# 参数
class Constants:
    def __init__(self, B, k):
        self.B = B
        self.k = k
        self._calc()
    
    def _calc(self):
        self.alpha = 2**(1/3) - 1
        self.dB = int(self.alpha * self.B) # bound of d
        # eps = (k mod 9)/3
        eps = (self.k%9)//3
        self.eps = eps if eps==1 else -1
        # 素数表
        self.is_primes, self.smallest_primes, self.primes = eratosthenes_sieve(self.dB)
        # P, M
        self.P = int(3*log(log(self.B))*log(log(log(self.B))))
        self.Mps = []
        for prime in self.primes:
            if prime>self.P:
                break
            if prime>=5:
                self.Mps.append(prime)
        if not self.Mps:
            self.Mps = [5]
        self.M = reduce(lambda a, b: a*b, self.Mps)
        # 4^(-1) mod M
        self.rev_4_mod_M = modulo_inv(4, self.M)
        # z^3 mod M
        self.z_cube_mod_M = [(z**3)%self.M for z in range(self.M)]
        # table of (i|p) for p in mps
        self.mps_quad_residu_table = get_quadratic_residue_table(self.Mps)
        # table of z^3=k(mod p^r)
        self.get_cube_table()

    def get_cube_table(self):
        '''cube_table[(r, p)] is solutions of z^3=k(mod p^r)'''
        self.cube_table = {}
        for p in self.primes:
            zs = cuberoots(self.k, p)
            self.cube_table[(1, p)] = zs
            pr, r = p*p, 2
            while pr<=self.dB:
                zs = hensel(p, self.k, r, zs)
                self.cube_table[(r,p)] = zs
                pr, r = pr*p, r+1

class DSlice:
    def __init__(self, cons, d):
        self.d = d
        self.cons = cons

    def run(self):
        if self.d%3==0:
            return
        self._calc_d()
        self._calc_z_mod_18()
        self._calc_z_mod_d()
        if self.d*self.cons.M<=self.cons.dB:
            self._calc_z_mod_M()
            self._calc_z_mod_dM18()
        else:
            self._calc_z_mod_d18()
        return self._is_square()
        
    def _calc_d(self):
        # (d|3)
        d_residu_3 = self.d%3
        self.d_residu_3 = -1 if d_residu_3==2 else d_residu_3
        # sgnz
        self.sgnz = self.cons.eps*self.d_residu_3
        # d = 2^s*q
        self.s = get_log_2(self.d)
        # c = (sgnz*k+rev_4*(d**3))%M
        self.c = (self.sgnz*self.cons.k+self.cons.rev_4_mod_M*(self.d**3))%self.cons.M
        # |z|>d/alpha
        self.z_limit = int(self.d/self.cons.alpha)+1
    
    def _calc_z_mod_M(self):
        '''计算 z mod M的剩余类'''
        self.z_mod_M = []
        for z in range(self.cons.M):
            # z3c = |z|^3-c
            z3c = self.sgnz*self.cons.z_cube_mod_M[(z%self.cons.M)]-self.c
            for p in self.cons.Mps:
                # dsc1 = (3d|p); dsc2 = (|z|^3-c|p)
                dsc1 = self.cons.mps_quad_residu_table[((3*self.d)%p, p)]
                dsc2 = self.cons.mps_quad_residu_table[(z3c%p, p)]
                dsc = dsc1*dsc2
                if dsc==-1:
                    break
            else:
                self.z_mod_M.append(z)
    
    def _calc_z_mod_18(self):
        '''计算 z mod 18的剩余类'''
        self.z_mod_18 = [(4*(self.cons.k//3)*(2-self.d**2)+9*(self.cons.k+self.d))%18]

    def _calc_z_mod_d(self):
        '''计算 z mod d的剩余类, 主要来自于z^3=k(mod d)'''
        if self.d==1:
            self.z_mod_d = [0]
        else:
            d_pf_counter = prime_factorization(self.d, self.cons.primes, self.cons.smallest_primes)
            n = []
            a = []
            for prime, count in d_pf_counter.items():
                n.append(prime**count)
                a.append(self.cons.cube_table[count, prime])
            self.z_mod_d, _ = multivalued_chinese_remainder(n, a, self.d)

            # self.z_mod_d, _ = multivalued_chinese_remainder(keys, 
            #     [self.cons.cube_table[count, prime] for prime, count in d_pf_counter.items()], prod=self.d)

    def _calc_z_mod_d18(self):
        self.z_mod_dM18, self.divisor = extended_multivalued_chinese_remainder(self.d, 18, self.z_mod_d, self.z_mod_18)

    def _calc_z_mod_dM18(self):
        # 剩余定理求解z_mod_m，z_mod_d,z_mod_18,
        # 由于M和18 coprime，可以用中国剩余定理先计算，然后计算18M和d
        z_mod_18M, M18 = multivalued_chinese_remainder([18, self.cons.M], [self.z_mod_18, self.z_mod_M])
        self.z_mod_dM18, self.divisor = extended_multivalued_chinese_remainder(self.d, M18, self.z_mod_d, z_mod_18M)

    def _is_square(self):
        for z_repr in self.z_mod_dM18:
            z0 = self.sgnz*self.z_limit + self.sgnz*((self.sgnz*z_repr-self.z_limit)%self.divisor)
            # z0 = self.sgnz*self.z_limit +  self.sgnz*((z_repr-self.sgnz*self.z_limit)%self.divisor)
            for z in range(z0, self.sgnz*(self.cons.B+1), self.sgnz*self.divisor):
                # 条件2、3的判定和平方根判定
                if (not self._condition_23(z)) or (not self._condition_delta(z)):
                    continue
                # 到这一步就说明是解了
                x,y,z = self.getXYZ(z)
                if self.testXYZ(x,y,z):
                    print(x,y,z)
                    return True                    

    def _condition_23(self, z):
        t = get_log_2(z**3-self.cons.k)
        return t<=3*self.s and (self.s-t)%2==0
    
    def _condition_delta(self, z):
        delta = 3*self.d*(4*self.sgnz*(z**3-self.cons.k)-self.d**3)
        return is_square(delta)

    def getXYZ(self, z):
        '''获得方程的解
        由于3不能整除d,delta为平方数，由于delta能整除3,所以delta能整除9，
        意味着delta/3d能整除3，所以根号里面必然是整数解，因此不需要额外的判断'''
        sq = round(((4*abs(self.cons.k-z**3)-self.d**3)//(3*self.d))**0.5)
        x, y = -self.sgnz*(self.d+sq)//2, -self.sgnz*(self.d-sq)//2
        if abs(x)<abs(y):
            x, y = y, x
        return x, y, z   

    def testXYZ(self, x, y, z):
        return x**3+y**3+z**3==self.cons.k 
        
# =======素数和素数分解====================
def eratosthenes_sieve(n):
    '''获取n以内对素数和最小素数'''
    is_primes = [False] * (n+1) # 是否是素数（>=2以上成立）
    smallest_primes = [1] * (n+1)    # 该数的最小素数（>=2以上成立）
    # 主程序
    for i in range(2, n+1, 2):
        if i>2:
            is_primes[i] = True
        smallest_primes[i] = 2
    for i in range(3, n+1, 2):
        if not is_primes[i]:
            smallest_primes[i] = i
            for j in range(i**2, n+1, 2*i):
                if not is_primes[j]:
                    is_primes[j] = True
                    smallest_primes[j] = i
    # 获取素数表
    primes = [i for i, is_prime in enumerate(is_primes) if i>1 and is_prime is False]
    return is_primes, smallest_primes, primes

def prime_factorization(n, primes, smallest_primes):
    '''获取数n的整数分解'''
    assert n>1
    prime_factor_counter = Counter()
    quotient = n
    while quotient>1:
        prime_factor_counter[smallest_primes[quotient]] += 1
        quotient //= smallest_primes[quotient]
    return prime_factor_counter

# =======广义欧拉算法，modulo求逆，gcd等等=========
def extended_euclidean_algorithm(a, b):
    """return (g, x, y) such that a*x + b*y = g = gcd(a, b) where a,b>=1"""
    x0, x1, y0, y1 = 0, 1, 1, 0
    while a != 0:
        q, b, a = b // a, a, b % a
        y0, y1 = y1, y0 - q * y1
        x0, x1 = x1, x0 - q * x1
    return b, x0, y0

def modulo_inv(a, b):
    """return a^(-1) mod b 
    or x such that  (x * a) == 1 (mod b)
    or x*a + b*y = 1 
    still need a,b>=1"""
    g, x, _ = extended_euclidean_algorithm(a, b)
    if g==1:
        return x % b
    return False

# =====中国剩余定理，多值中国剩余定理，广义中国剩余定理，广义多值中国剩余定理=======
def multivalued_chinese_remainder(n, a, prod = None):
    '''x == (a_i = [...]) (mod n_i) 
        求x
        x = sum(ga_i*s_i*(prod/n_i)), where
        prod = n_1* n_2*..., ga为a_i笛卡尔积的一项，ga_i为ga的第i项，
        s_i = (prod/n_i)^-1 (mod n_i)，程序中
        coeff = s_i*(prod/n_i)
    '''
    prod = reduce(lambda a, b: a*b, n) if prod is None else prod
    for a_i in a:
        if not len(a_i):
            return [], 1

    coeffs = []
    for n_i in n:
        coeff = prod // n_i
        coeffs.append(modulo_inv(coeff, n_i)*coeff)

    xs = [sum([a_perm_i * coeff for a_perm_i, coeff in zip(a_perm, coeffs)])%prod 
        for a_perm in itertools.product(*a)]
    return xs, prod

def extended_multivalued_chinese_remainder(n1, n2, a1, a2, prod = None):
    '''x == (a1 = [...]) (mod n1)
       x == (a2 = [...]) (mod n2) 
       求x
       x = a1*(n2/gcd)*q+a2*(n1/gcd)*p,where n1*p+n2*q = gcd
    '''
    # prod = reduce(lambda a, b: a*b, n) if prod is None else prod
    if len(a1)==0 or len(a2)==0:
        return [], 1
    gcd, p, q = extended_euclidean_algorithm(n1, n2)
    prod = n1*n2//gcd
    coeff1, coeff2 = (n2//gcd)*q, (n1//gcd)*p
    xs = []
    for a1_, a2_ in itertools.product(a1, a2):
        if (a1_-a2_)%gcd == 0:
            xs.append((a1_*coeff1+a2_*coeff2)%prod)
    return xs, prod

# ===========获取二次剩余类=====================
def get_quadratic_residue_table(ps):
    '''ps = [p...] return table[i, p] = (i|p)'''
    table = {(i, p):0 if not i else -1 for p in ps for i in range(p)}
    for p in ps:
        for i in get_quadratic_residue_set(p):
            table[(i, p)] = 1
    return table

def get_quadratic_residue_set(p):
    '''p的二次剩余类'''
    return {(i*i)%p for i in range(1, p)}

def get_quadratic_residue(i, p, table):
    '''取得是否在二次剩余类的值'''
    return table[(i, p)]

# ======判断一个数能整除2的最大幂次方=====
def get_log_2(a):
    if a == 0:
        return 0
    if a<0:
        a = -a
    '''返回一个数的2的最大幂次方'''
    cnt = 0
    while a&1==0:
        a>>=1
        cnt += 1
    return cnt

# =======判断是否是平方根=============
def is_square(a):
    '''if a is a perfect square, return a^0.5, else return False'''
    if a<0:
        return False
    if a==0 or a==1:
        return a
    x = a // 2
    seen = set([x])
    while x * x != a:
        x = (x + (a // x)) // 2
        if x in seen:
            return False
        seen.add(x)
    return x

# ======= 求解z^3 = k (mod p) ================
def cuberoots(k, p):
    '''获取z^3 = k (mod p)的所有解，
    一些特殊情况下，可以简单获得，否则，调用Tonelli-Shanks算法'''
    #Non-trivial solution of x^r=1
    def onemod(p,r):
        t=p-2
        while pow(t,(p-1)//r,p)==1:
            t-=1
        return pow(t,(p-1)//r,p)

    def solution(p,root): 
        g=onemod(p,3)
        return [root%p,(root*g)%p,(root*g*g)%p]

    #---MAIN---
    k=k%p

    if p in [2,3] or k==0:
        return [k]

    if p%3 == 2:
        return [pow(k,(2*p - 1)//3, p)] #Eén oplossing

    #There are 3 or no solutions 

    #No solution
    if pow(k,(p-1)//3,p)>1:
        return []

    if p%9 == 4:                                #[13, 31, 67]
        root = pow(k,(2*p + 1)//9, p)  
        if pow(root,3,p) == k:
            return solution(p,root)        
        else:
            return []

    if p%9 == 7:                                #[7, 43, 61, 79, 97    
        root = pow(k,(p + 2)//9, p)
        if pow(root,3,p) == k:
            return solution(p,root)
        else:
            return []

    if p%27 == 10:                              #[37, 199]
        root = pow(k,(2*p +7)//27, p)         
        h=onemod(p,9)
        for _ in range(9):
            if pow(root,3,p) == k:
                return solution(p,root)                
            root*=h
        return []        

    if p%27 == 19:                              #[19, 73, 127, 181]
        root = pow(k,(p + 8)//27, p)
        h=onemod(p,9)
        for _ in range(9):
            if pow(root,3,p) == k:
                return solution(p,root)
            root*=h
        return []        

    #We need an algorithm for the remaining cases
    return tonelli3(k,p,True)

def tonelli3(k,p,many=False):
    '''assumes p prime, it returns all cube roots of a mod p using Tonelli-Shanks algorithm'''
    def solution(p,root):
        g=p-2
        while pow(g,(p-1)//3,p)==1:
            g-=1  #Non-trivial solution of x^3=1
        g=pow(g,(p-1)//3,p)
        return [root%p,(root*g)%p,(root*g*g)%p]

    #---MAIN---
    k=k%p
    if p in [2,3] or k==0:
        return [k]
    if p%3 == 2:
        return [pow(k,(2*p - 1)//3, p)] #Eén oplossing

    #No solution
    if pow(k,(p-1)//3,p)>1:
        return []

    #p-1=3^s*t
    s=0
    t=p-1
    while t%3==0:
        s+=1
        t//=3

    #Cubic nonresidu b
    b=p-2
    while pow(b,(p-1)//3,p)==1:
        b-=1

    c,r=pow(b,t,p),pow(k,t,p) 
  
    c1,h=pow(c,3^(s-1),p),1    
    c=pow(c,p-2,p) #c=inverse_mod(Integer(c),p)

    for i in range(1,s):
        d=pow(r,3^(s-i-1),p)
        if d==c1:
            h,r=h*c,r*pow(c,3,p)
        elif d!=1:
            h,r=h*pow(c,2,p),r*pow(c,6,p)           
        c=pow(c,3,p)

    if (t-1)%3==0:
        k=(t-1)//3
    else:
        k=(t+1)//3

    r=pow(k,k,p)*h
    if (t-1)%3==0:
        r=pow(r,p-2,p) #r=inverse_mod(Integer(r),p)

    if pow(r,3,p)==k: 
        if many: 
            return solution(p,r)
        else:
            return [r]
    else: return []

# ======= 求解z^3 = k (mod p^r) ================

def hensel(p, k, r, zs_lower):
    '''r>1; zs_lower is all solutions to z^3=k mod p^{r-1}; return all z for z^3 = k mod p^r = 0'''
    zs = []
    for z_lower in zs_lower:
        zs.extend(_hensel(p, k, r, z_lower))
    return zs

def _hensel(p, k, r, z_lower):
    '''
    r>1; z_lower a solution to z^3=k mod p^{r-1}; return all z for z^3 = k mod p^r
    Hensel's Lemma: 
    let f(x) be a polynomial with integer coefficients, p a prime number, and k>=2, 
    and consider equation f(x) mod p^k = 0. Suppose r is a solution to f(x) mod p^{k-1} = 0. Then:
    I) if f'(r) mod p != 0, then there is a unique integer t, with 0<=t<p, such that 
        f(r+tp^(k-1)) mod p^k = 0, given by t=-(f'(r)^(-1))(f(r)/p^(k-1)) mod p.
    II) if f'(r) mod p = 0 and f(r) mod p^k = 0, then f(r+tp^(k-1)) mod p^k = 0 for all integers t.
    III) if f'(r) mod p = 0 and f(r) mod p^k != 0, then f(x) mod p^k = 0 has no solutions with x = r mod p^(k-1).
    '''
    fr = z_lower**3-k
    dfr = 3*z_lower**2
    if dfr%p != 0:
        t= (-modulo_inv(dfr, p)*(fr//p**(r-1)))%p
        zs = [z_lower+t*p**(r-1)]
    elif fr % p**r == 0:
        zs = [z_lower+t*p**(r-1) for t in range(p)]
    else:
        zs = []
    return zs

# ======== main和test ==================
def main(B, k):
    cons = Constants(B, k)
    for d in range(1, cons.dB):
        if DSlice(cons, d).run():
            break

if __name__ == "__main__":
    main(B = int(2e5), k = 39)
    # main(B = int(5e6), k = 75)
    # main(B = int(100), k = 78)
