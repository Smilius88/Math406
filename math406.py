# Returns a list of all primes less than limit, using a sieve.
def primes(limit):
    Primes = []
    isp = [False, False, True] + [True] * (limit-1)
    for i in xrange(3, limit+1, 2):
        if isp[i]:
            Primes.append(i)
            for j in xrange(i+i, limit, i):
                isp[j] = false
    return Primes


# Returns a tuple (x, y, m) such that ax+by == m, where m = gcd(a,b). 
# Implements the exteded euclidean algorithm. 
# Uncomment the two "print" lines and recompile to show all steps.
def extended_euclidean(a, b):
#	print (1, 0, a)
    def ee(line1, line2):
        c11, c12, r1 = line1
        c21, c22, r2 = line2
        mul, remainder = divmod(r1, r2)
        # print line2
        if remainder == 0:
            return line2
        else:
            return ee(line2, (c11-c21*mul, c12-c22*mul, remainder))
    return ee((1, 0, a), (0, 1, b))


# Returns an integer x, such that x == a mod b, for all tuples (a, b) in Moduli.
# Accepts a list of tuples, Moduli, such that x == a mod b for each one.
def chinese_remainder(Moduli):

    rems, mods = zip(*Moduli)
    total = 0
    N = reduce(int.__mul__, mods)
    for rem, mod in Moduli:
        n = N/mod
        a1, _, _ = extended_euclidean(n, mod)
        total += n*a1*rem
    return (total % N), N


#Polynomial is a class for polynomials (surprise!). It is very basic.
class Polynomial:

    def __init__(self, Codict):
        self.Coeffs = {}
        for key, value in Codict.iteritems():
            self.Coeffs[key] = value
        self.degree = max(Codict.keys())

    def derive(self):
        Derivative = {}
        for a in self.Coeffs.keys():
            if a > 0:
                Derivative[a-1] = a*self.Coeffs[a]
        return Polynomial(Derivative)

    def evaluate(self, x, mod=0):
        ans, temp = 0,  0
        for a in self.Coeffs.keys():
            temp = (x**a)*self.Coeffs[a]
            if mod > 0:
                temp %= mod
            ans += temp
        if mod > 0:
            ans %= mod
        return ans

# Returns a list of numbers n such that f(n) == 0 (mod p^k).
# Only works for primes. For arbitrary moduli, use the (slower) solve_poly_con
# Accepts a Polynomial (Funct), and tuple (p, k), representing p^kth power (Pmod).
# Implements Hensel's lemma. 
def lift(Funct, Pmod):
    prime, power = Pmod
    Derivative = Funct.derive()
    Solutions = []

    if power == 1:
        for i in xrange(prime):
            if (Funct.evaluate(i, prime) == 0):
                Solutions.append(i)

    else:
        Candidates = lift(Funct, (prime, power-1))
        for candidate in Candidates:
            f = Funct.evaluate(candidate, prime**power)
            fprime = Derivative.evaluate(candidate, prime)
            if fprime != 0:
                inv, _, _ = extended_euclidean(Derivative.evaluate(candidate), prime)
                t = -1*(inv*(Funct.evaluate(candidate)/(prime**(power-1)))) % prime
                soln = candidate+t*(prime**(power-1)) % (prime**power)
                Solutions.append(soln)
            elif f == 0 and fprime == 0:
                p, soln = prime**(power-1), candidate
                for i in xrange(prime):
                    Solutions.append(soln+ p*i)
    return Solutions
        

# returns a defaultdict of factors of the number N
def factor(N):
    from math import sqrt
    from collections import defaultdict
    bound, fact = int(sqrt(N))+1, 2
    factors = defaultdict(int)
    while N > 1:
            if N % fact == 0:
                factors[fact] += 1
                N /= fact
            else:
                fact += 1
    return factors

# Returns a list of all integers x such that f(x) == 0 (mod z)
# Accepts a dict of coefficients and a modulus.
# Solves a polynomial congruence for arbitrary numbers and polynomials.
#
# At some point I ran out of descriptive variable names. Guess where!
def solve_poly_con(coeffs, modulus):
    from itertools import repeat, product

    r = []
    p = Polynomial(coeffs)
    factors = factor(modulus)
    congruences = []
    for prime, power in factors.iteritems():
        mods = lift(p, (prime, power))
        mods = zip(mods, repeat(prime**power))
        congruences.append(mods)
    grist = product(*congruences)
    for l in grist:
        c, m = chinese_remainder(l)
        r.append(c)
    return r

