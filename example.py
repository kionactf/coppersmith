from sage.all import *

import time
from Crypto.Util.number import *
from random import randrange, randint, choices as random_choices
import string
import sys

from coppersmith_onevariable import coppersmith_onevariable
from coppersmith_linear import coppersmith_linear
from coppersmith_multivariate_heuristic import coppersmith_multivariate_heuristic
from lll import *
from logger import *


sys.set_int_max_str_digits(8000)


lllopt = {}
#lllopt = {'algorithm':FPLLL}
#lllopt = {'algorithm':FLATTER, 'use_pari_kernel':True}
#lllopt = {'algorithm':FPLLL_BKZ, 'blocksize':3}


def example_onevariable_linear():
    our_small_roots = lambda f, X: coppersmith_onevariable(f, [X], beta, **lllopt)
    bitsize = 2048
    while True:
        p = getPrime(bitsize//2)
        q = getPrime(bitsize//2)
        if p != q:
            break
    N = p*q

    beta = (1.0*(bitsize//2-1))/bitsize # (1024-1)/2048 >= 0.4995 (also 511/1024 >=0.499)

    # in my exp, worked for 500, but took much time...
    discardbitsizelst = [40, 256, 488]
    for discardbitsize in discardbitsizelst:
        print(discardbitsize)
        p0 = p>>discardbitsize
        P = PolynomialRing(Zmod(N), 'x')
        ## following P works for our_small_roots, but not Sagemath small_roots
        #P = PolynomialRing(Zmod(N), 1, 'x')
        x = P.gens()[0]
        f = (p0 << discardbitsize) + x
        result = our_small_roots(f, 2**discardbitsize)
        print(f"result:{result}, real:{p % (2**discardbitsize)}")

        # check output of sage small_roots
        if discardbitsize >= 488:
            epsilon = 0.02
        else:
            epsilon = 0.05

        sage_st = time.time()
        print(f"sage result:{f.small_roots(X=2**discardbitsize, beta=beta, epsilon=epsilon)}")
        sage_ed = time.time()
        logger.debug("sage comp elapsed time: %f", sage_ed - sage_st)

        # sometimes works with small beta (but 488 not works)
        print(f"sage result (small beta):{f.small_roots(X=2**discardbitsize, beta=0.1)}")


def example_twovariable_linear():
    our_small_roots = lambda f, bounds: coppersmith_linear(f, bounds, beta, **lllopt)
    bitsize = 2048
    while True:
        p = getPrime(bitsize//2)
        q = getPrime(bitsize//2)
        if p != q:
            break
    N = p*q

    beta = (1.0*(bitsize//2-1))/bitsize # (1024-1)/2048 >= 0.4995 (also 511/1024 >=0.499)

    # it seems severe for over (160, 160) (t=3 needs much time)
    discardbitpointlst = [[(756, 20),(256, 20)], [(756, 135), (256, 135)], [(756, 160), (256, 160)]]
    for discardbitpoint in discardbitpointlst:
        print(discardbitpoint)
        p0 = 0
        real = []
        for i, ele in enumerate(discardbitpoint):
            if i == 0:
                p0 += (p>>ele[0]) << ele[0]
            else:
                p0 += ((p % (2**(discardbitpoint[i-1][0]-discardbitpoint[i-1][1]))) >> ele[0]) << ele[0]
            real.append((p % (2**ele[0]))>>(ele[0]-ele[1]))
        p0 += p % (2**(discardbitpoint[-1][0] - discardbitpoint[-1][1]))

        P = PolynomialRing(Zmod(N), 2, 'xy')
        P_vars = P.gens()
        bounds = []
        f = p0
        for i,ele in enumerate(discardbitpoint):
            f += (2**(ele[0]-ele[1]))*P_vars[i]
            bounds.append(2**ele[1])
        result = our_small_roots(f, bounds)
        print(f"result:{result}, real:{real}")


def example_threevariable_linear():
    our_small_roots = lambda f, bounds: coppersmith_linear(f, bounds, beta, **lllopt)
    bitsize = 2048
    while True:
        p = getPrime(bitsize//2)
        q = getPrime(bitsize//2)
        if p != q:
            break
    N = p*q

    beta = (1.0*(bitsize//2-1))/bitsize # (1024-1)/2048 >= 0.4995 (also 511/1024 >=0.499)

    # it seems severe for over (160, 160) (t=3 needs much time)
    discardbitpointlst = [[(756, 20),(512, 20),(256, 20)], [(756, 40),(512, 40),(256,40)], [(756, 72), (512,72), (256, 72)]]
    for discardbitpoint in discardbitpointlst:
        print(discardbitpoint)
        p0 = 0
        real = []
        for i, ele in enumerate(discardbitpoint):
            if i == 0:
                p0 += (p>>ele[0]) << ele[0]
            else:
                p0 += ((p % (2**(discardbitpoint[i-1][0]-discardbitpoint[i-1][1]))) >> ele[0]) << ele[0]
            real.append((p % (2**ele[0]))>>(ele[0]-ele[1]))
        p0 += p % (2**(discardbitpoint[-1][0] - discardbitpoint[-1][1]))

        P = PolynomialRing(Zmod(N), 3, 'xyz')
        P_vars = P.gens()
        bounds = []
        f = p0
        for i,ele in enumerate(discardbitpoint):
            f += (2**(ele[0]-ele[1]))*P_vars[i]
            bounds.append(2**ele[1])
        result = our_small_roots(f, bounds)
        print(f"result:{result}, real:{real}")


def example_shortpad_attack():
    # example of Coppersmith's short-pad attack; non-monic univariate polynomial case
    bitsize = 2048
    padbytelen = 24
    while True:
        p = getPrime(bitsize//2)
        q = getPrime(bitsize//2)
        N = p * q
        e = 3
        phi = (p - 1) * (q - 1)
        if GCD(phi, e) == 1:
            d = pow(e, -1, phi)
            break
    charlist = string.ascii_uppercase + string.ascii_lowercase + string.digits
    M = ''.join(random_choices(charlist, k=115)) + '_' + ''.join(random_choices(charlist, k=115))
    pad = ''.join(random_choices(charlist, k=padbytelen))

    M_1 = bytes_to_long((M + '\x00' * padbytelen).encode())
    M_2 = bytes_to_long((M + pad).encode())

    C_1 = pow(M_1, e, N)
    C_2 = pow(M_2, e, N)

    # attack from here
    P_first = PolynomialRing(ZZ, 2, "xy")
    x, y = P_first.gens()

    ## x = (M + '\x00' * padbytelen), y = pad
    pol1 = x ** e - C_1
    pol2 = (x + y) ** e - C_2
    pol = pol1.resultant(pol2, x)

    pol_uni = pol.univariate_polynomial().change_ring(Zmod(N))
    sol = coppersmith_onevariable(pol_uni, [2**(8*padbytelen)], 1.0, **lllopt)[0]

    ## Franklin-Reiter related-message attack
    pol1_uni = pol1.univariate_polynomial().change_ring(Zmod(N))
    pol2_uni = pol2.subs({x:x, y:sol}).univariate_polynomial().change_ring(Zmod(N))

    def composite_gcd(f1, f2):
        if f2 == 0:
            return f1.monic()
        if f1.degree() < f2.degree():
            return composite_gcd(f2, f1)
        return composite_gcd(f2, f1 % f2)

    pol_gcd = composite_gcd(pol1_uni, pol2_uni)
    assert pol_gcd.degree() == 1

    degoneinv = (pol_gcd.monomial_coefficient(pol_gcd.parent().gens()[0]) ** (-1))
    found_M_N = -pol_gcd.constant_coefficient() * degoneinv
    found_M = long_to_bytes(int(found_M_N.lift())).split(b'\x00')[0]

    print(f"result:{found_M}, real:{M}")


def example_chronophobia():
    # chronophobia from idekCTF2022
    L = 200

    p = getPrime(512)
    q = getPrime(512)
    n = p*q
    phi = (p-1) * (q-1)

    t = randint(0, n-1)
    d = randint(128, 256)
    r = pow(2, 1 << d, phi)

    ans1 = pow(t, r, n)
    u1 = int(str(ans1)[:L])
    L1down = len(str(ans1)[L:])
    ans2 = pow(pow(t, 2, n), r, n)
    u2 = int(str(ans2)[:L])
    L2down = len(str(ans2)[L:])

    P = PolynomialRing(Zmod(n), 2, ["x", "y"])
    x, y = P.gens()

    f = (u1 * (10**L1down) + x)**2 - (u2 * (10**L2down) + y)
    bounds = [10**L1down, 10**L2down]

    sol = coppersmith_multivariate_heuristic(f, bounds, 1.0, **lllopt)
    print(f"result:{sol}, real:{(int(str(ans1)[L:]), int(str(ans2)[L:]))}")


def example_bivariate_stereotyped_message_attack():
    # @Warri posted to cryptohack discord channel (#cryptography, May.23, 2023)
    bitsize = 1024
    part_M_first_size = 14
    part_M_second_size = 15
    while True:
        p = getPrime(bitsize//2)
        q = getPrime(bitsize//2)
        N = p * q
        e = 3
        phi = (p - 1) * (q - 1)
        if GCD(phi, e) == 1:
            d = pow(e, -1, phi)
            break
    charlist = string.ascii_uppercase + string.ascii_lowercase + string.digits
    part_M_first = ''.join(random_choices(charlist, k=part_M_first_size))
    part_M_second = ''.join(random_choices(charlist, k=part_M_second_size))
    prefix = ''.join(random_choices(charlist, k=40))
    midfix = ''.join(random_choices(charlist, k=30))
    suffix = ''.join(random_choices(charlist, k=20))

    M = bytes_to_long(
        (prefix + part_M_first + midfix + part_M_second + suffix).encode()
    )
    C = pow(M, e, N)

    # attack from here
    P = PolynomialRing(Zmod(N), 2, "xy")
    x, y = P.gens()

    f_p = bytes_to_long(suffix.encode())
    f_p += y * (2**(8*len(suffix)))
    f_p += bytes_to_long(midfix.encode()) * (2**(8*(part_M_second_size + len(suffix))))
    f_p += x * (2**(8*(len(midfix) + part_M_second_size + len(suffix))))
    f_p += bytes_to_long(prefix.encode()) * (2**(8*(part_M_first_size + len(midfix) + part_M_second_size + len(suffix))))
    f = f_p ** e - C

    bounds = (2**(8*part_M_first_size), 2**(8*part_M_second_size))
    sol = coppersmith_multivariate_heuristic(f, bounds, 1.0, **lllopt)

    found_part_M_first = long_to_bytes(int(sol[0][0]))
    found_part_M_second = long_to_bytes(int(sol[0][1]))
    print(f"result:{(found_part_M_first, found_part_M_second)}, real:{(part_M_first, part_M_second)}")


def _example_multivariate_heuristic_1():
    # from bivariate_example on https://github.com/josephsurin/lattice-based-cryptanalysis/blob/main/examples/problems/small_roots.sage
    N = random_prime(2**512) * random_prime(2**512)
    bounds = (2**164, 2**164) # N**0.16
    roots = tuple(randrange(bound) for bound in bounds)
    P = PolynomialRing(Zmod(N), 2, ["x", "y"])
    x, y = P.gens()
    monomials = [x, y, x*y, x**2, y**2]
    f = sum(randrange(N) * monomial for monomial in monomials)
    f -= f(*roots)
    sol = coppersmith_multivariate_heuristic(f, bounds, 1.0, **lllopt)
    print(f"result:{sol}, real:{roots}")


def _example_multivariate_heuristic_2():
    # from trivariate_example on https://github.com/defund/coppersmith/blob/master/examples.sage
    p = random_prime(2**1024)
    q = random_prime(2**1024)
    N = p*q
    bounds = (2**246, 2**246, 2**246) # N**0.12
    roots = tuple(randrange(bound) for bound in bounds)
    P = PolynomialRing(Zmod(N), 3, ["x", "y", "z"])
    x, y, z = P.gens()
    monomials = [x, y, x*y, x*z, y*z]
    f = sum(randrange(N)*monomial for monomial in monomials)
    f -= f(*roots)
    sol = coppersmith_multivariate_heuristic(f, bounds, 1.0, **lllopt)
    print(f"result:{sol}, real:{roots}")




def example_multivariate_heuristic():
    _example_multivariate_heuristic_1()
    _example_multivariate_heuristic_2()


if __name__ == '__main__':
    example_onevariable_linear()
    example_twovariable_linear()
    example_threevariable_linear()
    example_shortpad_attack()
    example_multivariate_heuristic()
    example_chronophobia()
    example_bivariate_stereotyped_message_attack()
