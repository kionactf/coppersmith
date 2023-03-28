from sage.all import *

import time
from Crypto.Util.number import *
from random import randrange, randint

from coppersmith_onevariable import coppersmith_onevariable
from coppersmith_linear import coppersmith_linear
from coppersmith_multivariate_heuristic import coppersmith_multivariate_heuristic
from logger import logger


def example_onevariable_linear():
    our_small_roots = lambda f, X: coppersmith_onevariable(f, [X], beta)
    bitsize = 2048
    while True:
        p = getPrime(bitsize//2)
        q = getPrime(bitsize//2)
        if p != q:
            break
    N = p*q

    beta = (1.0*(bitsize//2-1))/bitsize # (1024-1)/2048 >= 0.4995 (also 511/1024 >=0.499)

    # in my exp, worked for 500, but took much time...
    discardbitsizelst = [40, 256, 496]
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
        if discardbitsize >= 490:
            epsilon = 0.01
        else:
            epsilon = 0.05

        sage_st = time.time()
        print(f"sage result:{f.small_roots(X=2**discardbitsize, beta=beta, epsilon=epsilon)}")
        sage_ed = time.time()
        logger.debug("sage comp elapsed time: %f", sage_ed - sage_st)

        # sometimes works with small beta (but 496 not works)
        print(f"sage result (small beta):{f.small_roots(X=2**discardbitsize, beta=0.4)}")


def example_twovariable_linear():
    our_small_roots = lambda f, bounds: coppersmith_linear(f, bounds, beta)
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
    our_small_roots = lambda f, bounds: coppersmith_linear(f, bounds, beta)
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
    sol = coppersmith_multivariate_heuristic(f, bounds, 1.0)
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
    sol = coppersmith_multivariate_heuristic(f, bounds, 1.0)
    print(f"result:{sol}, real:{roots}")


def _example_multivariate_heuristic_3():
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

    sol = coppersmith_multivariate_heuristic(f, bounds, 1.0)
    print(f"result:{sol}, real:{(str(ans1)[L:], str(ans2)[L:])}")


def example_multivariate_heuristic():
    _example_multivariate_heuristic_1()
    _example_multivariate_heuristic_2()
    _example_multivariate_heuristic_3()


#example_onevariable_linear()
#example_twovariable_linear()
#example_threevariable_linear()
example_multivariate_heuristic()
