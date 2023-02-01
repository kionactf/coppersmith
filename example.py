from sage.all import *

import time
from Crypto.Util.number import *

from coppersmith_onevariable import coppersmith_onevariable
from coppersmith_linear import coppersmith_linear
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

#example_onevariable_linear()
#example_twovariable_linear()
example_threevariable_linear()