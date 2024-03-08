from sage.all import *

import time
import itertools
import traceback

from coppersmith_common import RRh, shiftpoly, genmatrix_from_shiftpolys, do_LLL, filter_LLLresult_coppersmith
from rootfind_ZZ import rootfind_ZZ
from contextclass import context
from logger import *


### multivariate linear coppersmith (herrmann-may)
def coppersmith_linear_core(basepoly, bounds, beta, t, m):
    logger.info("trying param: beta=%f, t=%d, m=%d", beta, t, m)
    basepoly_vars = basepoly.parent().gens()
    n = len(basepoly_vars)

    shiftpolys = []
    for i, basepoly_var in enumerate(basepoly_vars):
        try:
            # NOTE: for working not integral domain, write as basepoly * (1/val) (do not write as basepoly / val)
            basepoly_i = basepoly * (1 / basepoly.monomial_coefficient(basepoly_var))
        except:
            traceback.print_exc()
            logger.warning(f"maybe, facing an non-invertible monomial coefficient({basepoly_var}). still continuing...")
            ## TODO: continue procedures with hoping some monomial coefficients for at least one variable is invertible
            continue

        for k in range(m+1):
            for j in range(m-k+1):
                for xi_idx_sub in itertools.combinations_with_replacement(range(n-1), j):
                    xi_idx = [xi_idx_sub.count(l) for l in range(n-1)]
                    assert sum(xi_idx) == j
                    xi_idx.insert(i, 0)
                    # x2^i2 * ... * xn^in * f^k * N^max(t-k,0)
                    shiftpolys.append(shiftpoly(basepoly_i, k, max(t-k, 0), xi_idx))

    mat, m_lst = genmatrix_from_shiftpolys(shiftpolys, bounds)
    lll, _ = do_LLL(mat)
    result = filter_LLLresult_coppersmith(basepoly, beta, t, m_lst, lll, bounds)
    return result


def coppersmith_linear(basepoly, bounds, beta, maxmatsize=100, maxm=8):
    if type(bounds) not in [list, tuple]:
        raise ValueError("bounds should be list or tuple")

    if beta >= 1.0:
        raise ValueError("beta is invalid. (for beta=1.0, use normal lattice reduction method directly.)")

    N = basepoly.parent().characteristic()

    basepoly_vars = basepoly.parent().gens()
    n = len(basepoly_vars)
    if n == 1:
        raise ValueError("one variable poly")

    if not set(basepoly.monomials()).issubset(set(list(basepoly_vars)+[1])):
        raise ValueError("non linear poly")

    log_N_X = RRh(log(product(bounds), N))
    log_N_X_bound = 1-(1-RRh(beta))**(RRh(n+1)/n) - (n+1)*(1-(1-RRh(beta))**(RRh(1)/n)) * (1-RRh(beta))

    if log_N_X >= log_N_X_bound:
        raise ValueError("too much large bound")

    mestimate = (n*(-RRh(beta)*ln(1-beta) + ((1-RRh(beta))**(-0.278465))/pi)/(log_N_X_bound - log_N_X))/(n+1.5)
    tau = 1 - (1-RRh(beta))**(RRh(1)/n)
    testimate = int(mestimate * tau + 0.5)

    logger.debug("testimate: %d", testimate)
    t = max(testimate, 1)

    while True:
        if t == 1:
            break
        m = int(t/tau+0.5)
        if binomial(n+1+m-1, m) <= maxmatsize:
            break
        t -= 1

    whole_st = time.time()

    curfoundpols = []
    while True:
        m0 = int(t/tau+0.5)
        if binomial(n+1+m0-1, m0) > maxmatsize:
            raise ValueError(f"maxmatsize exceeded: {binomial(n+1+m0-1, m0)}")
        for m_diff in range(0, maxm+1):
            m = m0 + m_diff
            if binomial(n+1+m-1, m) > maxmatsize:
                break
            foundpols = coppersmith_linear_core(basepoly, bounds, beta, t, m)
            if len(foundpols) == 0:
                continue

            curfoundpols += foundpols
            curfoundpols = list(set(curfoundpols))
            sol = rootfind_ZZ(curfoundpols, bounds, **context.rootfindZZopt)
            if sol != [] and sol is not None:
                whole_ed = time.time()
                logger.info("whole elapsed time: %f", whole_ed-whole_st)
                return sol

            polrate = (1.0 * len(curfoundpols))/n
            if polrate > 1.0:
                logger.warning(f"polrate is over 1.0 (you might have inputted wrong pol): {polrate}")
                whole_ed = time.time()
                logger.info("whole elapsed time (not ended): %f", whole_ed-whole_st)
        t += 1
    # never reached here
    return None
