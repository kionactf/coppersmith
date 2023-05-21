from sage.all import *

import time
import itertools

from coppersmith_common import RRh, shiftpoly, genmatrix_from_shiftpolys, do_LLL, filter_LLLresult_coppersmith
from rootfind_ZZ import rootfind_ZZ
from logger import logger


def gen_set_leading_monomials(basepoly):
    if basepoly.is_constant():
        return [basepoly.parent()(1)]

    lmset = [basepoly.parent()(1)]
    for monomial in basepoly.monomials():
        newlmfound = True
        newlmset = []
        for lmsetele in lmset:
            if monomial % lmsetele != 0:
                newlmset.append(lmsetele)
            if lmsetele % monomial == 0:
                newlmfound = False
        if newlmfound:
            newlmset.append(monomial)
        lmset = newlmset[:]
    return lmset


def generate_M_with_ExtendedStrategy(basepoly, lm, t, d):
    basepoly_vars = basepoly.parent().gens()
    n = len(basepoly_vars)

    M = {}
    basepoly_pow_monos = (basepoly ** t).monomials()
    for k in range(t+1):
        M[k] = set()
        basepoly_powk_monos = (basepoly ** (t - k)).monomials()
        for monos in basepoly_pow_monos:
            if monos // (lm ** k) in basepoly_powk_monos:
                for extra in itertools.product(range(d), repeat=n):
                    g = monos * prod([v ** i for v, i in zip(basepoly_vars, extra)])
                    M[k].add(g)
    M[t+1] = set()
    return M


### multivariate coppersmith with some heuristic (jochemsz-may)
def coppersmith_multivariate_heuristic_core(basepoly, bounds, beta, t, d, lm, maxmatsize=100):
    logger.info("trying param: beta=%f, t=%d, d=%d, lm=%s", beta, t, d, str(lm))
    basepoly_vars = basepoly.parent().gens()
    n = len(basepoly_vars)

    basepoly_i = basepoly / basepoly.monomial_coefficient(lm)

    M = generate_M_with_ExtendedStrategy(basepoly_i, lm, t, d)
    shiftpolys = []
    for k in range(t+1):
        for mono in M[k] - M[k+1]:
            curmono = (mono // (lm ** k))
            xi_idx = curmono.exponents()[0]
            shiftpolys.append(shiftpoly(basepoly_i, k, t - k, xi_idx))

    mat = genmatrix_from_shiftpolys(shiftpolys, bounds)
    if mat.ncols() > maxmatsize:
        logger.warning("maxmatsize exceeded(on coppersmith_multivariate_heuristic)")
        return []

    lll, trans = do_LLL(mat)
    result = filter_LLLresult_coppersmith(basepoly, beta, t, shiftpolys, lll, trans)
    return result


def coppersmith_multivariate_heuristic(basepoly, bounds, beta, maxmatsize=100, maxd=8):
    if type(bounds) not in [list, tuple]:
        raise ValueError("bounds should be list or tuple")

    N = basepoly.parent().characteristic()

    basepoly_vars = basepoly.parent().gens()
    n = len(basepoly_vars)
    if n == 1:
        raise ValueError("one variable poly")

    # dealing with all candidates of leading monomials
    lms = gen_set_leading_monomials(basepoly)

    t = 2

    whole_st = time.time()

    curfoundpols = []
    while True:
        d0 = t
        for d_diff in range(0, maxd+1):
            d = d0 + d_diff
            for lm in lms:
                foundpols = coppersmith_multivariate_heuristic_core(basepoly, bounds, beta, t, d, lm, maxmatsize=maxmatsize)
                if len(foundpols) == 0:
                    continue
                curfoundpols += foundpols
                curfoundpols = list(set(curfoundpols))
                sol = rootfind_ZZ(curfoundpols, bounds)
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
