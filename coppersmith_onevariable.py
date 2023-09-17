from sage.all import *

import time

from coppersmith_common import RRh, shiftpoly, genmatrix_from_shiftpolys, do_LLL, filter_LLLresult_coppersmith
from rootfind_ZZ import rootfind_ZZ
from logger import logger


### one variable coppersmith
def coppersmith_one_var_core(basepoly, bounds, beta, t, u, delta, **lllopt):
    logger.info("trying param: beta=%f, t=%d, u=%d, delta=%d", beta, t, u, delta)
    basepoly_vars = basepoly.parent().gens()
    basepoly = basepoly / basepoly.monomial_coefficient(basepoly_vars[0] ** delta)

    shiftpolys = []
    for i in range(u-1, -1, -1):
        # x^i * f(x)^t
        shiftpolys.append(shiftpoly(basepoly, t, 0, [i]))
    for i in range(1, t+1, 1):
        for j in range(delta-1, -1, -1):
            # x^j * f(x)^(t-i) * N^i
            shiftpolys.append(shiftpoly(basepoly, t-i, i, [j]))

    mat = genmatrix_from_shiftpolys(shiftpolys, bounds)
    lll, trans = do_LLL(mat, **lllopt)
    result = filter_LLLresult_coppersmith(basepoly, beta, t, shiftpolys, lll, trans)
    return result


def coppersmith_onevariable(basepoly, bounds, beta, maxmatsize=100, maxu=8, **lllopt):
    if type(bounds) not in [list, tuple]:
        bounds = [bounds]

    N = basepoly.parent().characteristic()

    basepoly_vars = basepoly.parent().gens()
    if len(basepoly_vars) != 1:
        raise ValueError("not one variable poly")
    try:
        delta = basepoly.weighted_degree([1])
    except:
        delta = basepoly.degree()

    log_N_X = RRh(log(bounds[0], N))
    if log_N_X >= RRh(beta)**2/delta:
        raise ValueError("too much large bound")

    testimate = int(1/(((RRh(beta)**2)/delta)/log_N_X - 1))//2

    logger.debug("testimate: %d", testimate)
    t = min([maxmatsize//delta, max(testimate, 1)])

    whole_st = time.time()

    curfoundpols = []
    while True:
        if t*delta > maxmatsize:
            raise ValueError(f"maxmatsize exceeded: {t*delta}")
        u0 = max([int((t+1)/RRh(beta) - t*delta), 0])
        for u_diff in range(0, maxu+1):
            u = u0 + u_diff
            if t*delta + u > maxmatsize:
                break
            foundpols = coppersmith_one_var_core(basepoly, bounds, beta, t, u, delta, **lllopt)
            if len(foundpols) == 0:
                continue

            curfoundpols += foundpols
            curfoundpols = list(set(curfoundpols))
            sol = rootfind_ZZ(curfoundpols, bounds)
            if sol != [] and sol is not None:
                whole_ed = time.time()
                logger.info("whole elapsed time: %f", whole_ed-whole_st)
                return sol
            elif len(curfoundpols) >= 2:
                whole_ed = time.time()
                logger.warning(f"failed. maybe, wrong pol was passed.")
                logger.info("whole elapsed time: %f", whole_ed-whole_st)
                return []
        t += 1
    # never reached here
    return None
