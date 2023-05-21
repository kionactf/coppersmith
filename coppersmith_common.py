from sage.all import *

import time

from logger import logger

from lll import do_lattice_reduction, FPLLL, FPLLL_BKZ, FLATTER, NTL, NTL_BKZ

# for large N, precision is problematic
RRh = RealField(prec=4096)


### common function for Coppersmith method variants


def shiftpoly(basepoly, baseidx, Nidx, varsidx_lst):
    N = basepoly.parent().characteristic()
    basepoly_ZZ = basepoly.change_ring(ZZ)
    vars_ZZ = basepoly_ZZ.parent().gens()
    if len(vars_ZZ) != len(varsidx_lst):
        raise ValueError("varsidx_lst len is invalid (on shiftpoly)")
    return (basepoly_ZZ ** baseidx) * (N ** Nidx) * prod([v ** i for v, i in zip(vars_ZZ, varsidx_lst)])


def monomialset(fis):
    m_set = set()
    for fi in fis:
        m_set = m_set.union(set(fi.monomials()))
    return m_set


def genmatrix_from_shiftpolys(shiftpolys, bounds):
    m_lst = list(monomialset(shiftpolys))
    vars_ZZ = shiftpolys[0].parent().gens()
    if len(vars_ZZ) != len(bounds):
        raise ValueError("bounds len is invalid (on genmatrix_from_shiftpolys)")
    matele = []
    for sftpol in shiftpolys:
        sftpol_sub_bound = sftpol.subs({vars_ZZ[i]: vars_ZZ[i]*bounds[i] for i in range(len(vars_ZZ))})
        matele += [sftpol_sub_bound.monomial_coefficient(m_lst[i]) for i in range(len(m_lst))]
    mat = matrix(ZZ, len(matele)//len(m_lst), len(m_lst), matele)
    return mat


def do_LLL(mat):
    #lll, trans = do_lattice_reduction(mat, algorithm=FLATTER, use_pari_kernel=True)
    lll, trans = do_lattice_reduction(mat, algorithm=FPLLL)

    return lll, trans


def filter_LLLresult_coppersmith(basepoly, beta, t, shiftpolys, lll, trans):
    N = basepoly.parent().characteristic()
    howgrave_bound = (RRh(N)**RRh(beta))**RRh(t)
    if len(shiftpolys) != lll.nrows() or len(shiftpolys) != trans.nrows():
        raise ValueError("lll or trans result is invalid (on filter_LLLresult_coppersmith)")
    # use vector (not use matrix norm, but vector norm)
    lll_vec = lll.rows()
    result = []
    for i, lll_vecele in enumerate(lll_vec):
        if all([int(lll_vecele_ele) == 0 for lll_vecele_ele in lll_vecele]):
            continue
        lll_l1norm = lll_vecele.norm(p=1)
        if lll_l1norm >= howgrave_bound:
            continue
        howgrave_ratio = int(((lll_l1norm/howgrave_bound)*(10**15))*(0.1**15))
        logger.debug("lll_l1norm/howgrave_bound: %s", str(howgrave_ratio) )
        pol = 0
        for j, shiftpolysele in enumerate(shiftpolys):
            pol += trans[i,j] * shiftpolysele
        result.append(pol)
    return result
