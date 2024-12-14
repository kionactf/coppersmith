from sage.all import *

from typing import Tuple, List
from subprocess import run as subprocess_run
from re import sub as re_sub
import os
import time
import traceback

# NTL (some code is copied from https://github.com/sagemath/sage/blob/develop/src/sage/matrix/matrix_integer_dense.pyx)
import sage.libs.ntl.all
import sage.libs.ntl.ntl_mat_ZZ
ntl_ZZ = sage.libs.ntl.all.ZZ
ntl_mat = lambda A_: sage.libs.ntl.ntl_mat_ZZ.ntl_mat_ZZ(A_.nrows(), A_.ncols(), [ntl_ZZ(z) for z in A_.list()])

from fpylll import IntegerMatrix, GSO, Pruning, Enumeration

from logger import *


_coppersmith_dir = os.path.dirname(__file__)
_fplll_path = os.path.join(_coppersmith_dir, 'fplll', 'fplll') # /usr/bin
_flatter_path = os.path.join(_coppersmith_dir, 'flatter', 'build', 'bin') # /usr/bin

fplll_path = os.environ.get('COPPERSMITHFPLLLPATH', _fplll_path)
flatter_path = os.environ.get('COPPERSMITHFLATTERPATH', _flatter_path)

# algorithm
FPLLL = 0
FPLLL_BKZ = 1
FLATTER = 2
NTL = 3
NTL_BKZ = 4


# fplll option
## fplll_version ('fast' is only double, 'proved' cannot be used with early reduction)
WRAPPER = 'wrapper'
HEURISTIC = 'heuristic'

# on flatter, call kermat on pari
pari.allocatemem(1024*1024*1024)


def _from_sagematrix_to_fplllmatrix(mat: matrix) -> str:
    return '[' + re_sub(
        r'\[ ',
        r'[',
        re_sub(r' +', r' ', str(mat))
    ) + ']'


def _fplllmatrix_to_sagematrix(matrixstr: str) -> matrix:
    matlist = eval(matrixstr.replace(' ', ',').replace('\n', ','))
    return matrix(ZZ, matlist)


def _transformation_matrix(mat, lllmat, use_pari_matsol=False):
    # pari.matker() does not assure smallest kernel in Z (seems not call hermite normal form)
    # Sage kernel calls hermite normal form
    #
    # for computing ZZ transformation, use pari.matker, pari.matsolvemod
    # assume first kerdim vectors for lllmat are zero vector
    #
    # anyway, transformation computation after LLL/BKZ is slow.
    # instead, use builtin transformation computation on LLL/BKZ package

    if use_pari_matsol:
        mat_pari = pari.matrix(mat.nrows(), mat.ncols(), mat.list())
        ker_pari_t = pari.matker(pari.mattranspose(mat_pari), 1)
        kerdim = len(ker_pari_t)
        if kerdim == 0:
            # empty matrix
            trans = matrix(ZZ, 0, mat.nrows())
        else:
            trans = matrix(ZZ, pari.mattranspose(ker_pari_t).Col().list())

        mat_pari = pari.matrix(mat.nrows(), mat.ncols(), mat.list())
        for i in range(kerdim, lllmat.nrows(), 1):
            lllmat_pari = pari.vector(lllmat.ncols(), lllmat[i].list())
            trans_pari_t = pari.matsolvemod(
                pari.mattranspose(mat_pari), 0, pari.mattranspose(lllmat_pari)
            )
            transele = matrix(ZZ, trans_pari_t.mattranspose().Col().list())
            trans = trans.stack(transele)
    else:
        trans = mat.kernel().matrix()
        kerdim = trans.nrows()

        for i in range(kerdim, lllmat.nrows(), 1):
            transele = mat.solve_left(lllmat[i])
            trans = trans.stack(transele)

    return trans


def _xgcd_list(intlst: List[int]) -> Tuple[int, List[int]]:
    """
    extended gcd algorithm for a_0,...,a_k
    input: [a_0, ..., a_k]
    output: d_, [b_0, ..., b_k] s.t. gcd(a_0,...,a_k) = d_, sum(a_i*b_i for i) = d_
    """

    if len(intlst) == 1:
        if intlst[0] >= 0:
            return intlst[0], [1]
        else:
            return -intlst[0], [-1]

    d, a, b = xgcd(intlst[0], intlst[1])

    curgcd = d
    curlst = [a, b]
    for i in range(2, len(intlst)):
        d, a, b = xgcd(curgcd, intlst[i])
        curlst = list(map(lambda x: x*a, curlst)) + [b]
        curgcd = d
    return curgcd, curlst


def do_LLL_fplll(mat: matrix, **kwds) -> Tuple[matrix, matrix]:
    if 'transformation' not in kwds:
        kwds['transformation'] = True
    if 'use_siegel' not in kwds:
        kwds['use_siegel'] = True
    if 'fplll_version' not in kwds:
        kwds['fplll_version'] = WRAPPER
    if 'early_reduction' not in kwds:
        kwds['early_reduction'] = True
    transformation = kwds['transformation']
    use_siegel = kwds['use_siegel']
    fplll_version = kwds['fplll_version']
    early_reduction = kwds['early_reduction']

    matstr = _from_sagematrix_to_fplllmatrix(mat)
    if early_reduction:
        result = subprocess_run(
            [os.path.join(fplll_path, 'fplll'), '-l', str(1-int(use_siegel)), '-m', fplll_version, '-y', '-of', 'u'],
            input=matstr.encode(), cwd=fplll_path, capture_output=True
        )
    else:
        result = subprocess_run(
            [os.path.join(fplll_path, 'fplll'), '-l', str(1-int(use_siegel)), '-m', fplll_version, '-of', 'u'],
            input=matstr.encode(), cwd=fplll_path, capture_output=True
        )
    if result.returncode != 0:
        print(result.stderr)
        raise ValueError(f"LLL failed with return code {result.returncode}")

    trans = _fplllmatrix_to_sagematrix(result.stdout.decode().strip())
    lllmat = trans * mat

    if not(transformation):
        trans = None

    return lllmat, trans


def do_BKZ_fplll(mat: matrix, **kwds) -> Tuple[matrix, matrix]:
    if 'transformation' not in kwds:
        kwds['transformation'] = True
    if 'blocksize' not in kwds:
        kwds['blocksize'] = 10
    if 'bkzautoabort' not in kwds:
        kwds['bkzautoabort'] = True
    transformation = kwds['transformation']
    blocksize = kwds['blocksize']
    bkzautoabort = kwds['bkzautoabort']

    matstr = _from_sagematrix_to_fplllmatrix(mat)
    if bkzautoabort:
        result = subprocess_run(
            [os.path.join(fplll_path, 'fplll'), '-a', 'bkz', '-b', str(blocksize), '-bkzautoabort', '-of', 'u'],
            input=matstr.encode(), cwd=fplll_path, capture_output=True
        )
    else:
        result = subprocess_run(
            [os.path.join(fplll_path, 'fplll'), '-a', 'bkz', '-b', str(blocksize), '-of', 'u'],
            input=matstr.encode(), cwd=fplll_path, capture_output=True
        )
    if result.returncode != 0:
        print(result.stderr)
        raise ValueError(f"LLL failed with return code {result.returncode}")

    trans = _fplllmatrix_to_sagematrix(result.stdout.decode().strip())
    lllmat = trans * mat

    if not(transformation):
        trans = None

    return lllmat, trans


def do_LLL_flatter(mat: matrix, **kwds) -> Tuple[matrix, matrix]:
    if 'transformation' not in kwds:
        kwds['transformation'] = True
    if 'use_pari_kernel' not in kwds:
        kwds['use_pari_kernel'] = True
    if 'use_pari_matsol' not in kwds:
        kwds['use_pari_matsol'] = False
    transformation = kwds['transformation']
    use_pari_kernel = kwds['use_pari_kernel']
    use_pari_matsol = kwds['use_pari_matsol']

    kerproc_st = time.time()

    if mat == zero_matrix(ZZ, mat.nrows(), mat.ncols()):
        return mat, identity_matrix(ZZ, mat.nrows())

    # sage has integer_kernel(), but somehow slow. instead using pari.matker
    if use_pari_kernel:
        mat_pari = pari.matrix(mat.nrows(), mat.ncols(), mat.list())
        ker_pari_t = pari.matker(mat_pari.mattranspose(), 1)
        ker = matrix(ZZ, ker_pari_t.mattranspose().Col().list())
    else:
        ker = mat.kernel().matrix()

    kerdim = ker.nrows()
    matrow = mat.nrows()
    col = mat.ncols()
    if kerdim == matrow: # full kernel
        return zero_matrix(ZZ, matrow, col), ker
    if kerdim == 0:
        Hsub = mat
        U = identity_matrix(ZZ, matrow)
    else:
        # heuristic construction for unimodular matrix which maps zero vectors on kernel
        # searching unimodular matrix can be done by HNF
        # (echeron_form(algorithm='pari') calls mathnf()),
        # but it is slow and produces big elements
        #
        # instead, searching determinant of submatrix = 1/-1,
        # then the determinant of whole unimodular matrix is det(submatrix)*(-1)^j
        # assume kernel has good property for gcd (gcd of some row elements might be 1)
        found_choice = False
        ker_submat_rows = tuple(range(kerdim))
        ker_submat_cols = []
        pivot = matrow - 1
        # search submatrix of kernel assuming last column vectors are triangulate
        while len(ker_submat_cols) < kerdim:
            if ker[ker_submat_rows, tuple([pivot])] != zero_matrix(ZZ, kerdim, 1):
                ker_submat_cols.append(pivot)
            pivot -= 1
        ker_submat_cols = tuple(sorted(ker_submat_cols))
        ker_last_det = int(ker[ker_submat_rows, ker_submat_cols].determinant())
        if ker_last_det == 0:
            raise ValueError("no unimodular matrix found (cause ker_last_det=0)")
        for choice in range(pivot, -1, -1):
            # gcd check
            gcd_row = ker_last_det
            for i in range(kerdim):
                gcd_row = GCD(gcd_row, ker[i, choice])
            if abs(gcd_row) != 1:
                continue

            # choice pivot: last columes for kernel are triangulated and small
            kersubidxes = [choice] + list(ker_submat_cols)
            detlst = [ker_last_det]
            for i in range(1, kerdim+1, 1):
                ker_submat_rows = tuple(range(kerdim))
                ker_submat_cols = tuple(kersubidxes[:i] + kersubidxes[i+1:])
                detlst.append(ker[ker_submat_rows, ker_submat_cols].determinant())
                detlist_gcd, detlist_coef = _xgcd_list(detlst)
                if detlist_gcd == 1:
                    found_choice = True
                    break
            if not found_choice:
                continue
            detlist_coef = detlist_coef + [0] * ((kerdim + 1) - len(detlist_coef))
            break
        if not found_choice:
            raise ValueError("no unimodular matrix found")
        U_top_vec = [0 for _ in range(matrow)]
        for i in range(kerdim+1):
            U_top_vec[kersubidxes[i]] = (-1)**i * detlist_coef[i]
        U_sub = matrix(ZZ, 1, matrow, U_top_vec)
        not_kersubidxes = sorted(list(set(list(range(matrow))) - set(kersubidxes)))
        for j in range(kerdim+1, matrow):
            onevec = [0 for _ in range(matrow)]
            onevec[not_kersubidxes[j-(kerdim+1)]] = 1
            U_sub = U_sub.stack(vector(ZZ, matrow, onevec))
        Hsub = U_sub * mat
        U = ker.stack(U_sub)
        #assert abs(U.determinant()) == 1
    kerproc_ed = time.time()
    logger.info("processing kernel elapsed time: %f", kerproc_ed - kerproc_st)

    if Hsub.nrows() == 1:
        lllmat = Hsub
    else:
        matstr = _from_sagematrix_to_fplllmatrix(Hsub)
        result = subprocess_run(
            os.path.join(flatter_path, 'flatter'),
            input=matstr.encode(), cwd=flatter_path, capture_output=True
        )
        if result.returncode != 0:
            print(result.stderr)
            raise ValueError(f"LLL failed with return code {result.returncode}")
        lllmat = _fplllmatrix_to_sagematrix(result.stdout.decode().strip())

    if transformation:
        trans = _transformation_matrix(Hsub, lllmat, use_pari_matsol=use_pari_matsol)
    else:
        trans = None

    restrows = mat.nrows() - lllmat.nrows()
    final_lllmat = zero_matrix(ZZ, restrows, lllmat.ncols()).stack(lllmat)

    if transformation:
        middle_trans = identity_matrix(ZZ, restrows).augment(zero_matrix(ZZ, restrows, trans.ncols())).stack(
            zero_matrix(ZZ, trans.nrows(), restrows).augment(trans)
        )
        final_trans = middle_trans * U
        #assert abs(final_trans.determinant()) == 1
        #assert final_trans * mat == final_lllmat
    else:
        final_trans = None

    return final_lllmat, final_trans


def do_LLL_NTL(mat: matrix, **kwds) -> Tuple[matrix, matrix]:
    if 'transformation' not in kwds:
        kwds['transformation'] = True
    transformation = kwds['transformation']

    delta_lll = ZZ(99)/ZZ(100)
    a_lll = delta_lll.numer()
    b_lll = delta_lll.denom()

    A = ntl_mat(mat)

    # TODO: support various floating point precision, and use_givens option
    r, det2, U = A.LLL(a_lll, b_lll, return_U=transformation)

    lllmat = matrix(ZZ, mat.nrows(), mat.ncols(), [ZZ(z) for z in A.list()])
    if transformation:
        trans = matrix(ZZ, mat.nrows(), mat.nrows(), [ZZ(z) for z in U.list()])
    else:
        trans = None

    return lllmat, trans


def do_BKZ_NTL(mat: matrix, **kwds) -> Tuple[matrix, matrix]:
    if 'transformation' not in kwds:
        kwds['transformation'] = True
    if 'blocksize' not in kwds:
        kwds['blocksize'] = 10
    if 'prune' not in kwds:
        kwds['prune'] = 0
    transformation = kwds['transformation']
    blocksize = kwds['blocksize']
    prune = kwds['prune']

    delta_lll = 0.99
    A = ntl_mat(mat)
    U = ntl_mat(identity_matrix(ZZ, A.nrows()))

    # TODO: support various floating point precision, and use_givens option
    r = A.BKZ_RR(U=U, delta=delta_lll, BlockSize=blocksize, prune=prune)

    lllmat = matrix(ZZ, mat.nrows(), mat.ncols(), [ZZ(z) for z in A.list()])
    if transformation:
        trans = matrix(ZZ, mat.nrows(), mat.nrows(), [ZZ(z) for z in U.list()])
    else:
        trans = None

    return lllmat, trans


## wrapper function
def do_lattice_reduction(mat: matrix, **kwds) -> Tuple[matrix, matrix]:
    """
    LLL/BKZ reduction
    input: (mat, algorithm, **kwds)
            - mat: target lattice representation matrix for LLL/BKZ reduction
            - algorithm: int value which specify which algorithm will be used
              (FPLLL, FPLLL_BKZ, FLATTER, NTL, NTL_BKZ)
    output: (lllmat, trans)
            - lllmat: LLL/BKZ reduced basis matrix (might include zero-vectors)
            - trans: transformation matrix s.t. lllmat = trans * mat
    """
    if 'algorithm' not in kwds:
        kwds['algorithm'] = FLATTER
    algorithm = kwds['algorithm']

    logger.info("size of mat for lattice reduction: (%d, %d)", int(mat.nrows()), int(mat.ncols()))
    logger.debug(
        "lattice reduction param: algorithm=%s, param=%s",
        LLL_algorithm_str[algorithm], str(kwds)
    )
    logger.info("start lattice reduction")
    st = time.time()

    result = LLL_algorithm_dict[algorithm](mat, **kwds)

    ed = time.time()
    logger.info("end lattice reduction. elapsed %f", ed-st)

    return result


def babai(mat: matrix, target: vector, algorithm: int = FLATTER, **kwds) -> Tuple[vector, vector]:
    """
    Babai nearlest plain algorithm for solving CVP
    input: (mat, target, **kwds)
            - mat: lattice representation matrix for LLL/BKZ reduction
            - target: target integer vector for solving close point in lattice
            - algorithm: int value which specify which algorithm will be used
              (FPLLL, FPLLL_BKZ, FLATTER, NTL, NTL_BKZ)
    output: (diff, trans)
            - diff: subtract of target from lattice_point which is close for target
            - trans: transformation matrix s.t. lattice_point = trans * mat
    """
    kwds['transformation'] = True
    lll, trans = do_lattice_reduction(mat, algorithm, **kwds)
    # gram-schmidt process is slow. use solve_left in QQ
    sol_QQ = (lll.change_ring(QQ)).solve_left((target.change_ring(QQ)))
    sol_approx_ZZ_lst = [ZZ(QQ(sol_QQ_ele).round()) for sol_QQ_ele in sol_QQ.list()]
    sol_approx_ZZ = vector(ZZ, len(sol_approx_ZZ_lst), sol_approx_ZZ_lst)
    return target - sol_approx_ZZ * lll, sol_approx_ZZ * trans


def enumeration(mat: matrix, bound: int, target: vector = None, algorithm: int = FLATTER, **kwds):
    """
    Enumeration (SVP or CVP)
    input: (mat, target, **kwds)
            - mat: lattice representation matrix for LLL/BKZ reduction
            - target: None for SVP, target integer vector for solving close point in lattice for CVP
              bound: expected norm size estimation as bound = sqrt((L2-norm**2) / size)
            - algorithm: int value which specify which algorithm will be used
              (FPLLL, FPLLL_BKZ, FLATTER, NTL, NTL_BKZ)
    output: enumeration generator
    """
    kwds['transformation'] = True
    lll, trans = do_lattice_reduction(mat, algorithm, **kwds)
    lllele = []
    for i in range(0, lll.nrows()):
        lllele.append(lll[i].list())
    lll_fpylll = IntegerMatrix.from_matrix(lllele)
    MG = GSO.Mat(lll_fpylll)
    MG.update_gso()
    enum = Enumeration(MG)
    size = lll.ncols()
    answers = enum.enumerate(0, size, size * (bound ** 2), 0, target=target, pruning=None)
    for _, s in answers:
        v = vector(ZZ, size, list(map(int, s)))
        enumresult = v * trans
        yield (enumresult * mat, enumresult)


def test():
    testlst = [
        ("zerodim", [[0,0,0]]),
        ("onedim", [[1,2,3]]),
        ("twodim_indep", [[1,2,3],[4,5,6]]),
        ("twodim_dep", [[1,2,3],[2,4,6]]),
        ("threedim_indep", [[1,2,3],[4,5,6],[7,8,9]]),
        ("threedim_one_dep", [[1,2,3],[2,4,6],[8,9,10]]),
        ("threedim_two_dep", [[1,2,3],[2,4,6],[3,6,9]]),
        ("overdim", [[1,2,3],[4,5,6],[7,8,9],[10,11,12]]),
        ("overdim_onedep", [[1,2,3],[4,5,6],[3,6,9],[5,6,7]]),
        ("multiple_2_ker", [[-2,-4,-6],[1,2,3],[3,6,9]]),
    ]

    for LLL_algorithm in range(5):
        print(f"LLL_algorithm: {LLL_algorithm_str[LLL_algorithm]}")
        for testlstele in testlst:
            curmat = matrix(ZZ, testlstele[1])
            try:
                if LLL_algorithm == FLATTER:
                    lll, trans = LLL_algorithm_dict[LLL_algorithm](curmat, **{'use_pari_kernel':False})
                    #lll, trans = LLL_algorithm_dict[LLL_algorithm](curmat, use_pari_kernel=True)
                else:
                    lll, trans = LLL_algorithm_dict[LLL_algorithm](curmat, **{})
            except:
                traceback.print_exc()
                continue

            print(f"test {testlstele[1]}: {(trans * curmat == lll, abs(trans.determinant()) == 1)}")
            print((lll.rows(), trans.rows()))
            print("")


LLL_algorithm_dict = {
    FLATTER: do_LLL_flatter,
    FPLLL: do_LLL_fplll, FPLLL_BKZ: do_BKZ_fplll,
    NTL: do_LLL_NTL, NTL_BKZ: do_BKZ_NTL
}

LLL_algorithm_str = {
    FLATTER: 'FLATTER',
    FPLLL: 'FPLLL', FPLLL_BKZ: 'FPLLL_BKZ',
    NTL: 'NTL', NTL_BKZ: 'NTL_BKZ'
}


if __name__ == '__main__':
    test()
