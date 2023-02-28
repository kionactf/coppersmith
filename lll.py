from sage.all import *

from typing import Tuple, List
from subprocess import run as subprocess_run
from re import sub as re_sub
import time
import traceback

from logger import logger


coppersmith_dir = '/home/sage/coppersmith/'

fplll_path = coppersmith_dir + 'fplll/fplll' # /usr/bin
flatter_path = coppersmith_dir + 'flatter/build/bin' # /usr/bin

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


def _transformation_matrix(mat, lllmat, use_pari_kernel=True):
    # pari.matker() does not assure smallest kernel in Z (seems not call hermite normal form)
    # Sage kernel calls hermite normal form
    #
    # for computing ZZ transformation, use pari.matker, pari.matsolvemod
    # assume first kerdim vectors for lllmat are zero vector
    if use_pari_kernel:
        mat_pari = pari.matrix(mat.nrows(), mat.ncols(), mat.list())
        ker_pari_t = pari.matker(pari.mattranspose(mat_pari), 1)
        kerdim = len(ker_pari_t)
        if kerdim == 0:
            # empty matrix
            trans = matrix(ZZ, 0, mat.nrows())
        else:
            trans = matrix(ZZ, pari.mattranspose(ker_pari_t).Col().list())
    else:
        trans = mat.kernel().matrix()
        kerdim = trans.nrows()

    mat_pari = pari.matrix(mat.nrows(), mat.ncols(), mat.list())
    for i in range(kerdim, lllmat.nrows(), 1):
        lllmat_pari = pari.vector(lllmat.ncols(), lllmat[i].list())
        trans_pari_t = pari.matsolvemod(
            pari.mattranspose(mat_pari), 0, pari.mattranspose(lllmat_pari)
        )
        transele = matrix(ZZ, trans_pari_t.mattranspose().Col().list())
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


def do_LLL_fplll(
        mat: matrix,
        use_siegel : int = True, fplll_version : str = WRAPPER, early_reduction : bool = True,
        use_pari_kernel : bool = True
    ) -> Tuple[matrix, matrix]:
    matstr = _from_sagematrix_to_fplllmatrix(mat)
    if early_reduction:
        result = subprocess_run(
            ['./fplll', '-l', str(1-int(use_siegel)), '-m', fplll_version, '-y'],
            input=matstr.encode(), cwd=fplll_path, capture_output=True
        )
    else:
        result = subprocess_run(
            ['./fplll', '-l', str(1-int(use_siegel)), '-m', fplll_version],
            input=matstr.encode(), cwd=fplll_path, capture_output=True
        )
    if result.returncode != 0:
        print(result.stderr)
        raise ValueError(f"LLL failed with return code {result.returncode}")
    lllmat = _fplllmatrix_to_sagematrix(result.stdout.decode().strip())
    trans = _transformation_matrix(mat, lllmat, use_pari_kernel=use_pari_kernel)
    return lllmat, trans


def do_BKZ_fplll(
        mat: matrix,
        blocksize : int = 10, bkzautoabort : bool = True, use_pari_kernel : bool = True
    ) -> Tuple[matrix, matrix]:
    matstr = _from_sagematrix_to_fplllmatrix(mat)
    if bkzautoabort:
        result = subprocess_run(
            ['./fplll', '-a', 'bkz', '-b', str(blocksize), '-bkzautoabort'],
            input=matstr.encode(), cwd=fplll_path, capture_output=True
        )
    else:
        result = subprocess_run(
            ['./fplll', '-a', 'bkz', '-b', str(blocksize)],
            input=matstr.encode(), cwd=fplll_path, capture_output=True
        )
    if result.returncode != 0:
        print(result.stderr)
        raise ValueError(f"LLL failed with return code {result.returncode}")
    lllmat = _fplllmatrix_to_sagematrix(result.stdout.decode().strip())
    trans = _transformation_matrix(mat, lllmat, use_pari_kernel=use_pari_kernel)
    return lllmat, trans


def do_LLL_flatter(
        mat: matrix, use_pari_kernel: bool = True
    ) -> Tuple[matrix, matrix]:

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
            './flatter',
            input=matstr.encode(), cwd=flatter_path, capture_output=True
        )
        if result.returncode != 0:
            print(result.stderr)
            raise ValueError(f"LLL failed with return code {result.returncode}")
        lllmat = _fplllmatrix_to_sagematrix(result.stdout.decode().strip())

    trans = _transformation_matrix(Hsub, lllmat, use_pari_kernel=use_pari_kernel)

    restrows = mat.nrows() - lllmat.nrows()
    final_lllmat = zero_matrix(ZZ, restrows, lllmat.ncols()).stack(lllmat)
    middle_trans = identity_matrix(ZZ, restrows).augment(zero_matrix(ZZ, restrows, trans.ncols())).stack(
        zero_matrix(ZZ, trans.nrows(), restrows).augment(trans)
    )
    final_trans = middle_trans * U
    #assert abs(final_trans.determinant()) == 1
    #assert final_trans * mat == final_lllmat

    return final_lllmat, final_trans


def do_LLL_NTL(mat: matrix, use_pari_kernel: bool = True) -> Tuple[matrix, matrix]:
    lllmat = mat.LLL(algorithm="NTL:LLL")
    trans = _transformation_matrix(mat, lllmat, use_pari_kernel=use_pari_kernel)
    return lllmat, trans


def do_BKZ_NTL(
        mat: matrix,
        blocksize : int = 10, prune : int = 0, use_pari_kernel: bool = True
    ) -> Tuple[matrix, matrix]:
    lllmat = mat.BKZ(algorithm="NTL", block_size=blocksize, prune=prune)
    trans = _transformation_matrix(mat, lllmat, use_pari_kernel=use_pari_kernel)
    return lllmat, trans


## wrapper function
def do_lattice_reduction(mat: matrix, algorithm: int = FLATTER, **kwds) -> Tuple[matrix, matrix]:
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
    lll, trans = do_lattice_reduction(mat, algorithm, **kwds)
    # gram-schmidt process is slow. use solve_left in QQ
    sol_QQ = (lll.change_ring(QQ)).solve_left((target.change_ring(QQ)))
    sol_approx_ZZ_lst = [ZZ(QQ(sol_QQ_ele).round()) for sol_QQ_ele in sol_QQ.list()]
    sol_approx_ZZ = vector(ZZ, len(sol_approx_ZZ_lst), sol_approx_ZZ_lst)
    return target - sol_approx_ZZ * lll, sol_approx_ZZ * trans


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
                lll, trans = LLL_algorithm_dict[LLL_algorithm](curmat, use_pari_kernel=False)
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
