from sage.all import *

from typing import Tuple, List
from subprocess import run as subprocess_run
from re import sub as re_sub
import time

from logger import logger


fplll_path = './fplll/fplll' # /usr/bin
flatter_path = './flatter/build/bin' # /usr/bin

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


def _transformation_matrix(mat, lllmat):
    trans = mat.solve_left(lllmat)
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
        use_siegel : int = True, fplll_version : str = WRAPPER, early_reduction : bool = True
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
    trans = _transformation_matrix(mat, lllmat)
    return lllmat, trans


def do_BKZ_fplll(
        mat: matrix,
        blocksize : int = 10, bkzautoabort : bool = True
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
    trans = _transformation_matrix(mat, lllmat)
    return lllmat, trans


def do_LLL_flatter(
        mat: matrix
    ) -> Tuple[matrix, matrix]:

    kerproc_st = time.time()
    
    # sage has integer_kernel(), but somehow slow. instead using pari.matker
    mat_pari = pari.matrix(mat.nrows(), mat.ncols(), mat.list())
    ker_pari_t = pari.matker(mat_pari.mattranspose(), 1)
    ker = matrix(ZZ, ker_pari_t.mattranspose().Col().list())
    kerdim = ker.nrows()
    matrow = mat.nrows()
    col = mat.ncols()
    if kerdim == matrow: # full kernel
        return zero_matrix(ZZ, matrow, col), ker
    if kerdim == 0:
        Hsub = mat
        U = identity_matrix(ZZ, matrow)
    else: # heuristic construct unimodular matrix for mapping zero vectors on first kernel dimension
        # searching unimodular matrix can be done by HNF (echoron_form(algorithm='pari') calls mathnf(), but it is slow and produces big elements)
        # instead, seaching determinant of submatrix = 1/-1, then the determinant of whole unimodular matrix is det(submatrix)*(-1)^j
        # assume kernel has good property for gcd (gcd of some row elements might be 1)
        found_choice = False
        ker_submat_rows = tuple(range(kerdim))
        ker_submat_cols = tuple(range(matrow-kerdim, matrow, 1))
        ker_last_det = int(ker[ker_submat_rows, ker_submat_cols].determinant())
        if ker_last_det == 0:
            raise ValueError("no unimodular matrix found (cause ker_last_det=0)")
        for choice in range(matrow-kerdim-1, -1, -1):
            # gcd check
            gcd_row = ker_last_det
            for i in range(kerdim):
                gcd_row = GCD(gcd_row, ker[i, choice])
            if abs(gcd_row) != 1:
                continue

            # choice pivot: last columes for kernel are triangulated and small
            kersubidxes = [choice] + list(range(matrow-kerdim, matrow, 1))
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

    matstr = _from_sagematrix_to_fplllmatrix(Hsub)
    result = subprocess_run(
        './flatter',
        input=matstr.encode(), cwd=flatter_path, capture_output=True
    )
    if result.returncode != 0:
        print(result.stderr)
        raise ValueError(f"LLL failed with return code {result.returncode}")
    lllmat = _fplllmatrix_to_sagematrix(result.stdout.decode().strip())
    trans = _transformation_matrix(Hsub, lllmat)

    restrows = mat.nrows() - lllmat.nrows()
    final_lllmat = zero_matrix(ZZ, restrows, lllmat.ncols()).stack(lllmat)
    middle_trans = identity_matrix(ZZ, restrows).augment(zero_matrix(ZZ, restrows, trans.ncols())).stack(
        zero_matrix(ZZ, trans.nrows(), restrows).augment(trans)
    )
    final_trans = middle_trans * U
    #assert abs(final_trans.determinant()) == 1
    #assert final_trans * mat == final_lllmat

    return final_lllmat, final_trans


def do_LLL_NTL(mat: matrix) -> Tuple[matrix, matrix]:
    return mat.LLL(algorithm="NTL:LLL", transformation=True)


def do_BKZ_NTL(
        mat: matrix,
        blocksize : int = 10, prune : int = 0
    ) -> Tuple[matrix, matrix]:
    lllmat = mat.BKZ(algorithm="NTL", block_size=blocksize, prune=prune)
    trans = _transformation_matrix(mat, lllmat)
    return lllmat, trans


## wrapper function
def do_lattice_reduction(mat: matrix, algorithm: int = FLATTER, **kwds) -> Tuple[matrix, matrix]:
    algorithm_dict = {FLATTER: do_LLL_flatter, FPLLL: do_LLL_fplll, FPLLL_BKZ: do_BKZ_fplll, NTL: do_LLL_NTL, NTL_BKZ: do_BKZ_NTL}
    return algorithm_dict[algorithm](mat, **kwds)
