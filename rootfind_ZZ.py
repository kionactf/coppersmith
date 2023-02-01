from sage.all import *

from random import shuffle as random_shuffle
from itertools import product as itertools_product
import time

from logger import logger


def solve_root_onevariable(pollst, bounds):
    logger.info("start solve_root_onevariable")
    st = time.time()

    for f in pollst:
        f_x = f.parent().gens()[0]
        try:
            rt_ = f.change_ring(ZZ).roots()
            rt = [ele for ele, exp in rt_]
        except:
            f_QQ = f.change_ring(QQ)
            f_QQ_x = f_QQ.parent().gens()[0]
            rt_ = f_QQ.parent().ideal([f_QQ]).variety()
            rt = [ele[f_QQ_x] for ele in rt_]
        if rt != []:
            break
    result = []
    for rtele in rt:
        if any([pollst[i].subs({f_x: int(rtele)}) != 0 for i in range(len(pollst))]):
            continue
        if abs(int(rtele)) < bounds[0]:
            result.append(rtele)

    ed = time.time()
    logger.info("end solve_root_onevariable. elapsed %f", ed-st)

    return result


def solve_root_groebner(pollst, bounds):
    logger.info("start solve_root_groebner")
    st = time.time()

    # I heard degrevlex is faster computation for groebner basis, but idk real effect
    polrng_QQ = pollst[0].change_ring(QQ).parent().change_ring(order='degrevlex')
    vars_QQ = polrng_QQ.gens()
    G = Sequence(pollst, polrng_QQ).groebner_basis()
    try:
        # not zero-dimensional ideal raises error
        rt_ = G.ideal().variety()
    except:
        logger.warning("variety failed. not zero-dimensional ideal?")
        return None
    rt = [[int(ele[v]) for v in vars_QQ] for ele in rt_]

    vars_ZZ = pollst[0].parent().gens()
    result = []
    for rtele in rt:
        if any([pollst[i].subs({v: int(rtele[i]) for i, v in enumerate(vars_ZZ)}) != 0 for i in range(len(pollst))]):
            continue
        if all([abs(int(rtele[i])) < bounds[i] for i in range(len(rtele))]):
            result.append(rtele)
    
    ed = time.time()
    logger.info("end solve_root_groebner. elapsed %f", ed-st)
    return result


def solve_ZZ_symbolic_linear_internal(sol_coefs, bounds):
    mult = prod(bounds)
    matele = []
    for i, sol_coef in enumerate(sol_coefs):
        denom = 1
        for sol_coef_ele in sol_coef:
            denom = LCM(denom, sol_coef_ele.denominator())
        for sol_coef_ele in sol_coef:
            matele.append(ZZ(sol_coef_ele * denom * mult))
        matele += [0]*i + [-mult*denom] + [0] * (len(bounds)-i-1)
    for idx, bd in enumerate(bounds):
        matele += [0]*len(sol_coefs[0]) + [0] * idx + [mult//bd] + [0]*(len(bounds)-idx-1)
    # const term
    matele += [0]*(len(sol_coefs[0])-1) + [mult] + [0]*len(bounds)
    mat = matrix(ZZ, len(sol_coefs)+len(bounds)+1, len(sol_coefs[0])+len(bounds), matele)
    logger.debug(f"start LLL for solve_ZZ_symbolic_linear_internal")
    mattrans = mat.transpose()
    lll, trans = mattrans.LLL(transformation=True)
    logger.debug(f"end LLL")
    for i in range(trans.nrows()):
        if all([lll[i, j] == 0 for j in range(len(sol_coefs))]):
            if int(trans[i,len(sol_coefs[0])-1]) in [1,-1]:
                linsolcoef = [int(trans[i,j])*int(trans[i,len(sol_coefs[0])-1]) for j in range(len(sol_coefs[0]))]
                logger.debug(f"linsolcoef found: {linsolcoef}")
                linsol = []
                for sol_coef in sol_coefs:
                    linsol.append(sum([ele*linsolcoef[idx] for idx, ele in enumerate(sol_coef)]))
                return [linsol]
    return []


def solve_root_triangulate(pollst, bounds):
    logger.info("start solve_root_triangulate")
    st = time.time()

    polrng_QQ = pollst[0].change_ring(QQ).parent().change_ring(order='lex')
    vars_QQ = polrng_QQ.gens()
    G = Sequence(pollst, polrng_QQ).groebner_basis()
    if len(G) == 0:
        return []

    symbolic_vars = [var(G_var) for G_var in G[0].parent().gens()]
    try:
        sols = solve([G_ele(*symbolic_vars) for G_ele in G], symbolic_vars, solution_dict=True)
    except:
        return None

    logger.debug(f"found sol on triangulate: {sols}")

    result = []
    # solve method returns parametrized solution. We treat only linear equation
    # TODO: use solver for more general integer equations (such as diophautus solver, integer programming solver, etc.)
    for sol in sols:
        sol_args = set()
        for symbolic_var in symbolic_vars:
            sol_var = sol[symbolic_var]
            sol_args = sol_args.union(set(sol_var.args()))

        sol_args = list(sol_args)
        sol_coefs = []
        for symbolic_var in symbolic_vars:
            sol_var = sol[symbolic_var]
            sol_coefs_ele = []
            for sol_arg in sol_args:
                if sol_var.is_polynomial(sol_arg) == False:
                    logger.warning("cannot deal with non-polynomial equation")
                    return None
                if sol_var.degree(sol_arg) > 1:
                    logger.warning("cannot deal with high degree equation")
                    return None
                sol_var_coef_arg = sol_var.coefficient(sol_arg)
                if sol_var_coef_arg not in QQ:
                    logger.warning("cannot deal with multivariate non-linear equation")
                    return None
                sol_coefs_ele.append(QQ(sol_var_coef_arg))
            # constant term
            const = sol_var.subs({sol_arg: 0 for sol_arg in sol_args})
            if const not in QQ:
                return None
            sol_coefs_ele.append(const)

            sol_coefs.append(sol_coefs_ele)
        ZZsol = solve_ZZ_symbolic_linear_internal(sol_coefs, bounds)
        result += ZZsol

    ed = time.time()
    logger.info("end solve_root_triangulate. elapsed %f", ed-st)
    return result


def solve_root_jacobian_newton_internal(pollst, startpnt):
    # NOTE: Newton method's complexity is larger than BFGS, but for small variables Newton method converges soon.
    pollst_Q = Sequence(pollst, pollst[0].parent().change_ring(QQ))
    vars_pol = pollst_Q[0].parent().gens()
    jac = jacobian(pollst_Q, vars_pol)

    if all([ele == 0 for ele in startpnt]):
        # just for prepnt != pnt
        prepnt = {vars_pol[i]: 1 for i in range(len(vars_pol))}
    else:
        prepnt = {vars_pol[i]: 0 for i in range(len(vars_pol))}
    pnt = {vars_pol[i]: startpnt[i] for i in range(len(vars_pol))}

    maxiternum = 1024
    iternum = 0
    while True:
        if iternum >= maxiternum:
            logger.warning("failed. maybe, going wrong way.")
            return None

        evalpollst = [(pollst_Q[i].subs(pnt)) for i in range(len(pollst_Q))]
        if all([int(ele) == 0 for ele in evalpollst]):
            break
        jac_eval = jac.subs(pnt)
        evalpolvec = vector(QQ, len(evalpollst), evalpollst)
        try:
            pnt_diff_vec = jac_eval.solve_right(evalpolvec)
        except:
            logger.warning("pnt_diff comp failed.")
            return None

        prepnt = {key:value for key,value in prepnt.items()}
        pnt = {vars_pol[i]: round(QQ(pnt[vars_pol[i]] - pnt_diff_vec[i])) for i in range(len(pollst_Q))}

        if all([prepnt[vars_pol[i]] == pnt[vars_pol[i]] for i in range(len(vars_pol))]):
            logger.warning("point update failed. (converged local sol)")
            return None
        prepnt = {key:value for key,value in pnt.items()}
        iternum += 1
    return [int(pnt[vars_pol[i]]) for i in range(len(vars_pol))]


def solve_root_jacobian_newton(pollst, bounds):
    logger.info("start solve_root_jacobian newton")
    st = time.time()

    pollst_local = pollst[:]
    vars_pol = pollst[0].parent().gens()

    # not applicable to non-determined system
    if len(vars_pol) > len(pollst):
        return []

    for _ in range(10):
        random_shuffle(pollst_local)
        for signs in itertools_product([1, -1], repeat=len(vars_pol)):
            startpnt = [signs[i] * bounds[i] for i in range(len(vars_pol))]
            result = solve_root_jacobian_newton_internal(pollst_local[:len(vars_pol)], startpnt)
            # filter too much small solution
            if result is not None:
                if all([abs(ele) < 2**16 for ele in result]):
                    continue
                ed = time.time()
                logger.info("end solve_root_jacobian newton. elapsed %f", ed-st)
                return [result]


def _solve_root_GF_smallp(pollst, smallp):
    Fsmallp = GF(smallp)
    polrng_Fsmallp = pollst[0].change_ring(Fsmallp).parent().change_ring(order='degrevlex')
    vars_Fsmallp = polrng_Fsmallp.gens()
    fieldpolys = [varele**smallp - varele for varele in vars_Fsmallp]
    pollst_Fsmallp = [polrng_Fsmallp(ele) for ele in pollst]
    G = pollst_Fsmallp[0].parent().ideal(pollst_Fsmallp + fieldpolys).groebner_basis()
    rt_ = G.ideal().variety()
    rt = [[int(ele[v].lift()) for v in vars_Fsmallp] for ele in rt_]
    return rt


def solve_root_hensel_smallp(pollst, bounds, smallp):
    logger.info("start solve_root_hensel")
    st = time.time()

    vars_ZZ = pollst[0].parent().gens()
    smallp_exp_max = max([int(log(ele, smallp)+0.5) for ele in bounds]) + 1
    # firstly, compute low order
    rt_lows = _solve_root_GF_smallp(pollst, smallp)
    for smallp_exp in range(1, smallp_exp_max+1, 1):
        cur_rt_low = []
        for rt_low in rt_lows:
            evalpnt = {vars_ZZ[i]:(smallp**smallp_exp)*vars_ZZ[i]+rt_low[i] for i in range(len(vars_ZZ))}
            nextpollst = [pol.subs(evalpnt)/(smallp**smallp_exp) for pol in pollst]
            rt_up = _solve_root_GF_smallp(nextpollst, smallp)
            cur_rt_low += [tuple([smallp**smallp_exp*rt_upele[i] + rt_low[i] for i in range(len(rt_low))]) for rt_upele in rt_up]
        rt_lows = list(set(cur_rt_low))
        if len(rt_lows) >= 800:
            logger.warning("too much root candidates found")
            return None

    result = []
    for rt in rt_lows:
        rtele = [[ele, ele - smallp**(smallp_exp_max+1)][ele >= smallp**smallp_exp_max] for ele in rt]
        if any([pollst[i].subs({v: int(rtele[i]) for i, v in enumerate(vars_ZZ)}) != 0 for i in range(len(pollst))]):
            continue
        if all([abs(int(rtele[i])) < bounds[i] for i in range(len(rtele))]):
            result.append(rtele)

    ed = time.time()
    logger.info("end solve_root_hensel. elapsed %f", ed-st)
    return result


def solve_root_hensel(pollst, bounds):
    for smallp in [2, 3, 5]:
        result = solve_root_hensel_smallp(pollst, bounds, smallp)
        if result != [] and result is not None:
            return result
    return None


## wrapper function
def rootfind_ZZ(pollst, bounds):
    vars_pol = pollst[0].parent().gens()
    if len(vars_pol) != len(bounds):
        raise ValueError("vars len is invalid (on rootfind_ZZ)")

    # Note: match-case statement introduced on python3.10, but not used for backward compati
    if len(vars_pol) == 1:
        return solve_root_onevariable(pollst, bounds)
    else:
        # first numeric
        result = solve_root_jacobian_newton(pollst, bounds)
        if result != [] and result is not None:
            return result

        # next hensel (fast if the number of solutions mod smallp**a are small. in not case, cannot find solution)
        result = solve_root_hensel(pollst, bounds)
        if result != [] and result is not None:
            return result

        # last triangulate with groebner (slow, but sometimes solve when above methods does not work)
        #return solve_root_groebner(pollst, bounds)
        return solve_root_triangulate(pollst, bounds)

