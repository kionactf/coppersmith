from rootfind_ZZ import JACOBIAN, HENSEL, TRIANGULATE, GROEBNER, LINEAR_SIMPLE, LINEAR_NEAR_BOUNDS
from lll import do_lattice_reduction, FPLLL, FPLLL_BKZ, FLATTER, NTL, NTL_BKZ, WRAPPER, HEURISTIC


class ContextClass():
    def __init__(self):
        self.__rootfindZZopt = {}
        self.__lllopt = {}

    ## TODO: custom check for lllopt
    @property
    def lllopt(self):
        return self.__lllopt
    
    @lllopt.setter
    def lllopt(self, newlllopt):
        self.__lllopt = newlllopt

    ## TODO: custom check for rootfindZZopt
    @property
    def rootfindZZopt(self):
        return self.__rootfindZZopt
    
    @rootfindZZopt.setter
    def rootfindZZopt(self, newrootfindZZopt):
        self.__rootfindZZopt = newrootfindZZopt


context = ContextClass()


# NOTE: not depend on contextclass for specific coppersmith libraries
# (cause these libraries itself would be reused for other objectives)
# instead, we add options here, and call functions in these libraries with context options


## lll
def register_options_lll(_context):
    lllopt = {}

    ## general
    lllopt['algorithm'] = FPLLL # for coppersmith(if need, use FLATTER or FPLLL_BKZ)
    lllopt['transformation'] = False # for coppersmith (normally, True)

    ## BKZ(FPLLL, NTL)
    lllopt['blocksize'] = 10

    ## FPLLL
    lllopt['use_siegel'] = True
    lllopt['fplll_version'] = WRAPPER
    lllopt['early_reduction'] = True

    ## FPLLL_BKZ
    lllopt['bkzautoabort'] = True

    ## FLATTER
    lllopt['use_pari_kernel'] = True
    lllopt['use_pari_matsol'] = False

    ## NTL
    # (none)

    ## NTL_BKZ
    lllopt['prune'] = 0


    _context.lllopt = lllopt


## rootfind_ZZ
def register_options_rootfind_ZZ(_context):
    rootfindZZopt = {}
    # NOTE: we use dict for each libraries options (instead of property)
    # (cause avoiding too complicating debug and leave readability)

    ## general
    rootfindZZopt['algorithms'] = (JACOBIAN, HENSEL, TRIANGULATE)
    rootfindZZopt['monomial_order_for_variety'] = 'degrevlex'

    ## JACOBIAN
    rootfindZZopt['maxiternum'] = 1024
    rootfindZZopt['select_subpollst_loopnum'] = 10
    rootfindZZopt['search_near_positive_bounds_only'] = False
    rootfindZZopt['filter_small_solution_minbound'] = 2**16

    ## HENSEL
    rootfindZZopt['smallps'] = (2, 3, 5)
    rootfindZZopt['maxcands'] = 800

    ## TRIANGULATE
    rootfindZZopt['lllopt_symbolic_linear'] = {'algorithm': FPLLL_BKZ}
    rootfindZZopt['symbolic_linear_algorithm'] = LINEAR_NEAR_BOUNDS


    _context.rootfindZZopt = rootfindZZopt


register_options_lll(context)
register_options_rootfind_ZZ(context)
