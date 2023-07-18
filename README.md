# Coppersmith small roots

Coppersmith small roots method (solving polynomial equation over composite modulus on small bounds)

## Features

- multiple coppersmith equation support
  - univariate polynomial
  - multivariate linear polynomial (Herrmann-May method)
  - multivariate polynomial (Jochemsz-May heuristic method)

Mainly, dealing with theoretically established Coppersmith method applicable equation. We recommend using univariate or linear (instead of heuristic) if you know the type of equation.

- automated parameter tuning

Firstly choose suitable parameters based on bounds and $\beta$, and then increment parameters until finishing to find an integer equation.

- selectable internal LLL/BKZ method
  - [fplll](https://github.com/fplll/fplll) (LLL/BKZ)
  - [flatter](https://github.com/keeganryan/flatter/)
  - NTL (via Sagemath LLL/BKZ)

- selectable solving integer equation method
  - jacobian Newton method (numerical)
  - Hensel lifting
  - triangulation (Groebner basis + dealing with non-zero dimensional ideal partially)

## Installation

```bash
$ git clone https://github.com/kionactf/coppersmith --sync
```

If you use fplll or flatter, need to build these packages. if you already installed them, you can use this by setting these install paths to fplll_path or flatter_path on lll.py.

Or you can create a docker image. It is based on sagemath/sagemath:latest image and also downloads well-known lattice libraries. (If you know good other libraries, let us know for including these.)

```bash
$ docker build -t coppersmith .
$ docker run --rm -it coppersmith /bin/bash
```

## Usage

Call coppersmith_onevariable.coppersmith_onevariable or coppersmith_linear.coppersmith_linear with Sagemath PolynomialRing over Zmod(N), bounds, beta.

See `example.py`.

Also, you can use only `LLL `by calling lll.do_lattice_reduction with Sagemath matrix over ZZ and some optimization parameters. And lll.babai for CVP solver.

For `integer equations solver`, use rootfind_ZZ.rootfind_ZZ with a list of Sagemath polynomial over ZZ and bounds.

## Note (use_pari_kernel)

For computing LLL, we use pari.matker for eliminating linearly dependent vectors for defining lattice. The process needs to use flatter. Though pari.matker is fast and outputs correct results in many cases, it sometimes outputs `wrong` results. (You can check this by running lll.test().) You may disable to use pari.matker by setting the option `use_pari_kernel=False`, where it forces using Sagemath kernel (which do internally run computing hermite normal form (HNF).) Note that HNF tends to produce large elements, so it causes LLL process makes super slow.

## Background

See [Why we could not solve chronophobia… Analysis for Coppersmith method more](https://hackmd.io/pP-iS2FtSWevJBcE6MEjBg). (it might be an old article, though.)

## How to choose parameters?
Coppersmith small root method is to find a root of the following type equation:

$$f(r_1,\ldots,r_n)=0 \pmod{b}\ (|r_i|<X_i)$$
for known polynomial $f(r)$ and known $N$ such that $b\ |\ N$.

The package function requires the following parameters.

- basepoly: the polynomial $f(r)$ over Zmod($N$)
- bounds: the list $[X_1,\ldots,X_n]$ whose positive integers $X_1,\ldots,X_n$
- beta: a positive floating point number $\beta$ such that $b \ge N^\beta$

For determining $\beta$, we recommend the following guideline.

- If $b$ is known, then $\beta=\log_{N}(b)$
- If $b$ is unknown but bitsize of $b$ is known, then $\beta=(\text{bitsize}(p)-1)/(\text{bitsize}(N))$

For example, $N=pq$ and $\text{bitsize}(p)=\text{bitsize}(q)$, then $\beta\simeq 0.499$.

## Completeness of the Package

Q: Can you solve many modulus multivariate polynomial systems by the package?

A: Maybe no. It seems to be hard to create general modulus multivariate polynomial solver. It depends on monomials (such as \[ $f=ax^2 + by + c$ \] `v.s.` \[ $g=ax^2 + bxy + cy + d$ \]), coefficients (such as \[ $f=ax^2-by+c$ \] `v.s.` \[ $g=ax^2-ay+c$ \]), and bounds (such as \[ $X \simeq Y$ \] `v.s.` \[ $X \simeq Y^2$ \]). Many papers handle each specific cases.

Especially, we do not recommend to use heuristic method without understanding the polynomial. Heuristic method does not estimate bound condition (unlike univariate case or linear case), so you would be confused that the solver did not output as what you expected.

Alternatively, you may use linearization strategy and/or apply rootfind_ZZ directly. For example, if you want to solve $f=ax^2+by \pmod{N}$, first define ${f_\mathbb{Z}}=ax'+by+kN \in \mathbb{Z}[x',y,k]$, then apply rootfind_ZZ (or other integer polynomial solver) with bound $[X^2, Y]$ (we search $k$ with bruteforce). If we can assume $|k|$ is small, the system will be solved. Note that Coppersmith methods does not necessarily find algebraically independent polynomial sets. So you might have to solve same polynomial system $f_\mathbb{Z}$ even if you forced to apply Coppersmith method. rootfind_ZZ.solve_root_triangulate tries to solve non-zero dimensional linear system.

Note that rootfind_ZZ does not necessarily find all roots, but only a few roots. Finding roots over integer is not easy, so you should not use the package for multiple roots included system. You can devise some methods for avoiding multiple roots. Some method might be narrowing bounds range by converting variables. Some rootfind_ZZ internal functions assume that a root exist near bounds.

## Contribution

The package must not be perfect, we want you to be reviewed and improved this. Welcome to post any issues, pull requests and forks. And let us know `test cases` from any CTF or other resources. Failed test cases may make us to improve the package.

## Reference

Some of our codes are based on the following libraries.

- defund/coppersmith (https://github.com/defund/coppersmith)

multivariate coppersmith method: `coppersmith.sage`

- josephsurinlattice-based-cryptanalysis (https://github.com/josephsurin/lattice-based-cryptanalysis)

multivariate coppersmith method: `lbc_toolkit/problems/small_roots.sage`. Some solvers for integer polynomial roots are at `lbc_toolkit/common/systems_solvers.sage`.

- jvdsn/crypto-attacks (https://github.com/jvdsn/crypto-attacks)

various coppersmith methods (including Herrmann-May, Blomer-May, Boneh-Durfee, Jochemsz-May, etc.): `shared/small_roots/*.py`. Some solvers for integer polynomial roots are at `shared/small_roots/__init__.py`.

jvdsn/crypto-attacks is distributed under:

>MIT license (https://github.com/jvdsn/crypto-attacks/blob/master/LICENSE)
>
>(c) 2020 Joachim Vandersmissen.

## Copyright

This library is distributed under Apache 2.0 License. See LICENSE.

>(C) 2023 kiona
>
>https://github.com/kionactf/coppersmith

For redistribution, just say that it has been changed and note the link for our files.
