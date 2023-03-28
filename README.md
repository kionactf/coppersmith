# coppersmith

Coppersmith method (solving polynomial equation over composite modulus on small bounds)

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

## Contribution

The package must not be perfect, we want you to be reviewed and improved this. Welcome to post any issues, pull requests and forks. And let us know `test cases` from any CTF or other resources.

## Reference

Some of our codes are based on the following libraries.

- defund/coppersmith (https://github.com/defund/coppersmith)

multivariate coppersmith method: `coppersmith.sage`

- josephsurinlattice-based-cryptanalysis (https://github.com/josephsurin/lattice-based-cryptanalysis)

multivariate coppersmith method: `lbc_toolkit/problems/small_roots.sage`. Some solvers for integer polynomial roots are at `lbc_toolkit/common/systems_solvers.sage`.

- jvdsn/crypto-attacks (https://github.com/jvdsn/crypto-attacks)

various coppersmith methods (including Herrmann-May, Blomer-May, Boneh-Durfee, Jochemsz-May, etc.): `shared/small_roots/*.py`. Some solvers for integer polynomial roots are at `shared/small_roots/__init__.py`.

jvdsn/crypto-attacks is distributed under:

MIT license (https://github.com/jvdsn/crypto-attacks/blob/master/LICENSE)

(c) 2020 Joachim Vandersmissen.

## Copyright

This library is distributed under Apache 2.0 License. See LICENSE.

(C) 2023 kiona<kiona.ctf@gmail.com>

https://github.com/kionactf/coppersmith

For redistribution, just say that it has been changed and note the link for our files.
