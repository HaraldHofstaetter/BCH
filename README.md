# Fast computation of the Baker-Campbell-Hausdorff and similar series
An efficient  program written in the C programming language for the fast computation
of the Baker-Campbell-Hausdorff (BCH) and similar Lie series is provided.
The Lie series can be represented in the Lyndon basis, in the
classical Hall basis, or in the right-normed basis of 
E.S. Chibrikov.  In the Lyndon basis,
which proves to be particularly efficient for this purpose,
the computation of 111013 coefficients for the BCH series up to terms of degree 20
takes less than half a second on a common personal computer and requires negligible 11MB of memory.
Up to terms of degree 30, which is the maximum degree the program can handle, 
the computation of 74248451 coefficients takes 55 hours but still requires only a modest 5.5GB of  memory.

## Installation
On a Unix-like system with `gcc` compiler available just type
```
$ make
```
in a directory containing the source code, which causes 
the shared library `libbch.so` and the executable `bch`
to be created.
To use a different compiler, the `Makefile` has to be adapted accordingly.

