# Fast computation of the Baker-Campbell-Hausdorff and similar series
We provide an efficient  program written in the C programming language for the fast computation
of the Baker-Campbell-Hausdorff (BCH) and similar Lie series.
The Lie series can be represented in the Lyndon basis, in the
classical Hall basis, or in the right-normed basis of 
E.S. Chibrikov.  In the Lyndon basis,
which proves to be particularly efficient for this purpose,
the computation of 111013 coefficients for the BCH series up to terms of degree 20
takes less than half a second on an ordinary personal computer and requires negligible 11MB of memory.
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

## Documentation
Work in progress...

## Examples
```
$ ./bch
+1/1*A+1/1*B+1/2*[A,B]+1/12*[A,[A,B]]+1/12*[[A,B],B]+1/24*[A,[[A,B],B]]-1/720*[A
,[A,[A,[A,B]]]]+1/180*[A,[A,[[A,B],B]]]+1/360*[[A,[A,B]],[A,B]]+1/180*[A,[[[A,B]
,B],B]]+1/120*[[A,B],[[A,B],B]]-1/720*[[[[A,B],B],B],B]

```
### Output in tabular form:
```
$ ./bch table_output=1
0       1       0       0       1/1
1       1       1       0       1/1
2       2       0       1       1/2
3       3       0       2       1/12
4       3       2       1       1/12
5       4       0       3       0/1
6       4       0       4       1/24
7       4       4       1       0/1
8       5       0       5       -1/720
9       5       0       6       1/180
10      5       3       2       1/360
11      5       0       7       1/180
12      5       2       4       1/120
13      5       7       1       -1/720
```

### Computation up to higher degree:
```
$ ./bch N=8
+1/1*A+1/1*B+1/2*[A,B]+1/12*[A,[A,B]]+1/12*[[A,B],B]+1/24*[A,[[A,B],B]]-1/720*[A
,[A,[A,[A,B]]]]+1/180*[A,[A,[[A,B],B]]]+1/360*[[A,[A,B]],[A,B]]+1/180*[A,[[[A,B]
,B],B]]+1/120*[[A,B],[[A,B],B]]-1/720*[[[[A,B],B],B],B]-1/1440*[A,[A,[A,[[A,B],B
]]]]+1/720*[A,[[A,[A,B]],[A,B]]]+1/360*[A,[A,[[[A,B],B],B]]]+1/240*[A,[[A,B],[[A
,B],B]]]-1/1440*[A,[[[[A,B],B],B],B]]+1/30240*[A,[A,[A,[A,[A,[A,B]]]]]]-1/5040*[
A,[A,[A,[A,[[A,B],B]]]]]+1/10080*[A,[A,[[A,[A,B]],[A,B]]]]+1/3780*[A,[A,[A,[[[A,
B],B],B]]]]+1/10080*[[A,[A,[A,B]]],[A,[A,B]]]+1/1680*[A,[A,[[A,B],[[A,B],B]]]]+1
/1260*[A,[[A,[[A,B],B]],[A,B]]]+1/3780*[A,[A,[[[[A,B],B],B],B]]]+1/2016*[[A,[A,B
]],[A,[[A,B],B]]]-1/5040*[[[A,[A,B]],[A,B]],[A,B]]+13/15120*[A,[[A,B],[[[A,B],B]
,B]]]+1/10080*[[A,[[A,B],B]],[[A,B],B]]-1/1512*[[A,[[[A,B],B],B]],[A,B]]-1/5040*
[A,[[[[[A,B],B],B],B],B]]+1/1260*[[A,B],[[A,B],[[A,B],B]]]-1/2016*[[A,B],[[[[A,B
],B],B],B]]-1/5040*[[[A,B],B],[[[A,B],B],B]]+1/30240*[[[[[[A,B],B],B],B],B],B]+1
/60480*[A,[A,[A,[A,[A,[[A,B],B]]]]]]-1/15120*[A,[A,[A,[[A,[A,B]],[A,B]]]]]-1/100
80*[A,[A,[A,[A,[[[A,B],B],B]]]]]+1/20160*[A,[[A,[A,[A,B]]],[A,[A,B]]]]-1/20160*[
A,[A,[A,[[A,B],[[A,B],B]]]]]+1/2520*[A,[A,[[A,[[A,B],B]],[A,B]]]]+23/120960*[A,[
A,[A,[[[[A,B],B],B],B]]]]+1/4032*[A,[[A,[A,B]],[A,[[A,B],B]]]]-1/10080*[A,[[[A,[
A,B]],[A,B]],[A,B]]]+13/30240*[A,[A,[[A,B],[[[A,B],B],B]]]]+1/20160*[A,[[A,[[A,B
],B]],[[A,B],B]]]-1/3024*[A,[[A,[[[A,B],B],B]],[A,B]]]-1/10080*[A,[A,[[[[[A,B],B
],B],B],B]]]+1/2520*[A,[[A,B],[[A,B],[[A,B],B]]]]-1/4032*[A,[[A,B],[[[[A,B],B],B
],B]]]-1/10080*[A,[[[A,B],B],[[[A,B],B],B]]]+1/60480*[A,[[[[[[A,B],B],B],B],B],B
]]
```

### Output using right-normed Chibrikov basis:
```
$ ./bch basis=1
+1/1*A+1/1*B-1/2*[B,A]-1/12*[A,[B,A]]+1/12*[B,[B,A]]+1/24*[B,[A,[B,A]]]+1/720*[A
,[A,[A,[B,A]]]]-1/360*[B,[A,[A,[B,A]]]]+1/120*[A,[B,[A,[B,A]]]]-1/120*[B,[B,[A,[
B,A]]]]+1/360*[A,[B,[B,[B,A]]]]-1/720*[B,[B,[B,[B,A]]]]
```

### Compute symmetric BCH formula:
``` 
$./bch expression=1
+1/1*A+1/1*B-1/24*[A,[A,B]]+1/12*[[A,B],B]+7/5760*[A,[A,[A,[A,B]]]]-7/1440*[A,[A
,[[A,B],B]]]+1/360*[[A,[A,B]],[A,B]]+1/180*[A,[[[A,B],B],B]]+1/120*[[A,B],[[A,B]
,B]]-1/720*[[[[A,B],B],B],B]
```

