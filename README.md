M4RIE is a library for fast arithmetic with dense matrices over GF(2^e) for 2 ≤ e ≤ 16. The name
stems from the fact that is relies heavily on M4RI. M4RIE is part of the Sage mathematics software.
M4RIE is available under the General Public License Version 2 or later (GPLv2+).

[![Build Status](https://drone.io/bitbucket.org/malb/m4rie/status.png)](https://drone.io/bitbucket.org/malb/m4rie/latest)

# Main Features #

* basic arithmetic with dense matrices over GF(2^e) (addition, equality testing, stacking,
  augmenting, sub-matrices, etc.),

* asymptotically fast $O(n^{\log_2 7})$ matrix multiplication using the Karatsuba trick due to
  Boothby and Bradshaw,

* asymptotically fast $O(n^{\log_2 7})$ matrix multiplication using “Newton-John Table” based
  multiplication & the Strassen-Winograd algorithm,

* fast row echelon form computation and matrix inversion via a “Newton-John Table” based algorithm
  using only $O(n^2)$ field multiplications and $O(n^3)$ field additions,

* asymptotically fast row echelon form computation and matrix inversion via PLE decomposition,

* asymptotically fast TRiangular System solving with Matrices (upper left, lower left),

* support for Linux, OpenSolaris and OS X (all GCC).

See [Further Reading](https://bitbucket.org/malb/m4rie/wiki/Further%20Reading) for implemented
algorithms.

# Prerequisites 

M4RIE depends heavily on [M4RI](https://bitbucket.org/malb/m4ri).

# Performance

See [Performance](http://malb.bitbucket.org/m4ri-e-website-2008-2015/performance.html).

# Install #

If you downloaded M4RIE as a compressed tarball from its website, installation instructions are
contained in the file INSTALL. If you downloaded M4RI by cloning the mainline tree at

https://bitbucket.org/malb/m4ri

you need to first run the following command:

    autoreconf --install

Then do the usual

    ./configure
    make
    make check

For details see the instructions in the file INSTALL. Note that if your system does not have
`autoreconf` installed, you can download release tarballs from
[drone.io](https://drone.io/bitbucket.org/malb/m4rie/files) which can be built and installed without
it.

# Documentation #

To build the reference manual, ensure that you have Doxygen installed. The HTML version of the
reference manual can be built as follows:

    cd src/
    doxygen

The built documentation is contained under the doc subdirectory of m4ri/. Once the HTML version is
built, you can build the PDF version as follows:

    cd doc/latex/
    make

The documentation is also available [here](http://malb.bitbucket.org/m4rie/).

# Contributors

The following people have contributed to the M4RIE library.

* **[Martin Albrecht](http://martinralbrecht.wordpress.com)**

We are grateful to **[William Stein](http://modular.math.washington.edu/)** for providing our
hosting and general infrastructure in the past.

# Citing M4RIE

If you use our libraries in a non-trivial part of your research please consider citing them as
follows:

	@manual{M4RI,
	    key          = "M4RIE",
	    author       = "Martin Albrecht",
	    organization = "The M4RIE~Team",
	    title        = "{The M4RIE Library -- Version XXXXYYZZ}",
	    year         = 2015,
	    url          = "\url{https://bitbucket.org/malb/m4rie}",
	}

and cite the appropriate publications mentioned in
[Further Reading](https://bitbucket.org/malb/m4rie/wiki/Further%20Reading).

# History

* **2015/0417** Our hosting for http://m4ri.sagemath.org at University of Washington. is
  discontinued and we’re moving everything over to https://bitbucket.org/malb/m4ri. A copy of the
  old website (except for large files) is available at
  http://malb.bitbucket.org/m4ri-e-website-2008-2015/.

* **2014/09/14** A new version of M4RI and M4RIE is available for
  [download](https://bitbucket.org/malb/m4ri/downloads). The biggest change is that `A->offset` was
  dropped. Also, various small (multicore) performance improvements were implemented. The update for
  M4RIE is to maintain compatibility with M4RI. A few improvements were implemented for the mzd_poly
  module as well.

* **2012/06/13** New versions of both M4RI and M4RIE are available for
  [download](https://bitbucket.org/malb/m4ri/downloads). A detailed changlog are available
  [here](https://bitbucket.org/malb/m4ri/wiki/M4RI-20120613) for M4RI.

* **2012/04/13** New versions of both M4RI and M4RIE are available for
  [download](https://bitbucket.org/malb/m4ri/downloads). Detailed changlogs are available
  [here](https://bitbucket.org/malb/m4ri/wiki/M4RI-20120415) for M4RI and
  [here](https://bitbucket.org/malb/m4rie/wiki/M4RIE-20120415) for M4RIE.

* **2011/12/04** New versions of both M4RI and M4RIE are available for
  [download](https://bitbucket.org/malb/m4ri/downloads). The highlight of this version for M4RI is support
  for reading and writing 1-bit PNG images. The highlight of this release of M4RIE is much improved
  performance for $4 < e \leq 8$. Detailed changlogs are available
  [here](https://bitbucket.org/malb/m4ri/wiki/M4RI-20111203) for M4RI and
  [here](https://bitbucket.org/malb/m4rie/wiki/M4RIE-20111203) for M4RIE.

* **2011/11/30** A [technical report](http://arxiv.org/abs/1111.6900) by Martin R. Albrecht is available describing the M4RIE library. In particular, Newton-John tables are introduced and our implementation of Karatsuba based matrix-matrix multiplication is described:

  > **The M4RIE library for dense linear algebra over small fields with even characteristic**
  >  
  > *Abstract:* In this work, we present the M4RIE library which implements efficient algorithms for
  > linear algebra with dense matrices over GF(2^e) for 2 ≤ e ≤ 10. As the name of the library
  > indicates, it makes heavy use of the M4RI library both directly (i.e., by calling it) and
  > indirectly (i.e., by using its concepts). We provide an open-source GPLv2+ C library for
  > efficient linear algebra over GF(2^e) for e small. In this library we implemented an idea due to
  > Bradshaw and Boothby which reduces matrix multiplication over GF(p^k) to a series of matrix
  > multiplications over GF(p). Furthermore, we propose a caching technique - Newton-John tables -
  > to avoid finite field multiplications which is inspired by Kronrod's method ("M4RM") for matrix
  > multiplication over GF(2). Using these two techniques we provide asymptotically fast triangular
  > solving with matrices (TRSM) and PLE-based Gaussian elimination. As a result, we are able to
  > significantly improve upon the state of the art in dense linear algebra over $F(2^e) with 2 ≤ e
  > ≤ 10.

