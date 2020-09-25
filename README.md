# General Relativistic Atomic Structure Package

![Tests][tests-badge]
[![][doxygen-badge]][doxygen-url]
[![][manual-badge]][manual-pdf]

The General Relativistic Atomic Structure Package (GRASP) is a set of Fortran 90
programs for performing fully-relativistic electron structure calculations of
atoms.

## Installation

> **Please note:**
> The installation instructions here are for the _development version_ on the
> `master` branch.
>
> To install the _latest published release_ (2018-12-03), go to the
> ["Releases" page](https://github.com/compas/grasp/releases/tag/2018-12-03),
> download the tarball from there and refer to the instructions in the README in
> the tarball.

To compile and install GRASP, first clone this Git repository:

```sh
git clone https://github.com/compas/grasp.git
```

There are two ways to build GRASP: either via [CMake](https://cmake.org/) or via the
`Makefile`s in the source tree. Either works and you end up with the GRASP binaries in the
`bin/` directory.

CMake is the recommended way to build GRASP. The `Makefile`-based workflow is still there to
make smoother to transition from `Makefile`s to a modern build system.

### CMake-based build

The first step with CMake is to create a separate out-of-source build directory. The
`configure.sh` script can do that for you:

```sh
cd grasp/ && ./configure.sh
```

This will create a `build/` directory with the default _Release_ build
configuration. However, `configure.sh` is just a simple wrapper around a `cmake`
call and if you need more control over the build, you can always invoke `cmake`
yourself (see [CMake documentation](https://cmake.org/documentation/) for more
information).

To then compile GRASP, you need to go into the out-of-source build directory and
simply call `make`:

```sh
cd build/ && make install
```

Remarks:

* Running `make install` instructs CMake to actually _install_ the resulting binaries into
  the conventional `bin/` directory at the root of the repository.

  When you run just `make`, the resulting binaries will end up under the `build/` directory
  (specifically in `build/bin/`). This is useful when developing and debugging, as it allows
  you to compile many versions of the binaries from the same source tree with different
  compilation options (e.g. build with debug symbols enabled) by using several out of source
  build directories.

* With CMake, GRASP also supports parallel builds, which can be enabled by passing the `-j`
  option to `make` (e.g. `make -j4 install` to build with four processes).

* The CMake-based build allows running the (non-comprehensive) test suite by calling `ctest`
  in the `build/` directory. The configuration and source files for the tests are under
  `test/`/

### `Makefile`-based build

The legacy `Makefile`-based build can be performed by simply calling the `make` in the top
level directory:

```sh
make
```

In this case, the compilation of each of the libraries and programs happens in their
respective directory under `src/` and the build artifacts are stored in the source tree.
The resulting binaries and libraries will directly get installed under the `bin/` and `lib/`
directories.

To build a specific library or binary you can pass the path to the source directory as the
Make target:

```sh
# build libmod
make src/lib/libmod
# build the rci_mpi binary
make src/appl/rci90_mpi
```

Note that any necessary library dependencies will also get built automatically.

**WARNING:** the `Makefile`s do not know about the dependencies between the source files, so
parallel builds (i.e. calling `make` with the `-j` option) does not work.

#### Customizing the build

By default the `Makefile` is designed to use `gfortran`. The variables affecting GRASP
builds are defined and documented at the beginning of the `Makefile`.

For the user it should never be necessary to modify the `Makefile` itself. Rather, a
`Make.user` file can be create next to the main `Makefile` where the build variables can be
overridden. E.g. to use the Intel Fortran compiler instead, you may want to create the
following `Make.user` file:

```make
export FC = ifort
export FC_FLAGS = -O2 -save
export FC_LD = -mkl=sequential
export FC_MPI = mpiifort
```

As another example, to set up a linker search path for the BLAS or LAPACK libraries, you can
set up `Make.user` as follows:

```make
export FC_LD = -L /path/to/blas
```

## About GRASP

This version of GRASP is a major revision of the previous GRASP2K package by [P.
Jonsson, G. Gaigalas, J. Bieron, C. Froese Fischer, and I.P. Grant Computer
Physics Communication, 184, 2197 - 2203 (2013)][grasp2k-2013] written in FORTRAN
77 style with COMMON and using Cray pointers for memory management.  The present
version is a FORTRAN95 translation using standard FORTRAN for memory management.
In addition, COMMONS have been replaced with MODULES, with some COMMONS merged.
Some algorithms have been changed to improve performance for large cases and
efficiently.

The previous package, was an extension and modification of GRASP92 by [Farid
Parpia, Charlotte Froese Fischer, and Ian Grant. Computer Physics Communication,
94, 249-271 (1996)][grasp92-1996].

This version of GRASP has been published in:

> C. Froese Fischer, G. Gaigalas, P. Jönsson, J. Bieroń,
> "GRASP2018 — a Fortran 95 version of the General Relativistic Atomic Structure Package",
> Computer Physics Communications, 237, 184-187 (2018),
> https://doi.org/10.1016/j.cpc.2018.10.032

Development of this package was performed largely by:

|                           | email                         |
| ------------------------- | ------------------------------|
| Charlotte Froese Fischer  | cff@cs.ubc.ca                 |
| Gediminas Gaigalas        | Gediminas.Gaigalas@tfai.vu.lt |
| Per Jönsson               | per.jonsson@mau.se            |
| Jacek Bieron              | jacek.bieron@uj.edu.pl        |

Please contact one of these authors if you have questions.

Supporters include:

|                           | email                         |
| ------------------------- | ------------------------------|
| Jörgen Ekman              | jorgen.ekman@mah.se           |
| Ian Grant                 | iangrant15@btinternet.com     |



## Structure of the Package

The package has the structure shown below where executables, after successful
compilation, reside in the `bin` directory. Compiled libraries are in the `lib`
directory. Scripts for example runs and case studies are in folders under
`grasptest`. Source code is in the `src` directory and divided into applications
in the `appl` directory, libraries in the `lib` directory and tools in the
`tool` directory.

```
   |-bin
   |-grasptest
   |---case1
   |-----script
   |---case1_mpi
   |-----script
   |-----tmp_mpi
   |---case2
   |-----script
   |---case2_mpi
   |-----script
   |-----tmp_mpi
   |---case3
   |-----script
   |---example1
   |-----script
   |---example2
   |-----script
   |---example3
   |-----script
   |---example4
   |-----script
   |-------tmp_mpi
   |---example5
   |-----script
   |-lib
   |-src
   |---appl
   |-----HF
   |-----jj2lsj90
   |-----jjgen90
   |-----rangular90
   |-----rangular90_mpi
   |-----rbiotransform90
   |-----rbiotransform90_mpi
   |-----rci90
   |-----rci90_mpi
   |-----rcsfgenerate90
   |-----rcsfinteract90
   |-----rcsfzerofirst90
   |-----rhfs90
   |-----rmcdhf90
   |-----rmcdhf90_mpi
   |-----rnucleus90
   |-----rtransition90
   |-----rtransition90_mpi
   |-----rwfnestimate90
   |-----sms90
   |---lib
   |-----lib9290
   |-----libdvd90
   |-----libmcp90
   |-----libmod
   |-----librang90
   |-----mpi90
   |---tool
```


## Program Guide and Compilation

The software is distributed with a practical guide to [GRASP2018 in PDF-format
(click here to download)][manual-pdf]. The guide, which is under Creative
Commons Attribution 4.0 International (CC BY 4.0) license, contains full
information on how to compile and install the package.


## Acknowledgements

This work was supported by the Chemical Sciences, Geosciences and Biosciences
Division, Office of Basic Energy Sciences, Office of Science, U.S. Department of
Energy who made the Pacific Sierra translator available and the National
Institute of Standards and Technology. Computer resources were made available by
Compute Canada.  CFF had research support from the Canadian NSERC Discovery
Grant 2017-03851.  JB acknowledges financial support of the European Regional
Development Fund in the framework of the Polish Innovation Economy Operational
Program (Contract No. POIG.02.01.00-12-023/08).


## Copyright & license

The code in this repository is distributed under the [MIT license](LICENSE).
The accompanying guide  "A practical guide to GRASP2018" is licensed separately
under [the CC-BY-4.0 (Creative Commons Attribution 4.0 International) license][cc-by].

[manual-pdf]: https://github.com/compas/grasp2018/releases/download/2018-12-03/GRASP2018-manual.pdf
[manual-badge]: https://img.shields.io/badge/manual-pdf-blue.svg
[doxygen-url]: https://compas.github.io/grasp/
[doxygen-badge]: https://img.shields.io/badge/documentation-doxygen-blue.svg
[tests-badge]: https://github.com/compas/grasp/workflows/Tests/badge.svg
[grasp92-1996]: https://doi.org/10.1016/0010-4655(95)00136-0
[grasp2k-2013]: https://doi.org/10.1016/j.cpc.2013.02.016
[cc-by]: https://creativecommons.org/licenses/by/4.0/legalcode
