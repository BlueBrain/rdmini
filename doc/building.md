---
title: Building rdmini
---

## Code organisation

The top-level directories in the `rdmini` project are organised as follows:

`build/`
  ~ Out-of-source build directory, with Makefile.

`demo/`
  ~ Source for demonstration and informal test code.

`doc/`
  ~ Documentation source, primarily in Pandoc-flavoured Markdown.

`gtest-1.7.0/`
  ~ Google test sources.

`rdmini/`
  ~ Source and headers for `rdmini` core functionality.

`scripts/`
  ~ Miscellaneous scripts used in the build process or for parsing results.

`test/`
  ~ Source for unit tests.

`yaml-0.1.6/`
  ~ Subset of libyaml sources.
  
## Building the executables

The standard procedure for building the demo and test executables is:

````
cd build;
make
````

Tests can be run with `make test`, which will also write the test results
in xUnit XML format to the subdirectory `test-xml`.

The `rdmini` source is written in C++11, and thus requires a C++11 capable
compiler. Everything should build 'out of the box' with gcc version 5 or
later.

The C and C++ compilers used in the Makefile are the default from the
environment, and can be specified via the `CC` and `CXX` environment
variables, or overriden on the `make` command-line. The corresponding
compiler options though are not automatically deduced, and must
be specified explicitly by setting the variable `COMPILER` in the
Makefile directly, or via the `make` command-line. Options are included
in the Makefile for the cases when `COMPILER` is `gnu`, `clang`,
`pgi` or `cray`.

Inter-source dependencies are determined automatically with `gcc`
and `clang`, but for `pgi` and `cray` these must be built
explicitly with `make depend`.

## Building the documentation

The documentation tree is built with `make doc`, which will compile
the markdown documentation into an HTML hierarchy under `www/` in
the build directory.

This build process uses [Hakyll](https://jaspervdj.be/hakyll/) and
[Pandoc](http://pandoc.org). Hakyll is a Haskell-based library
for the compilation and deployment of HTML sites; building the
documentation requires first the compilation of the site builder
program in `doc/src/site.hs`, and then running it to produce
the final HTML.

Building the documentaion thus requires a Haskell compiler and
the Hakyll package and its dependencies. The current `site.hs`
has been verified to build with [GHC](https://www.haskell.org/ghc/)
7.6.3 and Hakyll version 4.7.5.1. Installation of Hakyll
is probably best performed using Haskell
[Cabal](https://www.haskell.org/cabal/).

## Publishing the documentation

The [online documentation](https://bluebrain.github.io/rdmini/)
consists of the contents of the `gh-pages` branch of the
[BlueBrain repository](https://github.com/bluebrain/rdmini/).

The `gh-pages` make target builds the documentation and then
uses the included `commit-dir` script to commit the compiled
documentation tree into a local `gh-pages` branch. The script
constructs the commit as a merge commit between the previous
state of the `gh-pages` branch, and the latest commit on
the current branch.

Updating the online documentation to match the latest
commit on master then consists of:

1. Ensure that the local repository master and gh-pages
branches match those of `bluebrain/rdmini`.

2. Build and update the gh-pages branch locally with
`make gh-pages`.

3. Either push the updated gh-pages branch back to `bluebrain/rdmini`,
or push the branch to one's GitHub repository and submit
a pull request.


## Cleanup

Intermediate object files and temporary files will be removed
with `make clean`.

All executables, libraries, dependency files and
generated documentation will be removed with `make realclean`.



