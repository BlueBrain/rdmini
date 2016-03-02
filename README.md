# rdmini

Test platform for stochastic reaction-diffusion implementations.

## License

`rdmini` source and documentation is licensed under the
[GPL v2 license](http://www.gnu.org/licenses/gpl-2.0.txt)
(see the `LICENSE.txt` file.)

The CSL style file `doc/style.csl` is the APA 5th edition style from the
[CSL style repository](http://citationstyles.org/styles/)
and is licensed under the Creative Commons
[Attribution-ShareAlike license](http://creativecommons.org/licenses/by-sa/3.0/).

The `libyaml` YAML parsing library is licensed under the MIT
license; refer to the included `yaml-0.1.6/LICENSE` file for details.

`gtest` is licensed under the Revised BSD license; refer to the
included `gtest-1.7.0/LICENSE` file for details.

## Building

The `build` directory contains a Makefile that will build documentation
and test and demo executables in this directory.

**Build all executables**
 
    cd build
    make

**Build and run tests**

    cd build
    make test

**Build HTML documentation**

    cd build 
    make doc

Temporary build artefacts are removed with `make clean`, and all
built objects are removed with `make realclean`.

The executables and tests will (should) build without any
tinkering with gcc version 5 or later. Other compilers may
require modification of the Makefile; more information can
be found in the [documentation](https://bluebrain.github.io/rdmini/).

# Demo software

The main demonstration program is `demo_sim`, which will run
a model for a fixed number of SSA steps or simulation time
and report the results in a CSV format suitable for further
analysis or processing. The `demo/models.yaml` file contains
a number of demonstration models for use with the `rdmini`
code; refer to `doc/demo_models.md` for details.
