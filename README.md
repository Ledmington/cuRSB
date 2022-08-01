# cuRSB - CUDA Recursive Sparse Blocks
This project is a CUDA extension of the `librsb-1.3.0.1` C/C++ API. The website of the original project is available [here](http://librsb.sourceforge.net/). The original source code of the `1.3.0.1` version is available [here](https://sourceforge.net/projects/librsb/).

## Naming conventions
This project tries to follow the same function/file naming conventions of the `librsb` project, as defined by the original author Michele Martone in the `./librsb/README[.md]` file, which are briefly reported below:
 - all functions named `rsb_*` or `RSB_*` are intended for external users and actually define the `librsb` interface
 - all functions named `rsb__*` or `RSB__*` are considered internals and therefore are not meant for external users

The naming conventions for the `cuRSB` project are reported below:
 - all functions named `rsb_cuda_*` are intended for external users and actually define the `cuRSB` interface
 - all functions named `rsb_cuda__*` are considered internals and therefore are not meant for external users