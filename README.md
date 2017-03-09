# BCMTools

* Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
All rights reserved.

* Copyright (c) 2012-2016 Advanced Institute for Computational Science (AICS), RIKEN.
All rights reserved.

* Copyright (c) 2017 Research Institute for Information Technology (RIIT), Kyushu University.
All rights reserved.


## OUTLINE



## SOFTWARE REQUIREMENT
- Cmake
- MPI library
- TextParser (required to build examples)
- Polylib
- HDF library (option)
- Silo library (option)


## INGREDIENTS
~~~
BCMViewer/        Simple viewer for checking
cmake/            cmake modules
doc/              Documents
doxygen/          doxygen files
example/          Example source codes
include/          Header files
src/              Source codes
Utils/            Utility codes
ChangeLog.md      History of development
License.txt       License to apply
Readme.md         This document, including the description of build
~~~


## HOW TO BUILD

### Build

~~~
$ export CC=c_compiler CXX=c++_compiler F90=fortran_compiler FC=fortran_compiler
$ export BCM_HOME=/hogehoge
$ mkdir build
$ cd build
$ cmake [options] ..
$ make
$ sudo make install
~~~


### Options

`-D INSTALL_DIR=` *Install_directory*

>  Specify the directory that this library will be installed. Built library is installed at `install_directory/lib` and the header files are placed at `install_directory/include`. The default install directory is `/usr/local/BCMTools`.

`-D with_MPI=` {yes | no}

>  If you use an MPI library, specify `with_MPI=yes`, the default is yes.

`-D with_TP =` *TextParser_directory*

> IBuilding examples, specify the directory path that TextParser is installed.

`-D with_PL=` *Polylib_directory*

> Specify the directory path that Polylib is installed.

`-D with_SILO=` *silo_directory*

> Specify the directory path that Cutlib is installed, if required.

`-D with_CUT=` *Cutlib_directory*

> Specify the directory path that Cutlib is installed, if required.

`-D with_HDF5=` *HDF5_directory*

> Specify the directory path that HDF5 is installed, if required.

`-D with_example=` {no | yes}

>  This option turns on compiling sample codes. The default is no.

`-D enable_LARGE_BLOCK=` {no | yes}

>  Specify optimization of buffer size to accelerate file I/O performance. The default is no.
>  BCMブロックが大きく（ブロック内セル数が多い），数が少ない場合に高速化が期待できる．BCMブロックが小さく，数が多い場合は，逆効果．

`-D enable_OPENMP=` {no | yes}

> Enable OpenMP directives.

`-D real_type=` {float | double}

>  Specify the type of floating point. If this option is omitted, the default is float.

The default compiler options are described in `cmake/CompilerOptionSelector.cmake` file. See BUILD OPTION section in CMakeLists.txt in detail.



## Configure Examples

`$ export BCM_HOME=hogehoge`

In following exsmples, assuming that TextParser and Polylib are installed under the BCM_HOME directory. If not, please specify applicable directory paths.

### INTEL/GNU compiler

~~~
$ cmake -DINSTALL_DIR=${BCM_HOME}/BCMTools -Denable_OPENMP=yes -Dwith_MPI=yes -Dwith_example=yes -Dreal_type=double -Dwith_TP=${BCM_HOME}/TextParser -Dwith_PL=${BCM_HOME}/Polylib  -Denable_LARGE_BLOCK=no ..
~~~


### FUJITSU compiler / FX10, FX100, K on login nodes (Cross compilation)

~~~
$ cmake -DINSTALL_DIR=${BCM_HOME}/BCMTools \
            -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_fx10.cmake \
            -Denable_OPENMP=yes \
            -Dwith_MPI=yes \
            -Dwith_example=no \
            -Dreal_type=double \
            -Dwith_TP=${BCM_HOME}/TextParser \
            -Dwith_PL=${BCM_HOME}/Polylib \
            -Denable_LARGE_BLOCK=no ..

$ cmake -DINSTALL_DIR=${BCM_HOME}/BCMTools \
            -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_fx100.cmake \
            -Denable_OPENMP=yes \
            -Dwith_MPI=yes \
            -Dwith_example=no \
            -Dreal_type=double \
            -Dwith_TP=${BCM_HOME}/TextParser \
            -Dwith_PL=${BCM_HOME}/Polylib \
            -Denable_LARGE_BLOCK=no ..

$ cmake -DINSTALL_DIR=${BCM_HOME}/BCMTools \
            -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_K.cmake \
            -Denable_OPENMP=yes \
            -Dwith_MPI=yes \
            -Dwith_example=no \
            -Dreal_type=double \
            -Dwith_TP=${BCM_HOME}/TextParser \
            -Dwith_PL=${BCM_HOME}/Polylib \
            -Denable_LARGE_BLOCK=no ..
~~~


##### Note
- On Fujitsu machines(fx10, K, fx100), confirm appropriate directrory path for compiler environment.
- Before building, execute following command for clean. `$ make distclean`


## EXAMPLES

* If you specify the test option by `-Denable_example=yes`, you can
execute the intrinsic tests by;

	`$ make test` or `$ ctest`

* The detailed results are written in `BUILD/Testing/Temporary/LastTest.log` file.
Meanwhile, the summary is displayed for stdout.




## CONTRIBUTORS

* Kenji     Ono      _keno@{cc.kyushu-u.ac, riken}.jp_
* Soichiro  Suzuki   soichiro.suzuki@riken.jp
* Junya     Onishi   jonishi@iis.u-tokyo.ac.jp
* Syoyo     Fujita
