# For GCC
CXX        = g++
OPT_FLAGS  = -O3 
OMP_FLAGS  = -fopenmp
MPI_DIR    = /usr/lib/openmpi
DEFINES    = #-DDEBUG
DEFINES   += -D_BLOCK_IS_LARGE
CFLAGS     = $(OPT_FLAGS) $(OMP_FLAGS) $(DEFINES) -I$(MPI_DIR)/include
AR     = ar cru
RANLIB = ranlib
RM     = \rm -f

# For K/FX10
CXX        = mpiFCCpx
OPT_FLAGS  = -Kfast,ocl,preex,simd=2,uxsimd,array_private,parallel,optmsg=1 -V -Nsrc -Xg
OMP_FLAGS  = -Kopenmp
MPI_DIR    =
DEFINES    = #-DDEBUG
DEFINES   += -D_BLOCK_IS_LARGE
CFLAGS     = $(OPT_FLAGS) $(OMP_FLAGS) $(DEFINES) -I$(MPI_DIR)/include
AR     = ar cr
RANLIB = ranlib
RM     = \rm -f

# For Intel compiler
CXX        = icpc
OPT_FLAGS  = -O3 #-xAVX
OMP_FLAGS  = -openmp
MPI_DIR    = /usr/opt/openmpi
DEFINES    = #-DDEBUG
DEFINES   += #-D_BLOCK_IS_LARGE
CFLAGS     = $(OPT_FLAGS) $(OMP_FLAGS) $(DEFINES) -I$(MPI_DIR)/include
AR     = ar cru
RANLIB = ranlib
RM     = \rm -f

# マクロ: _BLOCK_IS_LARGE_
#  BCMブロックが大きく（ブロック内セル数が多い），数が少ない場合に，
#  高速化が期待できる．
#  BCMブロックが小さく，数が多い場合は，逆効果．
