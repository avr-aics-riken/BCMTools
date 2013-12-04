
DEFINES = -DDEBUG

#CXX = mpic++
#OPTIONS = $(DEFINES) -DTIMING -O3 -ftree-vectorize -g -Wall -Wno-sign-compare -Wunknown-pragmas -fopenmp

CXX = icpc
OPTIONS = $(DEFINES) -DTIMING -O3 -openmp -I/usr/lib/openmpi/include -L/usr/lib/openmpi/lib -lmpi -lmpi_cxx -lmpi_f90 -lmpi_f77 -lifport -lifcore

#CXX = mpiFCCpx
#OPTIONS = $(DEFINES) -Xg -Kfast,openmp,ocl,optmsg=2 -Nsrc


#可視化用出力を利用するには，次の設定も必要
#  SILO_INC: Siloライブラリのインクルードファイルディレクトリ指定
#  SILO_LIB: Siloライブラリのリンクオプション

#SILO_INC = -I$(HOME)/silo-4.8-bsd/include
#SILO_LIB = -L$(HOME)/silo-4.8-bsd/lib -L$(HOME)/hdf5-1.8.9/hdf5/lib -lsiloh5 -lhdf5 -lz

#SILO_INC = -I$(HOME)/include
#SILO_LIB = -L$(HOME)/lib -lsiloh5 -lhdf5 -lz

#SILO_INC = -I$(HOME)/FX10/include
#SILO_LIB = -L$(HOME)/FX10/lib -lsiloh5 -lhdf5 -lz

#SILO_INC = -I$(HOME)/K/include
#SILO_LIB = -L$(HOME)/K/lib -lsiloh5 -lhdf5 -lz

SILO_INC = -I/usr/include
SILO_LIB = -L/usr/lib -lsiloh5 -lhdf5 -lz
