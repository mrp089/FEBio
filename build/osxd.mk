# Make include file for FEBio on Mac

CC = icpc

# Remove -DHAVE_LEVMAR and $(LEV_LIB) from LIBS if not linking with the Lourakis levmar routine.
DEF = -DPARDISO -DMKL_ISS -DHAVE_LEVMAR -DHAVE_GSL -DSVN -DHYPRE -DUSE_MPI -DHAS_MMG -DNDEBUG

# Comment this out during development
DEF += -DNDEBUG

FLG = -O3 -qopenmp -fPIC -static-intel -no-intel-extensions -std=c++11

# Pardiso solver
INTELROOT = $(subst /mkl,,$(MKLROOT))
INTEL_INC = $(INTELROOT)/compiler/include
INTEL_LIB = $(INTELROOT)/compiler/lib/
MKL_PATH = $(MKLROOT)/lib/
MKL_LIB = $(MKL_PATH)libmkl_intel_lp64.a $(MKL_PATH)libmkl_intel_thread.a $(MKL_PATH)libmkl_core.a \
	$(INTEL_LIB)libiomp5.a -pthread -lz

#Levmar library
LEV_LIB = -llevmar

GSL_LIB = /usr/local/lib/libgsl.a

#HYPRE library
HYPRE_LIB = -lHYPRE
HYPRE_PATH = /usr/local/hypre-master/src/hypre

#OMP library
OMP_PATH = /usr/local/opt/libomp

#MPI library
MPI_LIB = -lmpi

#MMG library
MMG_LIB = -lmmg3d

LIBS = -L$(FEBDIR)build/lib -L$(HYPRE_PATH)/lib $(HYPRE_LIB) $(LEV_LIB) $(GSL_LIB) $(MKL_LIB) $(MPI_LIB) $(MMG_LIB)

INC = -I$(LLVM_PATH)/include -I$(FEBDIR) -I$(FEBDIR)build/include -I$(HYPRE_PATH)/include -I$(OMP_LIB)/include

