# Paths for libraries
FFTW_DIR =/opt/homebrew/Cellar/fftw/3.3.10_2
SUITE_SPARSE_DIR = /opt/homebrew/Cellar/suite-sparse/7.10.3/include/suitesparse
SUITE_SPARSE_LIB_DIR = /opt/homebrew/Cellar/suite-sparse/7.10.3/lib
MPI_DIR = /opt/homebrew/opt/open-mpi
PETSC_DIR = /opt/homebrew/opt/petsc
PETSC_INCLUDE = $(PETSC_DIR)/include
PETSC_LIB = $(PETSC_DIR)/lib
PETSC_CFLAGS = $(shell pkg-config --cflags petsc)
PETSC_LIBS   = $(shell pkg-config --libs petsc)
METIS_DIR = /opt/homebrew/Cellar/metis/5.1.0/lib

# Target executable
all: program.exe

# Clean command
clean: 
	rm -f *.o program.exe

# Object files list
OBJ_FILES = mesh.o wlog.o gmshtools.o gelement.o gelement3D.o triplet.o Mesh_GetBndrSide.o csr.o integral.o saving.o sparse_1.o triplet_1.o vector_1.o dynamics.o problems.o Poisson.o Nonlinear.o

CFLAGS = -g -Wall -fsanitize=address $(PETSC_CFLAGS) -I$(MPI_DIR)/include -I$(SUITE_SPARSE_DIR) -I$(FFTW_DIR)/include -I$(PETSC_INCLUDE)

# Linker and compiler flags for program.exe (gcc -O3 -g -Wall program.c $(IBJ_FILES) for debugging change to gcc -O0 -g -Wall program.c $(IBJ_FILES) )
program.exe: program.c $(OBJ_FILES)
	gcc -O3 -g -Wall -fsanitize=address $(PETSC_CFLAGS) \
		program.c $(OBJ_FILES) -o program.exe \
		-Wl,-stack_size,0x10000000 \
		-I$(FFTW_DIR)/include -I$(SUITE_SPARSE_DIR) -I$(MPI_DIR)/include \
		-L$(FFTW_DIR)/lib -lfftw3 -lm \
		-L$(SUITE_SPARSE_LIB_DIR) -lsuitesparseconfig -lcholmod -lumfpack -lamd -lcolamd -lcamd \
		-L$(METIS_DIR) -lmetis \
		-L$(MPI_DIR)/lib -lmpi \
		$(PETSC_LIBS)


# Individual compilation rules for each object file, with necessary include paths
mesh.o: mesh.c mesh.h
	gcc -c $(CFLAGS) mesh.c 

wlog.o: wlog.c wlog.h
	gcc -c $(CFLAGS) wlog.c 

gmshtools.o: gmshtools.c gmshtools.h
	gcc -c $(CFLAGS) gmshtools.c 

gelement.o: gelement.c gelement.h
	gcc -c $(CFLAGS) gelement.c

gelement3D.o: gelement3D.c gelement3D.h
	gcc -c $(CFLAGS) gelement3D.c

triplet.o: triplet.c triplet.h
	gcc -c $(CFLAGS) triplet.c

Mesh_GetBndrSide.o: Mesh_GetBndrSide.c Mesh_GetBndrSide.h
	gcc -c $(CFLAGS) Mesh_GetBndrSide.c 

csr.o: csr.c csr.h
	gcc -c $(CFLAGS) csr.c

integral.o: integral.c integral.h
	gcc -c $(CFLAGS) integral.c

saving.o: saving.c saving.h
	gcc -c $(CFLAGS) saving.c

triplet_1.o: triplet_1.c triplet_1.h
	gcc -c $(CFLAGS) triplet_1.c

sparse_1.o: sparse_1.c sparse_1.h
	gcc -c $(CFLAGS) sparse_1.c

vector_1.o: vector_1.c vector_1.h
	gcc -c $(CFLAGS) vector_1.c

dynamics.o: dynamics.c dynamics.h
	gcc -c $(CFLAGS) dynamics.c

problems.o: problems.c problems.h
	gcc -c $(CFLAGS) problems.c

Poisson.o: Poisson.c Poisson.h
	gcc -c $(CFLAGS) Poisson.c

Nonlinear.o: Nonlinear.c Nonlinear.h
	gcc -c $(CFLAGS) Nonlinear.c
