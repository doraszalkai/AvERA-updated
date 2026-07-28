USING_CUDA = YES
SINGLE_PRECISION = 0

CXX = g++
#CXX = $(HOME)/avera_env/bin/x86_64-conda-linux-gnu-g++
BUILD_NUMBER_LDFLAGS  = -Xlinker --defsym -Xlinker __BUILD_DATE=$$(date +'%Y%m%d')
#CFLAGS = -lm -Wall -ansi -pedantic -O3 -std=c++11 -fopenmp #-ffast-math
#CFLAGS = -lm -Wall -ansi -pedantic -O0 -fno-fast-math -ffloat-store -std=c++11 -fopenmp
CFLAGS = -lm -Wall -ansi -pedantic -O0 -fno-fast-math -ffloat-store -std=c++11 -fopenmp -DDOUBLE -DNO_DIM=3
#voro++ library and DTFE library
VORO=${CURDIR}/../voro++-0.4.6/src
DTFE=${CURDIR}/../DTFE_1.1.1
E_INC=-I$(VORO)
D_INC=-I$(DTFE)/DTFE_include
E_LIB=-L$(VORO)
D_LIB=-L$(DTFE)/DTFE_lib

#LDFLAGS += $(D_LIB) -Wl,-R$(DTFE)/DTFE_lib

ifeq ($(USING_CUDA), YES)
# Location of the CUDA Toolkit
CUDA_PATH       ?= /usr


NVCC = $(CUDA_PATH)/bin/nvcc -ccbin
#CUDAFLAGS = -lineinfo --compiler-options --std=c++11 --compiler-options -Wall --compiler-options -ansi -O3 -Xcompiler -fopenmp -lm
#CUDAFLAGS = -lineinfo --compiler-options --std=c++11 --compiler-options -Wall --compiler-options -ansi -O0 -Xcompiler -fopenmp -lm
CUDAFLAGS = -lineinfo --compiler-options --std=c++11 --compiler-options -Wall --compiler-options -ansi -O0 -Xcompiler -fopenmp -lm -DDOUBLE -DNO_DIM=3
CUDA_INC = -I$(CUDA_PATH)/include

LDFLAGS += $(D_LIB) -Wl,-R$(DTFE)/DTFE_lib
CUDALDFLAGS += $(D_LIB) -Xlinker "-R$(DTFE)/DTFE_lib"

## DTFE-hez
SRC = main.cc forces_cuda.cu ewald_space.cc step.cc read_paramfile.cc friedmann_solver.cc read_gadget_ic.cc nonis_friedmann.cc nonis_friedmann_voronoi.cc DTFE_density.cc
OBJ = main.o forces_cuda.o ewald_space.o step.o read_paramfile.o friedmann_solver.o read_gadget_ic.o nonis_friedmann.o nonis_friedmann_voronoi.o DTFE_density.o

## CIC-hez
#SRC = main.cc forces_cuda.cu ewald_space.cc step.cc read_paramfile.cc friedmann_solver.cc read_gadget_ic.cc nonis_friedmann.cc nonis_friedmann_voronoi.cc CIC_density.cc
#OBJ = main.o forces_cuda.o ewald_space.o step.o read_paramfile.o friedmann_solver.o read_gadget_ic.o nonis_friedmann.o nonis_friedmann_voronoi.o CIC_density.o


DEPS = global_variables.h $(VORO)/voro++.hh $(DTFE)/DTFE.h
PROG = CCLEA_CUDA

ifeq ($(SINGLE_PRECISION), 1)
$(PROG): $(OBJ)
	$(NVCC) $(CXX) $(CUDAFLAGS) $(CUDALDFLAGS) $(CUDA_INC) $(BUILD_NUMBER_LDFLAGS) $(E_INC) $(D_INC) $(E_LIB) $(D_LIB) -DUSE_SINGLE_PRECISION=$(SINGLE_PRECISION) -o $(PROG) $(OBJ) -lvoro++ -lDTFE -lCGAL -lgmp -lmpfr -lgsl -lgslcblas -lboost_program_options -lboost_system -lboost_filesystem -lboost_thread -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lcurl
%.o: %.cc
	$(CXX) $(CFLAGS) $(CUDA_INC) $(LDFLAGS) $(E_INC) $(D_INC) $(E_LIB) $(D_LIB) -DUSE_SINGLE_PRECISION=$(SINGLE_PRECISION) -o $@ -c $< -lDTFE -lCGAL -lgmp -lmpfr -lgsl -lgslcblas -lboost_program_options -lboost_system -lboost_filesystem -lboost_thread -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lcurl

%.o: %.cu
	$(NVCC) $(CXX) $(CUDA_INC) $(BUILD_NUMBER_LDFLAGS) $(E_INC) $(D_INC) $(E_LIB) $(D_LIB) -DUSE_SINGLE_PRECISION=$(SINGLE_PRECISION) -o $@ -c $<
else
$(PROG): $(OBJ)
	$(NVCC) $(CXX) $(CUDAFLAGS) $(CUDALDFLAGS) $(CUDA_INC) $(BUILD_NUMBER_LDFLAGS) $(E_INC) $(D_INC) $(E_LIB) $(D_LIB) -o $(PROG) $(OBJ) -lvoro++ -lDTFE -lCGAL -lgmp -lmpfr -lgsl -lgslcblas -lboost_program_options -lboost_system -lboost_filesystem -lboost_thread -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lcurl
%.o: %.cc
	$(CXX) $(CFLAGS) $(CUDA_INC) $(LDFLAGS) $(E_INC) $(D_INC) $(E_LIB) $(D_LIB) -o $@ -c $< -lDTFE -lCGAL -lgmp -lmpfr -lgsl -lgslcblas -lboost_program_options -lboost_system -lboost_filesystem -lboost_thread -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lcurl

%.o: %.cu
	$(NVCC) $(CXX) $(CUDA_INC) $(BUILD_NUMBER_LDFLAGS) $(E_INC) $(D_INC) $(E_LIB) $(D_LIB) -o $@ -c $<
endif
clean:
	rm -fv $(OBJ) $(PROG)

else

LDFLAGS += $(D_LIB) -Wl,-R$(DTFE)/DTFE_lib

SRC = main.cc forces.cc ewald_space.cc step.cc read_paramfile.cc friedmann_solver.cc read_gadget_ic.cc nonis_friedmann.cc nonis_friedmann_voronoi.cc DTFE_density.cc
OBJ = main.o forces.o ewald_space.o step.o read_paramfile.o friedmann_solver.o read_gadget_ic.o nonis_friedmann.o nonis_friedmann_voronoi.o DTFE_density.o


#SRC = main.cc forces.cc ewald_space.cc step.cc read_paramfile.cc friedmann_solver.cc read_gadget_ic.cc nonis_friedmann.cc nonis_friedmann_voronoi.cc CIC_density.cc
#OBJ = main.o forces.o ewald_space.o step.o read_paramfile.o friedmann_solver.o read_gadget_ic.o nonis_friedmann.o nonis_friedmann_voronoi.o CIC_density.o

##Ezekben meg valamiert a forces_cuda.o bent volt, pedig cuda nelkuli a forditas:
#Ez DTFE-hez
#SRC = main.cc forces.cc ewald_space.cc step.cc read_paramfile.cc friedmann_solver.cc read_gadget_ic.cc nonis_friedmann.cc nonis_friedmann_voronoi.cc DTFE_density.cc
#OBJ = main.o forces_cuda.o ewald_space.o step.o read_paramfile.o friedmann_solver.o read_gadget_ic.o nonis_friedmann.o nonis_friedmann_voronoi.o DTFE_density.o

#### Ez CIC-hez
#SRC = main.cc forces.cc ewald_space.cc step.cc read_paramfile.cc friedmann_solver.cc read_gadget_ic.cc nonis_friedmann.cc nonis_friedmann_voronoi.cc CIC_density.cc
#OBJ = main.o forces_cuda.o ewald_space.o step.o read_paramfile.o friedmann_solver.o read_gadget_ic.o nonis_friedmann.o nonis_friedmann_voronoi.o CIC_density.o

DEPS = global_variables.h $(VORO)/voro++.hh $(DTFE)/DTFE.h
PROG = CCLEA
ifeq ($(SINGLE_PRECISION), 1)
$(PROG): $(OBJ)
	$(CXX) $(CFLAGS) $(LDFLAGS) $(BUILD_NUMBER_LDFLAGS) $(E_INC) $(D_INC) $(E_LIB) $(D_LIB) -DUSE_SINGLE_PRECISION=$(SINGLE_PRECISION) -o $(PROG) $(OBJ) -lvoro++ -lDTFE -lCGAL -lgmp -lmpfr -lgsl -lgslcblas -lboost_program_options -lboost_system -lboost_filesystem -lboost_thread -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lcurl

%.o: %.cc
	$(CXX) $(CFLAGS) $(LDFLAGS) $(E_INC) $(D_INC) $(E_LIB) $(D_LIB) -DUSE_SINGLE_PRECISION=$(SINGLE_PRECISION) -o $@ -c $< -lDTFE -lCGAL -lgmp -lmpfr -lgsl -lgslcblas -lboost_program_options -lboost_system -lboost_filesystem -lboost_thread -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lcurl
else
$(PROG): $(OBJ)
	$(CXX) $(CFLAGS) $(LDFLAGS) $(BUILD_NUMBER_LDFLAGS) $(E_INC) $(D_INC) $(E_LIB) $(D_LIB) -o $(PROG) $(OBJ) -lvoro++ -lDTFE -lCGAL -lgmp -lmpfr -lgsl -lgslcblas -lboost_program_options -lboost_system -lboost_filesystem -lboost_thread -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lcurl

%.o: %.cc
	$(CXX) $(CFLAGS) $(LDFLAGS) $(E_INC) $(D_INC) $(E_LIB) $(D_LIB) -o $@ -c $< -lDTFE -lCGAL -lgmp -lmpfr -lgsl -lgslcblas -lboost_program_options -lboost_system -lboost_filesystem -lboost_thread -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lcurl

endif

clean:
	rm -fv $(OBJ) $(PROG)

endif
