CXX = g++ 
MPICXX = mpicxx
OPT = -O3 
LIBS += -lboost_program_options -lfftw3 -lopenblas -llapack -lpthread -lgfortran 

all : ehrenfest_flat_mpi ehrenfest_tully1_mpi fssh_flat_mpi fssh_tully1_mpi fsshx_flat_mpi fsshx_tully1_mpi

ehrenfest_flat_mpi: ehrenfest_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_tully1_mpi: ehrenfest_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_flat_mpi: fssh_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_tully1_mpi: fssh_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshx_flat_mpi: fsshx_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshx_tully1_mpi: fsshx_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)
