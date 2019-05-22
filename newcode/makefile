CXX = g++ 
MPICXX = mpicxx
OPT = -O3 
LIBS += -lboost_program_options -lfftw3 -lopenblas -llapack -lpthread -lgfortran 

all : exact ehrenfest fssh

exact : exact_flat exact_tully1 
ehrenfest : ehrenfest_flat_mpi ehrenfest_tully1_mpi 
fssh : \
	fsshm1_flat_mpi fsshm1_tully1_mpi \
	fsshm2_flat_mpi fsshm2_tully1_mpi \
	fsshx_flat_mpi  fsshx_tully1_mpi \
	fsshp_flat_mpi  fsshp_tully1_mpi

exact_flat: exact_flat.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_tully1: exact_tully1.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_flat_mpi: ehrenfest_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_tully1_mpi: ehrenfest_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshm1_flat_mpi: fsshm1_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshm1_tully1_mpi: fsshm1_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshm2_flat_mpi: fsshm2_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshm2_tully1_mpi: fsshm2_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshx_flat_mpi: fsshx_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshx_tully1_mpi: fsshx_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshp_flat_mpi: fsshp_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshp_tully1_mpi: fsshp_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)
