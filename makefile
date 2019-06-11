CXX = g++ 
MPICXX = mpicxx
OPT = -O3 
LIBS += -lboost_program_options -lfftw3 -lopenblas -llapack -lpthread -lgfortran 
SRC = src

all : exact ehrenfest fssh 

exact : exact_flat exact_tully1 exact_tully3
ehrenfest : ehrenfest_flat_mpi ehrenfest_tully1_mpi  ehrenfest_tully3_mpi
fssh : \
	fssh_flat_mpi \
	fsshm1_flat_mpi fsshm1_tully1_mpi fsshm1_tully3_mpi \
	fsshm2_flat_mpi fsshm2_tully1_mpi fsshm2_tully3_mpi \
	fsshx_flat_mpi  fsshx_tully1_mpi fsshx_tully3_mpi \
	fsshp_flat_mpi  fsshp_tully1_mpi fsshp_tully3_mpi

exact_flat: $(SRC)/exact_flat.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_tully1: $(SRC)/exact_tully1.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_tully3: $(SRC)/exact_tully3.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_flat_mpi: $(SRC)/ehrenfest_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_tully1_mpi: $(SRC)/ehrenfest_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_tully3_mpi: $(SRC)/ehrenfest_tully3_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_flat_mpi: $(SRC)/fssh_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshm1_flat_mpi: $(SRC)/fsshm1_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshm1_tully1_mpi: $(SRC)/fsshm1_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshm1_tully3_mpi: $(SRC)/fsshm1_tully3_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshm2_flat_mpi: $(SRC)/fsshm2_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshm2_tully1_mpi: $(SRC)/fsshm2_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshm2_tully3_mpi: $(SRC)/fsshm2_tully3_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshx_flat_mpi: $(SRC)/fsshx_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshx_tully1_mpi: $(SRC)/fsshx_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshx_tully3_mpi: $(SRC)/fsshx_tully3_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshp_flat_mpi: $(SRC)/fsshp_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshp_tully1_mpi: $(SRC)/fsshp_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fsshp_tully3_mpi: $(SRC)/fsshp_tully3_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

clean: 
	rm *_mpi 
	rm exact_*
