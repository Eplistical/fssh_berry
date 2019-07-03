CXX = g++ 
MPICXX = mpicxx
OPT = -O3 
LIBS += -lboost_program_options -lfftw3 -lopenblas -llapack -lpthread -lgfortran 
SRC = src

all : exact ehrenfest fssh 

exact : exact_flat exact_flat2 exact_tully1 exact_tully12 exact_tully3 
ehrenfest : ehrenfest_flat_mpi ehrenfest_flat2_mpi ehrenfest_tully1_mpi  ehrenfest_tully12_mpi ehrenfest_tully3_mpi
fssh : \
	fssh_flat_mpi \
	fssh_flat2_mpi \
	fssh_tully1_mpi \
	fssh_tully12_mpi \
	fssh_flat_anal_mpi 

exact_flat: $(SRC)/exact_flat.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_flat2: $(SRC)/exact_flat2.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_tully1: $(SRC)/exact_tully1.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_tully12: $(SRC)/exact_tully12.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_tully3: $(SRC)/exact_tully3.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_flat_mpi: $(SRC)/ehrenfest_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_flat2_mpi: $(SRC)/ehrenfest_flat2_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_tully1_mpi: $(SRC)/ehrenfest_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_tully12_mpi: $(SRC)/ehrenfest_tully12_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_tully3_mpi: $(SRC)/ehrenfest_tully3_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_flat_mpi: $(SRC)/fssh_flat_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_flat2_mpi: $(SRC)/fssh_flat2_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_tully1_mpi: $(SRC)/fssh_tully1_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_tully12_mpi: $(SRC)/fssh_tully12_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_flat_anal_mpi: $(SRC)/fssh_flat_anal_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

check : $(SRC)/check_potential.cpp
	$(CXX) $(OPT) $< -o $@ $(LIBS)

clean: 
	rm *_mpi 
	rm exact_*

