CXX = g++ 
MPICXX = mpicxx
OPT = -O3
DEBUG_OPT = -O0 
LIBS += -lboost_program_options -lfftw3 -lopenblas -llapack -lpthread -lgfortran 
SRC = src

all : exact ehrenfest fssh pcfssh

exact : exact_flat exact_flat2 exact_flat3 exact_tully1 exact_tully12 exact_tully12cos  exact_tully22 exact_tully3 
ehrenfest : ehrenfest_flat_mpi ehrenfest_flat2_mpi ehrenfest_tully1_mpi  ehrenfest_tully12_mpi ehrenfest_tully12cos_mpi ehrenfest_tully22_mpi ehrenfest_tully3_mpi
fssh : \
	fssh_flat_mpi \
	fssh_flat2_mpi \
	fssh_tully1_mpi \
	fssh_tully12_mpi \
	fssh_tully12_avgFberry_mpi \
	fssh_flat_anal_mpi \
	fssh_flat2_anal_mpi \
	fssh_flat3_anal_mpi \
	fssh_marcus_mpi \
	fssh_helix_mpi \
	fssh_helix_ir_mpi \
	fssh_helix_ir2_mpi \
	pcfssh_yanze_mpi 

pcfssh: \
	pcfssh_flat2_anal_mpi \
	pcfssh_tully12_mpi \
	pcfssh_tully12_allinone_mpi \
	pcfssh_tully12cos_allinone_mpi \
	pcfssh_tully22_allinone_mpi 

exact_flat: $(SRC)/exact_flat.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_flat2: $(SRC)/exact_flat2.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_flat3: $(SRC)/exact_flat3.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_tully1: $(SRC)/exact_tully1.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_tully12: $(SRC)/exact_tully12.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_tully12cos: $(SRC)/exact_tully12cos.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

exact_tully22: $(SRC)/exact_tully22.cpp 
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

ehrenfest_tully12cos_mpi: $(SRC)/ehrenfest_tully12cos_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest_tully22_mpi: $(SRC)/ehrenfest_tully22_mpi.cpp 
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

fssh_tully12_avgFberry_mpi: $(SRC)/fssh_tully12_avgFberry_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_flat_anal_mpi: $(SRC)/fssh_flat_anal_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_flat2_anal_mpi: $(SRC)/fssh_flat2_anal_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_flat3_anal_mpi: $(SRC)/fssh_flat3_anal_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_helix_mpi: $(SRC)/fssh_helix_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_helix_ir_mpi: $(SRC)/fssh_helix_ir_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_helix_ir2_mpi: $(SRC)/fssh_helix_ir2_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

fssh_marcus_mpi: $(SRC)/fssh_marcus_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

pcfssh_flat2_anal_mpi: $(SRC)/pcfssh_flat2_anal_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

pcfssh_tully12_mpi: $(SRC)/pcfssh_tully12_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

pcfssh_tully12_allinone_mpi: $(SRC)/pcfssh_tully12_allinone_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

pcfssh_tully12cos_allinone_mpi: $(SRC)/pcfssh_tully12cos_allinone_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

pcfssh_tully22_allinone_mpi: $(SRC)/pcfssh_tully22_allinone_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)

pcfssh_yanze_mpi: $(SRC)/pcfssh_yanze_mpi.cpp 
	$(MPICXX) $(OPT) $< -o $@ $(LIBS)


check : $(SRC)/check_potential.cpp
	$(CXX) $(DEBUG_OPT) $< -o $@ $(LIBS)

test : $(SRC)/test.cpp
	$(CXX) $(DEBUG_OPT) $< -o $@ $(LIBS)

clean: 
	rm *_mpi 
	rm exact_*

