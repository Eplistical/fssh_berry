CXX = g++ 
OPT = -O3 
LIBS = -lboost_program_options -lfftw3 -lopenblas -llapack -lpthread -lgfortran 

all : exact_2d fssh_2d fssh_2d_showprob fssh_2d_a2d2 fssh_2d_a2d2_rechop ehrenfest 

exact_2d: exact_2d.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)
	
fssh_2d: fssh_2d.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

fssh_2d_showprob: fssh_2d_showprob.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

fssh_2d_a2d2: fssh_2d_a2d2.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

fssh_2d_a2d2_rechop: fssh_2d_a2d2_rechop.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)

ehrenfest: ehrenfest.cpp 
	$(CXX) $(OPT) $< -o $@ $(LIBS)
