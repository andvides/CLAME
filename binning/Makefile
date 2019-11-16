CCP=g++
CFLAGS=-O3 -std=c++11 -fopenmp 
LFLAGS=-lsdsl -ldivsufsort -ldivsufsort64
CFLAGSD=$(CFLAGS) -D debug  -DNDEBUG
ILIB=
LLIB=
OBJ = binning
SOURCES= binning.cpp
all:
	$(CCP) $(CFLAGS) $(SOURCES) -o $(OBJ)  

debug:
	$(CCP) $(CFLAGSD) -I $(ILIB) -L $(LLIB) $(SOURCES) -o $(OBJ) $(LFLAGS) 

clean:
	rm $(OBJ)

