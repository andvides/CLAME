CPP=g++
FLAGS=-O3 -std=c++11 -fopenmp 
LFLAGS=-lsdsl -ldivsufsort -ldivsufsort64
FLAGSD=$(FLAGS) -D debug  -DNDEBUG
ILIB=-I ~/include
LLIB=-L ~/lib
OBJ = genFm9
SOURCES= genFm9.cpp genFm9_lite.cpp
all:
	$(CPP) $(FLAGS) $(ILIB) $(LLIB) $(SOURCES) -o $(OBJ)  $(LFLAGS)

debug:
	$(CPP) $(FLAGSD) $(ILIB) $(LLIB) $(SOURCES) -o $(OBJ) $(LFLAGS) 

clean:
	rm $(OBJ)









