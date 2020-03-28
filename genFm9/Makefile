CPP=g++
CFLAGS=-O3 -std=c++11 -fopenmp 
LFLAGS=-lsdsl -ldivsufsort -ldivsufsort64
CFLAGSD=$(CFLAGS) -D debug  -DNDEBUG
ILIB=-I ~/include
LLIB=-L ~/lib
OBJ = genFm9
SOURCES= genFm9.cpp genFm9_lite.cpp
all:
	$(CPP) $(CFLAGS) $(ILIB) $(LLIB) $(SOURCES) -o $(OBJ)  $(LFLAGS)

debug:
	$(CPP) $(CFLAGSD) -I $(ILIB) -L $(LLIB) $(SOURCES) -o $(OBJ) $(LFLAGS) 

clean:
	rm $(OBJ)






