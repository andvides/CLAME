CPP=g++
CFLAGS=-O3 -std=c++11 -fopenmp 
LFLAGS=-lsdsl -ldivsufsort -ldivsufsort64
CFLAGSD=$(CFLAGS) -D debug  -DNDEBUG
ILIB=~/include
LLIB=~/lib
OBJ = clame 
SOURCES= clame.cpp clame_lite.cpp
all:
	$(CPP) $(CFLAGS) -I $(CPPFLAGS) -L $(LLIB) $(SOURCES) -o $(OBJ)  $(LFLAGS)

debug:
	$(CPP) $(CFLAGSD) -I $(CPPFLAGS=) -L $(LLIB) $(SOURCES) -o $(OBJ) $(LFLAGS) 

clean:
	rm $(OBJ)
