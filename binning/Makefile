CPP= g++
CFLAGS= -O3 -std=c++11 -fopenmp 
LFLAGS= -lsdsl -ldivsufsort -ldivsufsort64
CFLAGSD=$(CFLAGS) -D debug  -DNDEBUG
ILIB=
LLIB=
OBJ = binning
SOURCES= binning.cpp
all:
	$(CPP) $(CFLAGS) $(SOURCES) -o $(OBJ)  

debug:
	$(CPP) $(CFLAGSD) -I $(ILIB) -L $(LLIB) $(SOURCES) -o $(OBJ) $(LFLAGS) 

clean:
	rm $(OBJ)

