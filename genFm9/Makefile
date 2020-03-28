CPP=g++

LFLAGS=-lsdsl -ldivsufsort -ldivsufsort64
CFLAGSD=$(CFLAGS) -D debug  -DNDEBUG
ILIB=~/include
LLIB=~/lib
OBJ = genFm9
SOURCES= genFm9.cpp genFm9_lite.cpp
all:
	$(CPP) $(CFLAGS) $(CFLAGS) $(LDFLAGS) $(SOURCES) -o $(OBJ)  $(LFLAGS)

debug:
	$(CPP) $(CFLAGSD) -I $(CPPFLAGS) -L $(LLIB2) $(SOURCES) -o $(OBJ) $(LFLAGS) 

clean:
	rm $(OBJ)



