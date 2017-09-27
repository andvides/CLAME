all:
        g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib clame.cpp clame_lite.cpp -o clame -lsdsl -ldivsufsort -ldivsufsort64 -fopenmp
debug:
        g++ -D debug -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib clame.cpp clame_lite.cpp -o clame -lsdsl -ldivsufsort -ldivsufsort64 -fopenmp


