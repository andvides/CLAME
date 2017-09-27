all:
	g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib clame.cpp clame_full.cpp -o clame -lsdsl -ldivsufsort -ldivsufsort64 -fopenmp
full:
	g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib clame.cpp clame_full.cpp -o clame -lsdsl -ldivsufsort -ldivsufsort64 -fopenmp
little:
	g++ -D little -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib clame.cpp clame_little.cpp -o clame -lsdsl -ldivsufsort -ldivsufsort64 -fopenmp


