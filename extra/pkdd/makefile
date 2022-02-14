all: main

CPP = g++ -std=c++11
CC = gcc
CPPFLAGS = -O3 

OBJ = shorten_Z.o main.o functions.o heuristic.o

shorten_Z.o: shorten_Z.cpp shorten_Z.hpp
	$(CPP) $(CPPFLAGS) -c $<

functions.o: functions.c functions.h
	$(CPP) $(CPPFLAGS) -c $<	

heuristic.o: heuristic.cpp heuristic.h
	$(CPP) $(CCPFLAGS) -c $<

main.o: main.cpp
	$(CPP) $(CPPFLAGS) -c -I/usr/local/include/igraph -L/usr/local/lib -ligraph  $<

main: $(OBJ)
	export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib
	$(CPP) $(CPPFLAGS) $(OBJ) -I/usr/local/include/igraph -L/usr/local/lib -ligraph -o $@ 
output: main
	./main testfile sen_stringfile k
clean:
	rm *.o main *~
