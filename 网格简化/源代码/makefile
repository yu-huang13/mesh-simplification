all: main.exe clear

heap.o: heap.cpp heap.h
	g++ -c heap.cpp
	
point.o: point.cpp point.h
	g++ -c point.cpp

triangle.o: triangle.cpp triangle.h
	g++ -c triangle.cpp

pair.o: pair.cpp pair.h
	g++ -c pair.cpp
	
vector3.o: vector3.cpp vector3.h
	g++ -c vector3.cpp

main.exe: main.cpp  heap.o point.o triangle.o pair.o vector3.o 
	g++ -o main.exe -std=c++11 -O2 main.cpp heap.o point.o triangle.o pair.o vector3.o 

clear: 
	del *.bak *.o
