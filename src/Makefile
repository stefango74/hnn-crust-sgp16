all: TestReconstruct2D

TestReconstruct2D: TestReconstruct2D.o Reconstruct2D.o circle.o
	g++ -g -o TestReconstruct2D TestReconstruct2D.o Reconstruct2D.o circle.o -lglut -lpng -lann -lGL -lGLU -lGLEW

TestReconstruct2D.o: TestReconstruct2D.cpp
	g++ -std=c++11 -g -c -I. TestReconstruct2D.cpp

Reconstruct2D.o: Reconstruct2D.cpp
	g++ -std=c++11 -g -c -I. Reconstruct2D.cpp

circle.o: circle.cpp
	g++ -std=c++11 -g -c -I. circle.cpp

clean:
	rm -rf *o TestReconstruct2D
