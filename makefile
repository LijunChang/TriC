CFLAGS=-c -I /import/ravel/1/ljchang/sigmod09/sparsehash/include/

all: trilist

trilist: .obj/main.o .obj/Graph.o .obj/HashSet.o .obj/Timer.o .obj/Utility.o
	g++ .obj/main.o .obj/Graph.o .obj/HashSet.o .obj/Timer.o .obj/Utility.o -o trilist
	rm .obj/*

.obj/main.o: main.cpp
	g++ -O3 ${CFLAGS} -o .obj/main.o main.cpp

.obj/Graph.o: Graph.cpp
	g++ -O3 ${CFLAGS} -o .obj/Graph.o Graph.cpp

.obj/HashSet.o: HashSet.cpp
	g++ -O3 ${CFLAGS} -o .obj/HashSet.o HashSet.cpp

.obj/Timer.o: Timer.cpp
	g++ -O3 ${CFLAGS} -o .obj/Timer.o Timer.cpp
	
.obj/Utility.o: Utility.cpp
	g++ -O3 ${CFLAGS} -o .obj/Utility.o Utility.cpp

clean:
	rm -rf *o .obj/
	mkdir .obj
