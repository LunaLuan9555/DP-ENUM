Objs:=main.o tools.o reduction.o psucc.o gsa.o cellenum.o simulator.o dpenum.o
# Srcs:=main_time.cpp cellenum.cpp reduction.cpp 
Defs:= tools.h reduction.h psucc.h gsa.h cellenum.h simulator.h dpenum.h
main: $(Objs)
	g++ -o ../DPENUM $(Objs) -lm -lntl -lgmp -lmpfr -lfplll
main.o: main.cpp $(Defs)
	g++ -c -O3 -std=c++11 main.cpp
tools.o: tools.cpp tools.h
	g++ -c -O3 -std=c++11 tools.cpp
reduction.o: reduction.cpp tools.h reduction.h
	g++ -c -O3 -std=c++11 reduction.cpp
psucc.o: psucc.cpp psucc.h
	g++ -c -O3 -std=c++11 psucc.cpp
cellenum.o: cellenum.cpp cellenum.h reduction.h
	g++ -c -O3 -std=c++11 cellenum.cpp
gsa.o: gsa.cpp $(Defs)
	g++ -c -O3 -std=c++11 gsa.cpp
simulator.o: simulator.cpp $(Defs)
	g++ -c -O3 -std=c++11 simulator.cpp
dpenum.o: dpenum.cpp $(Defs)
	g++ -c -O3 -std=c++11 dpenum.cpp
clean:
	rm $(Objs)