CXXFLAGS=-std=c++11

nbsim: r3vec.o particle.o engine.o universe.o minifb/src/x11/X11MiniFB.o
	g++ -std=c++11 -o nbsim engine.o particle.o\
		 r3vec.o universe.o minifb/src/x11/X11MiniFB.o  -lX11 -I.
