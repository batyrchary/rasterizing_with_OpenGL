
all:
	g++ batyr.cpp -o app -std=c++0x -lGLEW -lGL -lglut -lGLU -w -pthread
clean:
	rm -f *.o app
