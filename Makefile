FLAGS = -Wall -std=c++11

all: tubeplanner

tubeplanner: tubeplanner.o
	g++ $(FLAGS) -o $@ $^

%.o: %.cc tubeplanner.h
	g++ $(FLAGS) -c $<

clean:
	rm -f tubeplanner *.o
