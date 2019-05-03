CPP=g++
#CPP=clang++
CPPFLAGS= -ggdb -Wall -Wextra -Wno-sign-compare -Ofast -std=c++14

bin/PIC: main.o particle.o vose.o utils.o
	$(CPP) $^ -o $@

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $<

clean:
	rm -f *.o
