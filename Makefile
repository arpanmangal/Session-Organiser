 all: main.cpp Neighboursingle.cpp Neighbour.cpp Track.cpp Session.cpp Conference.cpp SessionOrganizer.cpp
	  g++ --std=c++11 -o main main.cpp Neighboursingle.cpp Neighbour.cpp Track.cpp Session.cpp Conference.cpp SessionOrganizer.cpp

clean: 
	rm -f main
