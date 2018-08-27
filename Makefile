 all: main.cpp Neighbour.cpp Track.cpp Session.cpp Conference.cpp SessionOrganizer.cpp
	  g++ -o main main.cpp Neighbour.cpp Track.cpp Session.cpp Conference.cpp SessionOrganizer.cpp

clean: 
	rm -f main
