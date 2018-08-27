#ifndef NEIGHBOUR_SINGLE_H
#define NEIGHBOUR_SINGLE_H

#include <iostream>
using namespace std;

class Neighboursingle
{
  private:
	// The type of neighbour, parallel or time
	bool type;

	// Indices of parallel tracks
	int trackA, trackB;

	// Indices of time slots
	int timeSlotA, timeSlotB;

	// Paper Indices
	int paperIdxA, paperIdxB;

	// Profit in goodness
	double goodnessIncrement;

  public:

	/**
	 * @constructor Pass the neighbour type, the index along the type,
	 * and the two swaping indices. Also provide the exchange size and 
	 * goodness increment.
	 */
	Neighboursingle (int trkA, int trkB, int timeA, int timeB, int paperA, int paperB, double goodInc);

	/**
	 * Getter functions
	 */
	bool getType ();

	int getTrackA();
	int getTimeA();
	int getTrackB();
	int getTimeB();
	
	int getPaperIdxA ();
	int getPaperIdxB ();
	double getGoodInc ();

	void printNeighbour ();
};

#endif /* NEIGHBOUR_SINGLE_H */