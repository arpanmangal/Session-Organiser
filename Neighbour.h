#ifndef NEIGHBOUR_H
#define NEIGHBOUR_H

#include <iostream>
using namespace std;

class Neighbour
{
  private:
	// The type of neighbour, parallel or time
	bool type;

	// Indices of parallel tracks
	int trackA, trackB;

	// Indices of time slots
	int timeSlotA, timeSlotB;

	// Exchange Size
	int exchangeSize;

	// Profit in goodness
	double goodnessIncrement;

  public:

	/**
	 * @constructor Pass the neighbour type, the index along the type,
	 * and the two swaping indices. Also provide the exchange size and 
	 * goodness increment.
	 */
	Neighbour (bool neighbourType, int trkA, int trkB, int timeA, int timeB, int exSize, double goodInc);

	/**
	 * Getter functions
	 */
	bool getType ();

	int getTrackA();
	int getTimeA();
	int getTrackB();
	int getTimeB();
	
	int getExSize ();
	double getGoodInc ();
};

#endif /* NEIGHBOUR_H */