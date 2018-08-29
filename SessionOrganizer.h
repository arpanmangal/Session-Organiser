/* 
 * File:   SessionOrganizer.h
 * Author: Kapil Thakkar
 *
 */

#ifndef SESSIONORGANIZER_H
#define SESSIONORGANIZER_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "Conference.h"
#include "Track.h"
#include "Session.h"
#include "Neighbour.h"
#include "Neighboursingle.h"

using namespace std;

/**
 * SessionOrganizer reads in a similarity matrix of papers, and organizes them
 * into sessions and tracks.
 * 
 * @author Kapil Thakkar
 *
 */
class SessionOrganizer
{
  private:
    double **distanceMatrix;

    int parallelTracks;  // p
    int papersInSession; // k
    int sessionsInTrack; // t
    int totalPapers;

    Conference *conference;
    Conference *bestConference; // The Maximum found by our algorithm
    double maxGoodness; // The maximum goodness found by our algorithm

    double neighbourRowSelectionProb;
    int numberNeighbours;
    double neighbourGeoProb;

    double processingTimeInMinutes;
    double tradeoffCoefficient; // the tradeoff coefficient

    /**
     * Initialize the conference according to the sorting rule 
     */
    // void initializeConference();

    /** 
     * Get the neighbours of this state
     */
    vector<Neighbour> getNeighbours ();
    vector<Neighboursingle> getNeighbours_nc2 (bool& typ, int& prob);

    /**
     * Make and return a neighbour with given parameters
     */
    Neighbour getNeighbour (bool neighbourType, int trkA, int trkB, int timeA, int timeB, int exSize);
    Neighboursingle getNeighbour_nc2 (int trkA, int trkB, int timeA, int timeB, int paperIdxA, int paperIdxB);

    /**
     * Get the change in goodness w.r.t. A, when swapping exSize elements with B
     */
    double sessionExchangeGoodness(int trkA, int trkB, int timeA, int timeB, int exSize);
    double sessionExchangeGoodness_nc2(int trkA, int trkB, int timeA, int timeB, int paperIdxA, int paperIdxB);
    // double timeSessionExchangeGoodness (int )
    /**
     * Change to go the specified neighbour
     */
    void gotoNeighbour (Neighbour ngh);
    void gotoNeighbour_nc2 (Neighboursingle ngh);


    /**
     * Do the local search on this state
     */
    // void localSearch ();


  public:
    // SessionOrganizer();
    SessionOrganizer(string filename);
    void initializeConference();
    void initializeConferenceFullSort();
    void initializeConferenceRandomly();
    void localSearch ();
    void localSearch_nc2();

    void updateMaximum (double newMaxGoodness);

    /**
     * Read in the number of parallel tracks, papers in session, sessions
     * in a track, and the similarity matrix from the specified filename.
     * @param filename is the name of the file containing the matrix.
     * @return the similarity matrix.
     */
    void readInInputFile(string filename);

    /**
     * Organize the papers according to some algorithm.
     */
    void organizePapers();

    /**
     * Get the distance matrix.
     * @return the distance matrix.
     */
    double **getDistanceMatrix();

    /**
     * Score the organization.
     * @return the score.
     */
    double scoreOrganization();

    void printSessionOrganiser(char *);
};

#endif /* SESSIONORGANIZER_H */
