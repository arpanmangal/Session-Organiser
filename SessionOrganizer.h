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

    double neighbourRowSelectionProb;
    int numberNeighbours;
    double neighbourGeoProb;

    Conference *conference;

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
    // vector<Neighbour> getNeighbours_nc2 ();

    /**
     * Make and return a neighbour with given parameters
     */
    Neighbour getNeighbour (bool neighbourType, int trkA, int trkB, int timeA, int timeB, int exSize);
    // Neighbour getNeighbour_nc2 ();

    /**
     * Get the change in goodness w.r.t. A, when swapping exSize elements with B
     */
    double sessionExchangeGoodness(int trkA, int trkB, int timeA, int timeB, int exSize);
    // double timeSessionExchangeGoodness (int )
    /**
     * Change to go the specified neighbour
     */
    void gotoNeighbour (Neighbour ngh);

    /**
     * Do the local search on this state
     */
    // void localSearch ();


  public:
    SessionOrganizer();
    SessionOrganizer(string filename);
    void initializeConference();
    void localSearch ();

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
