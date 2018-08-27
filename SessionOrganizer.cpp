/* 
 * File:   SessionOrganizer.cpp
 * Author: Kapil Thakkar
 * 
 */

#include "SessionOrganizer.h"
#include "Util.h"
#include <algorithm>
#include <ctime>
#include <math.h>

SessionOrganizer::SessionOrganizer()
{
    parallelTracks = 0;
    papersInSession = 0;
    sessionsInTrack = 0;
    processingTimeInMinutes = 0;
    tradeoffCoefficient = 1.0;
}

SessionOrganizer::SessionOrganizer(string filename)
{
    readInInputFile(filename);
    conference = new Conference(parallelTracks, sessionsInTrack, papersInSession);

    // Initialising various probabilities
    neighbourRowSelectionProb = 0.7;
    numberNeighbours = 10;
    neighbourGeoProb = 1 / sqrt(papersInSession);
}

void SessionOrganizer::organizePapers()
{
    int paperCounter = 0;
    for (int i = 0; i < conference->getSessionsInTrack(); i++)
    {
        for (int j = 0; j < conference->getParallelTracks(); j++)
        {
            for (int k = 0; k < conference->getPapersInSession(); k++)
            {
                conference->setPaper(j, i, k, paperCounter);
                paperCounter++;
            }
        }
    }
}

void SessionOrganizer::readInInputFile(string filename)
{
    vector<string> lines;
    string line;
    ifstream myfile(filename.c_str());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            //cout<<"Line read:"<<line<<endl;
            lines.push_back(line);
        }
        myfile.close();
    }
    else
    {
        cout << "Unable to open input file";
        exit(0);
    }

    if (6 > lines.size())
    {
        cout << "Not enough information given, check format of input file";
        exit(0);
    }

    processingTimeInMinutes = atof(lines[0].c_str());
    papersInSession = atoi(lines[1].c_str());     // k
    parallelTracks = atoi(lines[2].c_str());      // p
    sessionsInTrack = atoi(lines[3].c_str());     // t
    tradeoffCoefficient = atof(lines[4].c_str()); // C

    int n = lines.size() - 5;
    double **tempDistanceMatrix = new double *[n];
    for (int i = 0; i < n; ++i)
    {
        tempDistanceMatrix[i] = new double[n];
    }

    for (int i = 0; i < n; i++)
    {
        string tempLine = lines[i + 5];
        string elements[n];
        splitString(tempLine, " ", elements, n);

        for (int j = 0; j < n; j++)
        {
            tempDistanceMatrix[i][j] = atof(elements[j].c_str());
        }
    }
    distanceMatrix = tempDistanceMatrix;

    int numberOfPapers = n;
    int slots = parallelTracks * papersInSession * sessionsInTrack;
    if (slots != numberOfPapers)
    {
        cout << "More papers than slots available! slots:" << slots << " num papers:" << numberOfPapers << endl;
        exit(0);
    }
}

double **SessionOrganizer::getDistanceMatrix()
{
    return distanceMatrix;
}

void SessionOrganizer::printSessionOrganiser(char *filename)
{
    conference->printConference(filename);
}

double SessionOrganizer::scoreOrganization()
{
    // Sum of pairwise similarities per session.
    double score1 = 0.0;
    for (int i = 0; i < conference->getParallelTracks(); i++)
    {
        Track tmpTrack = conference->getTrack(i);
        for (int j = 0; j < tmpTrack.getNumberOfSessions(); j++)
        {
            Session tmpSession = tmpTrack.getSession(j);
            for (int k = 0; k < tmpSession.getNumberOfPapers(); k++)
            {
                int index1 = tmpSession.getPaper(k);
                for (int l = k + 1; l < tmpSession.getNumberOfPapers(); l++)
                {
                    int index2 = tmpSession.getPaper(l);
                    score1 += 1 - distanceMatrix[index1][index2];
                }
            }
        }
    }

    // Sum of distances for competing papers.
    double score2 = 0.0;
    for (int i = 0; i < conference->getParallelTracks(); i++)
    {
        Track tmpTrack1 = conference->getTrack(i);
        for (int j = 0; j < tmpTrack1.getNumberOfSessions(); j++)
        {
            Session tmpSession1 = tmpTrack1.getSession(j);
            for (int k = 0; k < tmpSession1.getNumberOfPapers(); k++)
            {
                int index1 = tmpSession1.getPaper(k);

                // Get competing papers.
                for (int l = i + 1; l < conference->getParallelTracks(); l++)
                {
                    Track tmpTrack2 = conference->getTrack(l);
                    Session tmpSession2 = tmpTrack2.getSession(j);
                    for (int m = 0; m < tmpSession2.getNumberOfPapers(); m++)
                    {
                        int index2 = tmpSession2.getPaper(m);
                        score2 += distanceMatrix[index1][index2];
                    }
                }
            }
        }
    }
    double score = score1 + tradeoffCoefficient * score2;
    return score;
}

void SessionOrganizer::initializeConference()
{
    // choose a random row
    totalPapers = parallelTracks * papersInSession * sessionsInTrack;
    // srand(time(NULL));
    int row = rand() / ((RAND_MAX + 1u) / totalPapers);

    if (row < 0 || row >= totalPapers)
    {
        cout << "Index out of range - SessionOrganizer::initializeConference" << endl;
        exit(0);
    }

    vector<pair<double, int>> arr;
    for (int j = 0; j < totalPapers; j++)
    {
        arr.push_back(make_pair(distanceMatrix[row][j], j));
    }
    sort(arr.begin(), arr.end());

    int paperCounter = 0;
    for (int j = 0; j < parallelTracks; j++)
    {
        for (int i = 0; i < sessionsInTrack; i++)
        {
            for (int k = 0; k < papersInSession; k++)
            {
                conference->setPaper(j, i, k, arr[paperCounter].second);
                paperCounter++;
            }
        }
    }
}

vector<Neighbour> SessionOrganizer::getNeighbours()
{
    // srand(time(NULL));
    int random = rand();
    bool rowNeigh = (double(random) / ((RAND_MAX + 1u)) < neighbourRowSelectionProb) ? true : false;
    // cout << "rowNeigh: " << rowNeigh << " | " << random << endl;

    // srand(time(NULL));
    if (rowNeigh)
    {
        // Select neighbours from a random row
        int row = rand() / ((RAND_MAX + 1u) / sessionsInTrack);

        // Generate `numNeighbours` neighbours
        vector<Neighbour> neighbours;
        int max_papers = min(numberNeighbours, parallelTracks * (parallelTracks - 1) / 2);
        for (int nh = 0; nh < max_papers; nh++)
        {
            // choose two track numbers at random
            int trackA = 0;
            int trackB = 0;

            while (trackA == trackB)
            {
                // find a new pair at random
                trackA = rand() / ((RAND_MAX + 1u) / parallelTracks);
                trackB = rand() / ((RAND_MAX + 1u) / parallelTracks);
            }

            // choose the exchange size with geometric distribution
            int exSize = 1;
            while (exSize < papersInSession)
            {
                // toss a coin
                bool shouldStop = (double(rand()) / ((RAND_MAX + 1u)) < neighbourGeoProb) ? false : true;
                if (shouldStop)
                {
                    // got a tails
                    break;
                }
                else
                {
                    // got heads
                    exSize++;
                }
            }

            // create a neighbour with the obtained parameters
            Neighbour ngbour = getNeighbour(true, trackA, trackB, row, row, exSize);

            neighbours.push_back(ngbour);
        }

        return neighbours;
    }
    else
    {
        // Select neighbours from a random column
        int col = rand() / ((RAND_MAX + 1u) / parallelTracks);

        // Generate `numNeighbours` neighbours
        vector<Neighbour> neighbours;
        int max_papers = min(numberNeighbours, sessionsInTrack * (sessionsInTrack - 1) / 2);
        for (int nh = 0; nh < numberNeighbours; nh++)
        {
            // choose two time slots at random
            int timeA = 0;
            int timeB = 0;

            while (timeA == timeB)
            {
                // find a new pair at random
                timeA = rand() / ((RAND_MAX + 1u) / sessionsInTrack);
                timeB = rand() / ((RAND_MAX + 1u) / sessionsInTrack);
            }

            // choose the exchange size with geometric distribution
            int exSize = 1;
            while (exSize < papersInSession)
            {
                // toss a coin
                bool shouldStop = (double(rand()) / ((RAND_MAX + 1u)) < neighbourGeoProb) ? false : true;
                if (shouldStop)
                {
                    // got a tails
                    break;
                }
                else
                {
                    // got heads
                    exSize++;
                }
            }

            // create a neighbour with the obtained parameters
            Neighbour ngbour = getNeighbour(true, col, col, timeA, timeB, exSize);

            neighbours.push_back(ngbour);
        }

        return neighbours;
    }
}

Neighbour SessionOrganizer::getNeighbour(bool neighbourType, int trkA, int trkB, int timeA, int timeB, int exSize)
{
    // Compute the goodness increment of the neighbour, and return a neighbour object
    double goodnessChange = 0;
    goodnessChange += sessionExchangeGoodness(trkA, trkB, timeA, timeB, exSize);
    goodnessChange += sessionExchangeGoodness(trkB, trkA, timeB, timeA, exSize);

    // make and return a new neighbour
    Neighbour ng(neighbourType, trkA, trkB, timeA, timeB, exSize, goodnessChange);
    return ng;
}

double SessionOrganizer::sessionExchangeGoodness(int trkA, int trkB, int timeA, int timeB, int exSize)
{
    double delta = 0;
    int paperA, paperB;
    for (int p = 0; p < exSize; p++)
    {
        paperA = conference->getPaper(trkA, timeA, p);

        // add the goodness change due to leaving papers of same sessions
        for (int j = exSize; j < papersInSession; j++)
        {
            paperB = conference->getPaper(trkA, timeA, j);
            delta += (tradeoffCoefficient + 1) * distanceMatrix[paperA][paperB] - 1;
        }

        // add goodness change due to embracing papers of different sessions
        for (int j = exSize; j < papersInSession; j++)
        {
            paperB = conference->getPaper(trkB, timeB, j);
            delta -= (tradeoffCoefficient + 1) * distanceMatrix[paperA][paperB] - 1;
        }
    }

    // return the change
    return delta;
}

void SessionOrganizer::gotoNeighbour(Neighbour ngh)
{
    // Convert the current conference state to that represented by the given neighbour
    int exSize = ngh.getExSize();
    for (int p = 0; p < exSize; p++)
    {
        // Exchange the pth paper in sessions A and B
        int paperA = conference->getPaper(ngh.getTrackA(), ngh.getTimeA(), p);
        int paperB = conference->getPaper(ngh.getTrackB(), ngh.getTimeB(), p);

        conference->setPaper(ngh.getTrackA(), ngh.getTimeA(), p, paperB);
        conference->setPaper(ngh.getTrackB(), ngh.getTimeB(), p, paperA);
    }
}

void SessionOrganizer::localSearch()
{
    cout << sessionsInTrack << " | " << parallelTracks << " | " << papersInSession << " | " << totalPapers << endl;
    cout << neighbourRowSelectionProb << " | " << numberNeighbours << " | " << neighbourGeoProb << endl;
    int iter = 0;
    double score = scoreOrganization();
    cout << "score:" << score << endl;
    while (iter++ < 50)
    {
        vector<Neighbour> neighbours = getNeighbours();
        // cout << iter << " " << neighbours.size() << endl;
        if (neighbours.size() < 1)
        {
            // no neighbours
            continue;
        }

               
        for (int nh = 0; nh < neighbours.size(); nh++)
        {
            neighbours.at(nh).printNeighbour();
        }
        cout << endl; 
        

        int max_nh_idx = 0;
        for (int nh = 1; nh < neighbours.size(); nh++)
        {
            if (neighbours.at(nh).getGoodInc() > neighbours.at(max_nh_idx).getGoodInc())
            {
                max_nh_idx = nh;
            }
        }

        if (neighbours.at(max_nh_idx).getGoodInc() <= 0)
        {
            // on an local optima
            break;
        }
        else
        {
            // goto neighbour
            gotoNeighbour(neighbours.at(max_nh_idx));
        }

        double score = scoreOrganization();
        cout << "score:" << score << endl;
    }
    cout << "Took " << iter << " steps" << endl;
}