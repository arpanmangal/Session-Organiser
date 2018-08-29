/* 
 * File:   SessionOrganizer.cpp
 * Author: Kapil Thakkar
 * 
 */

//Lines 750 and 817.

#include "SessionOrganizer.h"
#include "Util.h"
#include <algorithm>
#include <ctime>
#include <math.h>

// SessionOrganizer::SessionOrganizer()
// {
//     parallelTracks = 0;
//     papersInSession = 0;
//     sessionsInTrack = 0;
//     processingTimeInMinutes = 0;
//     tradeoffCoefficient = 1.0;
// }

SessionOrganizer::SessionOrganizer(string filename)
{
    // Initialize time
    organizerStartTime = time(NULL) - 1; // current time - 1 to account for some delays

    readInInputFile(filename);
    conference = new Conference(parallelTracks, sessionsInTrack, papersInSession);
    bestConference = new Conference(parallelTracks, sessionsInTrack, papersInSession);
    maxGoodness = -1;

    processingTimeInSeconds = processingTimeInMinutes * 60;

    // Initialising various probabilities
    neighbourRowSelectionProb = 0.5;
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
    /* Print the best conference in the file */
    bestConference->printConference(filename);
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

void SessionOrganizer::initializeConferenceFullSort()
{
    totalPapers = parallelTracks * papersInSession * sessionsInTrack;
    vector<pair<double, int>> arr;
    for (int j = 0; j < totalPapers; j++)
    {
        arr.push_back(make_pair(0, j));
    }

    for (int p = 0; p < parallelTracks; p++)
    {
        int start = p * sessionsInTrack * papersInSession;

        // Choose a random no. between start and n-1, interval size = n - start
        int pivot = start + rand() / ((RAND_MAX + 1u) / (totalPapers - start));
        if (pivot < start || pivot >= totalPapers)
        {
            cout << "Pivot out of range - SessionOrganizer::initializeConference" << endl;
            pivot = start;
        }

        // Fill in the distances according to the pivot
        for (int j = start; j < totalPapers; j++)
        {
            int row = arr[pivot].second;
            int col = arr[j].second;
            arr[j].first = distanceMatrix[row][col];
        }

        // Sort the remaining distances
        sort(arr.begin() + start, arr.end());

        // fill the next set
        for (int t = 0; t < sessionsInTrack; t++)
        {
            for (int k = 0; k < papersInSession; k++)
            {
                conference->setPaper(p, t, k, arr[start + t * papersInSession + k].second);
            }
        }

        // cout << start << " | " << pivot << endl;
        // for (int j = 0; j < totalPapers; j++) {
        //     cout << arr[j].second << "|" << arr[j].first << ", ";
        // }
        // cout << endl;
    }
}

void SessionOrganizer::initializeConferenceRandomly()
{
    totalPapers = parallelTracks * papersInSession * sessionsInTrack;
    vector<int> arr;
    for (int j = 0; j < totalPapers; j++)
    {
        arr.push_back(j);
    }
    random_shuffle(arr.begin(), arr.end());

    int paperCounter = 0;
    for (int j = 0; j < parallelTracks; j++)
    {
        for (int i = 0; i < sessionsInTrack; i++)
        {
            for (int k = 0; k < papersInSession; k++)
            {
                conference->setPaper(j, i, k, arr[paperCounter]);
                paperCounter++;
            }
        }
    }
}

vector<Neighboursingle> SessionOrganizer::getNeighbours()
{
    vector<Neighboursingle> neighbours;

    int par_m = min(((parallelTracks)/4 +1),parallelTracks); 
    int ses_m = min(sessionsInTrack, (sessionsInTrack/4+1));
    int pap_nu = min(papersInSession, (papersInSession/4 +1));

    int arr1[par_m];
    int arr2[ses_m];
    int arr3[pap_nu];

    for(int i = 0; i< par_m;i++)
    {
        arr1[i] = rand() % parallelTracks;
    }

    for(int i = 0; i< ses_m ;i++)
    {
        arr2[i] = rand() % sessionsInTrack;
    }
    
    for(int i = 0; i< pap_nu;i++)
    {
        arr3[i] = rand() % papersInSession;
    }
    double goodnessChange = 0;

    for (int i = 0; i < par_m; i++)
    {
        for (int j = 0; j < ses_m; j++)
        {
            for (int k = 0; k < pap_nu; k++)
            {
                int i1, j1, k1;

                for (i1 = (i+1) ; i1 < par_m; i1++)
                {
                    for (j1 = (j+1); j1 < ses_m; j1++)
                    {
                        for (k1 = (k+1); k1 < pap_nu; k1++)
                        {
                            //int i1,arr[j1],arr[k1];
                            goodnessChange = 0;
                            goodnessChange += sessionExchangeGoodness_nc2(arr1[i],arr1[i1],arr2[j],arr2[j1],arr3[k],arr3[k1]);
                            goodnessChange += sessionExchangeGoodness_nc2(arr1[i1],arr1[i],arr2[j1],arr2[j],arr3[k1],arr3[k]);

                            // make and return a new neighbour
                            //Neighboursingle ng(trkA, trkB, timeA, timeB, paperIdxA, paperIdxB, goodnessChange);

                            if (goodnessChange > -2.0)
                            {
                                Neighboursingle ng(arr1[i],arr1[i1],arr2[j],arr2[j1],arr3[k],arr3[k1],goodnessChange);
                                neighbours.push_back(ng);
                            }
                            
                        }
                    }
                }
                
            }
        }
    }

    return neighbours;
}

double SessionOrganizer::sessionExchangeGoodness(int trkA, int trkB, int timeA, int timeB, int exSize)
{
    double delta = 0;
    int paperA, paperB;
    if (timeA == timeB)
    {
        // Efficient approach for row exchange type
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
    }
    else
    {
        for (int p = 0; p < exSize; p++)
        {
            paperA = conference->getPaper(trkA, timeA, p);

            for (int track = 0; track < parallelTracks; track++)
            {
                // subtract distances and similarities for same time slot
                if (track == trkA)
                {
                    // subtract similarities
                    for (int j = exSize; j < papersInSession; j++)
                    {
                        paperB = conference->getPaper(track, timeA, j);
                        delta -= (1 - distanceMatrix[paperA][paperB]);
                    }
                }
                else
                {
                    // subtract differences
                    for (int j = 0; j < papersInSession; j++)
                    {
                        paperB = conference->getPaper(track, timeA, j);
                        delta -= tradeoffCoefficient * distanceMatrix[paperA][paperB];
                    }
                }

                // add distances and similarities for new time slot
                if (track == trkB)
                {
                    // subtract similarities
                    for (int j = exSize; j < papersInSession; j++)
                    {
                        paperB = conference->getPaper(track, timeB, j);
                        delta += (1 - distanceMatrix[paperA][paperB]);
                    }
                }
                else
                {
                    // subtract differences
                    for (int j = 0; j < papersInSession; j++)
                    {
                        paperB = conference->getPaper(track, timeB, j);
                        delta += tradeoffCoefficient * distanceMatrix[paperA][paperB];
                    }
                }
            }
        }
    }
    return delta;

    // double delta = 0;

    // return the change
    return delta;
}

vector<Neighboursingle> SessionOrganizer::getNeighbours_nc2(bool &close, int &prob)
{
    vector<Neighboursingle> neighbours;
    // Iterate over all pairs of sessions, and pick all pairs of papers
    for (int trackA = 0; trackA < parallelTracks; trackA++)
    {
        //bool TrueFalse = true;
        bool TrueFalse = (rand() % 100) < (100 - prob);
        int timeA, limit;
        if (TrueFalse)
        {
            timeA = (rand() % sessionsInTrack);
            limit = timeA;
            close = false;
        }

        else
        {
            timeA = 0;
            limit = sessionsInTrack - 1;
            close = true;
        }

        for (timeA; timeA <= limit; timeA++)
        {

            for (int trackB = trackA; trackB < parallelTracks; trackB++)
            {
                int timeB, paperA, paperB;
                double goodnessChange = 0;
                //Neighboursingle ng;

                if (trackA == trackB)
                    timeB = timeA + 1;
                else
                    timeB = 0;

                for (timeB; timeB < sessionsInTrack; timeB++)
                {
                    // pick all pairs of papers
                    for (paperA = 0; paperA < papersInSession; paperA++)
                    {
                        for (paperB = 0; paperB < papersInSession; paperB++)
                        {
                            goodnessChange = 0;
                            goodnessChange += sessionExchangeGoodness_nc2(trackA, trackB, timeA, timeB, paperA, paperB);
                            goodnessChange += sessionExchangeGoodness_nc2(trackB, trackA, timeB, timeA, paperB, paperA);

                            // make and return a new neighbour
                            //Neighboursingle ng(trkA, trkB, timeA, timeB, paperIdxA, paperIdxB, goodnessChange);

                            if (goodnessChange > -8.0)
                            {
                                Neighboursingle ng(trackA, trackB, timeA, timeB, paperA, paperB, goodnessChange);
                                neighbours.push_back(ng);
                            }
                            //cout<< trackA <<"|"<<timeA <<"|"<<trackB <<"|"<<timeB <<"|"<<endl;
                        }
                    }
                }
            }
        }
    }

    return neighbours;
}

Neighboursingle SessionOrganizer::getNeighbour_nc2(int trkA, int trkB, int timeA, int timeB, int paperIdxA, int paperIdxB)
{
    // compute the incremental goodness
    double goodnessChange = 0;
    goodnessChange += sessionExchangeGoodness_nc2(trkA, trkB, timeA, timeB, paperIdxA, paperIdxB);
    goodnessChange += sessionExchangeGoodness_nc2(trkB, trkA, timeB, timeA, paperIdxB, paperIdxA);

    // make and return a new neighbour
    Neighboursingle ng(trkA, trkB, timeA, timeB, paperIdxA, paperIdxB, goodnessChange);
    return ng;
}

double SessionOrganizer::sessionExchangeGoodness_nc2(int trkA, int trkB, int timeA, int timeB, int paperIdxA, int paperIdxB)
{
    double delta = 0;
    int paperA, paperB;
    paperA = conference->getPaper(trkA, timeA, paperIdxA);
    if (timeA == timeB)
    {
        // process session B
        for (int p = 0; p < papersInSession; p++)
        {
            if (p == paperIdxB)
                continue;

            paperB = conference->getPaper(trkB, timeB, p);

            // add similarity
            //delta += 1 - distanceMatrix[paperA][paperB];
            delta += (1) - (tradeoffCoefficient + 1) * distanceMatrix[paperA][paperB];

            // subtract difference
            //delta -= tradeoffCoefficient * distanceMatrix[paperA][paperB];
        }

        // process session A
        for (int p = 0; p < papersInSession; p++)
        {
            if (p == paperIdxA)
                continue;

            paperB = conference->getPaper(trkA, timeA, p);

            // subtract similarity
            //delta -= 1 - distanceMatrix[paperA][paperB];
            delta += (-1) + (tradeoffCoefficient + 1) * distanceMatrix[paperA][paperB];

            // add difference
            //delta += tradeoffCoefficient * distanceMatrix[paperA][paperB];
        }
    }
    else
    {
        // both in different time slots

        // Process timeB
        for (int track = 0; track < parallelTracks; track++)
        {
            if (track == trkB)
            {
                for (int p = 0; p < papersInSession; p++)
                {
                    if (p == paperIdxB)
                        continue;
                    // add similarity
                    paperB = conference->getPaper(track, timeB, p);
                    delta += 1 - distanceMatrix[paperA][paperB];
                }
            }
            else
            {
                for (int p = 0; p < papersInSession; p++)
                {
                    // add distance
                    paperB = conference->getPaper(track, timeB, p);
                    delta += tradeoffCoefficient * distanceMatrix[paperA][paperB];
                }
            }
        }

        // Process timeA
        for (int track = 0; track < parallelTracks; track++)
        {
            if (track == trkA)
            {
                for (int p = 0; p < papersInSession; p++)
                {
                    if (p == paperIdxA)
                        continue;
                    // subtract similarity
                    paperB = conference->getPaper(track, timeA, p);
                    delta -= 1 - distanceMatrix[paperA][paperB];
                }
            }
            else
            {
                for (int p = 0; p < papersInSession; p++)
                {
                    // subtract distance
                    paperB = conference->getPaper(track, timeA, p);
                    delta -= tradeoffCoefficient * distanceMatrix[paperA][paperB];
                }
            }
        }
    }

    return delta;
}

void SessionOrganizer::gotoNeighbour_nc2(Neighboursingle ngh)
{
    // Convert the current conference state to that represented by the given neighbour
    int paperA = conference->getPaper(ngh.getTrackA(), ngh.getTimeA(), ngh.getPaperIdxA());
    int paperB = conference->getPaper(ngh.getTrackB(), ngh.getTimeB(), ngh.getPaperIdxB());

    conference->setPaper(ngh.getTrackA(), ngh.getTimeA(), ngh.getPaperIdxA(), paperB);
    conference->setPaper(ngh.getTrackB(), ngh.getTimeB(), ngh.getPaperIdxB(), paperA);
}

void SessionOrganizer::updateMaximum(double newMaxGoodness)
{
    // Update the bestConference to present conference
    for (int i = 0; i < sessionsInTrack; i++)
    {
        for (int j = 0; j < parallelTracks; j++)
        {
            for (int k = 0; k < papersInSession; k++)
            {
                int paper = conference->getPaper(j, i, k);
                bestConference->setPaper(j, i, k, paper);
            }
        }
    }

    // Update the maxGoodness value
    maxGoodness = newMaxGoodness;
}

bool SessionOrganizer::isOutOfTime()
{
    // Checks if our time is going to be over

    int timeElapsed = time(NULL) - organizerStartTime;
    int timeRemaining = processingTimeInSeconds - timeElapsed;

    if (totalPapers <= 200)
    {
        // Give extra 2 secs
        return (timeRemaining <= 2);
    }
    else if (totalPapers <= 750)
    {
        // Give extra 4 secs
        return (timeRemaining <= 4);
    }
    else if (totalPapers <= 1500)
    {
        // Give extra 7 secs
        return (timeRemaining <= 7);
    }
    else if (totalPapers <= 2500)
    {
        // Give extra 10 secs
        return (timeRemaining <= 10);
    }
    else
    {
        // Give extra 15 secs
        return (timeRemaining <= 15);
    }
}

void SessionOrganizer::localSearch_nc2()
{
    double min_Admissible_Val = 0.01;
    int iter = 1;
    double score = scoreOrganization();
    cout << "score:" << score << endl;

    int start_time = time(NULL);

    // double max_prev_itr = 0;
    int prob_gen;
    bool close = false;

    int max_ch = 25, cnt_ch = 0, max_score = 0;

    double min_val_change = -8.0;

    int hillCount = 1;

    int max_pap_limit = 700;

    bool type =  (totalPapers < max_pap_limit) ;
    //
    // The while loop needs to be made dependent on time.
    //

    while (iter++)
    {
        // Break if out of time
        if (isOutOfTime())
        {
            break;
        }

        prob_gen = min(iter / 100, 50);

        vector<Neighboursingle> neighbours;

        if (type) 
        {
            neighbours = getNeighbours_nc2(close, prob_gen);
        }
            
        else neighbours = getNeighbours();


        int ngh_size = neighbours.size();
        int arr[ngh_size + 1];

        for (int j = 0; j < ngh_size; j++)
        {
            arr[j] = j;
        }

        for (int i = ngh_size - 1; i >= 0; --i)
        {
            //generate a random number [0, i]
            int j = rand() % (i + 1);

            //swap the last element with element at random index
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }

        if (neighbours.size() < 1)
        {
            // no neighbours
            continue;
        }

        int max_nh_idx = arr[0];
        // bool ab = true;

        for (int nh = 0; nh < neighbours.size(); nh++)
        {

            if (neighbours.at(arr[nh]).getGoodInc() >= (neighbours.at(max_nh_idx).getGoodInc()) && neighbours.at(arr[nh]).getGoodInc() > min_Admissible_Val)
            {
                max_nh_idx = arr[nh];
                if(type)
                {
                    gotoNeighbour_nc2(neighbours.at(max_nh_idx));
                }
                //neighbours.at(nh).printNeighbour();
                // ab = false;                
            }

            if(!type) 
            {
                gotoNeighbour_nc2(neighbours.at(max_nh_idx));
            }


        }

        // max_prev_itr = neighbours.at(max_nh_idx).getGoodInc();

        if (neighbours.at(max_nh_idx).getGoodInc() <= min_Admissible_Val)
        {
            // if (close && cnt_ch < max_ch)
            if (close)
            {
                // Reached a hill
                double score = scoreOrganization();
                if (score > maxGoodness)
                {
                    // Update the global schedule
                    updateMaximum(score);
                }
                // cout << endl;
                // cout << "Hill: " << (hillCount++) << " | Iter: " << iter << " | score: " << score << " | max score: " << maxGoodness << " | total time: " << (time(NULL) - start_time) << endl;
                // cout << "iter: " << iter << ", score:" << score << " , increment: " << neighbours.at(max_nh_idx).getGoodInc() << " " << neighbours.at(max_nh_idx).getType() << endl;
                // cout << "Took " << iter << " steps in time: " << (time(NULL) - start_time) << " secs" << endl;

                // Make transitions to neighbours with changes greater than min_val_change
                for (int nh = 1; nh < neighbours.size(); nh++)
                {
                    if (neighbours.at(arr[nh]).getGoodInc() > min_val_change)
                    {
                        max_nh_idx = arr[nh];
                        gotoNeighbour_nc2(neighbours.at(max_nh_idx));
                    }
                }
                // cnt_ch++;
            }

            // else if (close)
            // break;
        }
        // else
        // {
        //     // goto neighbour
        //     if (ab)
        //     {
        //         //Seems useless but am not changing it since code working fine!

        //         //neighbours.at(max_nh_idx).printNeighbour();
        //         gotoNeighbour_nc2(neighbours.at(max_nh_idx));
        //     }
        // }

        score = scoreOrganization();
        // min_Admissible_Val = max(min_Admissible_Val, score/5000 );

        cout << "iter: " << iter << ", score:" << score << " , increment: " << neighbours.at(max_nh_idx).getGoodInc() << " " << neighbours.at(max_nh_idx).getType() << endl;
        cout << "Took " <<  iter << " steps in time: " << (time(NULL) - start_time) << " secs" << endl;
    }

    // Check for score
    score = scoreOrganization();
    if (score > maxGoodness)
    {
        // Update the global schedule
        updateMaximum(score);
    }

    // Print time taken
    int end_time = time(NULL);
    double timeConsumed = (end_time - organizerStartTime) / 60.0;
    cout << "Score: " << maxGoodness << " | Time: " << timeConsumed << " min / " << processingTimeInMinutes << " min :)"  << endl;
    // cout << "Took " << (iter - 1) << " steps in time: " << (end_time - start_time) << " secs" << endl;
}