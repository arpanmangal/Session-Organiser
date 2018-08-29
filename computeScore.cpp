// C++ file to compute score of output file
#include <bits/stdc++.h>
using namespace std;

int main(int argc, char **argv)
{
    ifstream file;
    string inputfile(argv[1]);
    file.open(inputfile.c_str());

    double processingTimeInMinutes;
    int papersInSession;
    int parallelTracks;
    int sessionsInTrack;
    double tradeoffCoefficient;
    file >> processingTimeInMinutes >> papersInSession >> parallelTracks >> sessionsInTrack >> tradeoffCoefficient;
    int totalPapers = papersInSession * parallelTracks * sessionsInTrack;

    double **distanceMatrix = new double*[totalPapers];
    for (int i = 0; i < totalPapers; i++)
    {
        distanceMatrix[i] = new double[totalPapers];
        for (int j = 0; j < totalPapers; j++)
        {
            file >> distanceMatrix[i][j];
        }
    }

    cout << processingTimeInMinutes << papersInSession << parallelTracks << sessionsInTrack << tradeoffCoefficient << endl;

    // Read output file
    ifstream ofile;
    string outputfile(argv[2]);
    ofile.open(outputfile.c_str());

    vector<vector<vector<int>>> conference(sessionsInTrack);

    for (int t = 0; t < sessionsInTrack; t++)
    {
        conference[t].resize(parallelTracks);
        for (int p = 0; p < parallelTracks; p++)
        {
            conference[t][p].resize(papersInSession);
            for (int k = 0; k < papersInSession; k++)
            {
                ofile >> conference[t][p][k];
                // cout << conference[t][p][k];
            }
            if (p != parallelTracks - 1)
            {
                char pipe;
                ofile >> pipe;
            }
        }
    }

    // for (int t = 0; t < sessionsInTrack; t++)
    // {
    //     for (int p = 0; p < parallelTracks; p++)
    //     {
    //         for (int k = 0; k < papersInSession; k++)
    //         {
    //             cout << conference[t][p][k] << " ";
    //         }
    //         cout << " | ";
    //     }
    //     cout << endl;
    // }

    // Compute score
    // Sum of pairwise similarities per session.
    // double score = 0.0;
    // for (int t = 0; t < sessionsInTrack; t++)
    // {
    //     for (int p = 0; p < parallelTracks; p++)
    //     {
    //         for (int k = 0; k < papersInSession - 1; k++)
    //         {
    //             for (int k2 = k + 1; k2 < papersInSession; k2++)
    //             {

    //                 int A = conference[t][p][k];
    //                 int B = conference[t][p][k2];
    //                 score += 1 - distanceMatrix[A][B];
    //             }
    //         }
    //     }
    // }

    // for (int t = 0; t < sessionsInTrack; t++)
    // {
    //     for (int p = 0; p < parallelTracks; p++)
    //     {
    //         for (int p2 = p + 1; p < parallelTracks; p++)
    //         {
    //             for (int k1 = 0; k1 < papersInSession; k1++)
    //             {
    //                 for (int k2 = 0; k2 < papersInSession; k2++)
    //                 {
    //                     int A = conference[t][p][k1];
    //                     int B = conference[t][p2][k2];
    //                     score += tradeoffCoefficient * distanceMatrix[A][B];
    //                 }
    //             }
    //         }
    //     }
    // }

    double score1 = 0.0;
    for (int i = 0; i < parallelTracks; i++)
    {
        // Track tmpTrack = conference->getTrack(i);
        for (int j = 0; j < sessionsInTrack; j++)
        {
            // Session tmpSession = tmpTrack.getSession(j);
            for (int k = 0; k < papersInSession; k++)
            {
                int index1 = conference[j][i][k];
                for (int l = k + 1; l < papersInSession; l++)
                {
                    int index2 = conference[j][i][l];
                    score1 += 1 - distanceMatrix[index1][index2];
                }
            }
        }
    }

    // Sum of distances for competing papers.
    double score2 = 0.0;
    for (int i = 0; i < parallelTracks; i++)
    {
        // Track tmpTrack1 = conference->getTrack(i);
        for (int j = 0; j < sessionsInTrack; j++)
        {
            // Session tmpSession1 = tmpTrack1.getSession(j);
            for (int k = 0; k < papersInSession; k++)
            {
                int index1 = conference[j][i][k];

                // Get competing papers.
                for (int l = i + 1; l < parallelTracks; l++)
                {
                    // Track tmpTrack2 = conference->getTrack(l);
                    // Session tmpSession2 = tmpTrack2.getSession(j);
                    for (int m = 0; m < papersInSession; m++)
                    {
                        int index2 = conference[j][l][m];
                        score2 += distanceMatrix[index1][index2];
                    }
                }
            }
        }
    }
    double score3 = score1 + tradeoffCoefficient * score2;
    // return score;

    cout  << " | " << score3 << endl;
}