/* 
 * File:   main.cpp
 * Author: Kapil Thakkar
 *
 */

#include <cstdlib>

#include "SessionOrganizer.h"

using namespace std;

/*
 * 
 */
int main(int argc, char **argv)
{
    // Parse the input.
    if (argc < 3)
    {
        cout << "Missing arguments\n";
        cout << "Correct format : \n";
        cout << "./main <input_filename> <output_filename>";
        exit(0);
    }
    string inputfilename(argv[1]);

    // Initialize the conference organizer.
    SessionOrganizer *organizer = new SessionOrganizer(inputfilename);
    srand(time(NULL));
    // Organize the papers into tracks based on similarity.
    // organizer->organizePapers ( );
     organizer->initializeConference();
     //organizer->initializeConferenceRandomly();
    // organizer->localSearch();
    organizer->localSearch_nc2();

    organizer->printSessionOrganiser(argv[2]);

    // Score the organization against the gold standard.
    // double score = organizer->scoreOrganization();
    // cout << "score:" << score << endl;

    return 0;
}
