#include "Neighboursingle.h"

Neighboursingle::Neighboursingle (int trkA, int trkB, int timeA, int timeB, int paperA, int paperB, double goodInc) {
    // Initialise the neighbour

    trackA = trkA;
    trackB = trkB;
    timeSlotA = timeA;
    timeSlotB = timeB;

    paperIdxA = paperA;
    paperIdxB = paperB;
    goodnessIncrement = goodInc;
}

bool Neighboursingle::getType () {
    return type;
}

int Neighboursingle::getTrackA () {
    return trackA;
}

int Neighboursingle::getTrackB () {
    return trackB;
}

int Neighboursingle::getTimeA () {
    return timeSlotA;
}

int Neighboursingle::getTimeB () {
    return timeSlotB;
}

int Neighboursingle::getPaperIdxA () {
    return paperIdxA;
}

int Neighboursingle::getPaperIdxB () {
    return paperIdxB;
}

double Neighboursingle::getGoodInc () {
    return goodnessIncrement;
}

void Neighboursingle::printNeighbour () {
    cout << timeSlotA << "," << trackA << "," << paperIdxA << " | " << timeSlotB << "," << trackB << "," << paperIdxB << " | " << goodnessIncrement << endl;
}