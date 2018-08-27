#include "Neighbour.h"

Neighbour::Neighbour (bool neighbourType, int trkA, int trkB, int timeA, int timeB, int exSize, double goodInc) {
    // Initialise the neighbour
    type = neighbourType;

    trackA = trkA;
    trackB = trkB;
    timeSlotA = timeA;
    timeSlotB = timeB;

    exchangeSize = exSize;
    goodnessIncrement = goodInc;
}

bool Neighbour::getType () {
    return type;
}

int Neighbour::getTrackA () {
    return trackA;
}

int Neighbour::getTrackB () {
    return trackB;
}

int Neighbour::getTimeA () {
    return timeSlotA;
}

int Neighbour::getTimeB () {
    return timeSlotB;
}

int Neighbour::getExSize () {
    return exchangeSize;
}

double Neighbour::getGoodInc () {
    return goodnessIncrement;
}