#include "Neighbour.h"

Neighbour::Neighbour (bool neighbourType, int typeIndex, int swapIdx1, int swapIdx2, int exSize, double goodInc) {
    // Initialise the neighbour
    type = neighbourType;

    if (type) {
        // A parallel neighbour, with same time slot
        timeSlotA = typeIndex;
        timeSlotB = typeIndex;
        trackA = swapIdx1;
        trackB = swapIdx2;
    } else {
        // A time neighbour, with same track
        timeSlotA = swapIdx1;
        timeSlotB = swapIdx2;
        trackA = typeIndex;
        trackB = typeIndex;
    }

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