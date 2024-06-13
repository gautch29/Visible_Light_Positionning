#include "Arduino.h"
#include "signal.h"

//Gnerate a square signale wit 50% duty cycle and a given periode
double* generateMask(int cycles, int periode){
    double* mask = new double[cycles*periode];
    for (int i = 0; i < cycles*periode; i++){
        if (i % periode < periode / 2){
            mask[i] = 1;
        }
        else{
            mask[i] = -1;
        }
    }

    return mask;
}

//Compute the correlation between two signals with a shift
double correlation(double* mask, double* signal, int length, int shift){
    double sum = 0;
    double maskMean = 0;
    double signalMean = 0;
    for (int i = 0; i < length; i++){
        maskMean += mask[i];
    }
    maskMean /= length;

    for (int i = 0; i < length; i++){
        signalMean += signal[i];
    }
    signalMean /= length;

    for (int i = 0; i < length; i++){
        sum += (mask[i] - maskMean) * (signal[(i + shift) % length] - signalMean);
    }

    return sum/length;
}

//Compute the correlation between two signals with a shift
double correlationShift(double* mask, double* signal, int length, int max_shift){
    double max = 0;
    for (int i = 0; i < max_shift; i++){
        double c = correlation(mask, signal, length, i);
        if (c > max){
            max = c;
        }
    }
    return max;
}