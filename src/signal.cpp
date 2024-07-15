#include "Arduino.h"
#include "signal.h"

//Gnerate a square signale wit 50% duty cycle and a given periode
int* generateMask(int cycles, int periode){
    int* mask = new int[cycles*periode];
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
double correlation(int* mask, int* signal, int length, int shift, int signalIndex){
    double sum = 0;
    double maskMean = 0;
    double signalMean = 0;
    for (int i = 0; i < length; i++){
        maskMean += mask[i];
    }
    maskMean /= length;

    for (int i = 0; i < length; i++){
        signalMean += signal[(i + signalIndex + 1) % SAMPLING_BUFFER_SIZE];
    }
    signalMean /= length;

    for (int i = 0; i < length; i++){
        sum += ((double)mask[i] - maskMean) * ((double)signal[(((i + signalIndex + 1) % SAMPLING_BUFFER_SIZE) + shift) % length] - signalMean);
    }

    return sum/length;
}

//Compute the correlation between two signals with a shift
double correlationShift(int* mask, int* signal, int length, int max_shift, int signalIndex){
    double max = 0;
    for (int i = 0; i < max_shift; i++){
        double c = correlation(mask, signal, length, i, signalIndex);
        if (c > max){
            max = c;
        }
    }
    return max;
}

void setupADC(ADC* adc){
    adc->adc0->setAveraging(1); // no averaging
    adc->adc0->setResolution(10); // 10 bits resolution
    adc->adc0->setConversionSpeed(ADC_CONVERSION_SPEED::VERY_HIGH_SPEED); // fastest conversion
    adc->adc0->setSamplingSpeed(ADC_SAMPLING_SPEED::VERY_HIGH_SPEED); // fastest sampling

    pinMode(A0, INPUT);
    pinMode(A1, INPUT);
    pinMode(A2, INPUT);
}