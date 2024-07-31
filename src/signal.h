#include <ADC.h>

#define CYCLES_IN_MASK 5
#define SAMPLING_BUFFER_SIZE 480 //(10 + 1) * 200 CYCLES_IN_MASK + 1 pour le shifting * SAMPLING_FREQUENCY/SIGNAL1_FREQUENCY

//Lowest signal frequency = 1kHz -> longest cycles is 200 samples

//Generate a square signale wit 50% duty cycle and a given periode
double* generateMask(int maskLength, double periodeLength);
double correlation(double* signal1, int* signal2, int length, int shift, int signalIndex);
double correlationShift(double* signal1, int* signal2, int length, int max_shift, int signalIndex);
void setupADC(ADC* adc);