#include <ADC.h>

#define CYCLES_IN_MASK 5
#define SAMPLING_BUFFER_SIZE 480 //(10 + 1) * 200 CYCLES_IN_MASK + 1 pour le shifting * SAMPLING_FREQUENCY/SIGNAL1_FREQUENCY

//Lowest signal frequency = 1kHz -> longest cycles is 200 samples

//Generate a square signale wit 50% duty cycle and a given periode
int* generateMask(int cycles, int periode);
double correlation(int* signal1, int* signal2, int length, int shift, int signalIndex);
double correlationShift(int* signal1, int* signal2, int length, int max_shift, int signalIndex);
void setupADC(ADC* adc);