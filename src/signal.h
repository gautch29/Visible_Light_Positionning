#define SAMPLING_FREQUENCY 16000
#define SIGNAL1_FREQUENCY 1000
#define SIGNAL2_FREQUENCY 2000
#define SIGNAL1_SAMPLES_PER_CYCLE SAMPLING_FREQUENCY/SIGNAL1_FREQUENCY
#define SIGNAL2_SAMPLES_PER_CYCLE SAMPLING_FREQUENCY/SIGNAL2_FREQUENCY
#define CYCLES_IN_MASK 10

#define SAMPLING_BUFFER_SIZE 176 //(10+1) * 16 CYCLES_IN_MASK + 1 pour le shifting * SAMPLING_FREQUENCY/SIGNAL1_FREQUENCY

//Gnerate a square signale wit 50% duty cycle and a given periode
double* generateMask(int cycles, int periode);
double correlation(double* signal1, double* signal2, int length, int shift);
double correlationShift(double* signal1, double* signal2, int length, int max_shift);
