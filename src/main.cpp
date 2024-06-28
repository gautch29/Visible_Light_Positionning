#include <math.h>
#include <Arduino.h>
#include "kalmanFilter.h"
#include "signal.h"

void readInputs();
void printCorrelations();

// System parameters ***************************************************************************************************
double h = 0.35;    // Height of the receiver
double D = 0.4;    // Distance between reference points
double be = 60 * PI / 180; // Beam elevation angle
double al[3] = {0, 120 * PI / 180, -120 * PI / 180}; // Beam azimuth angles
double p0[3] = {0,0,0};
double theta0 = 0;
double ka = 5.46; // Placeholder for constant 'ka'
double m = 2; // Placeholder for constant 'm'
double P = 1; // Placeholder for constant 'P'
double p[3][K] = { // Reference point positions
    {0, D, 0, D},
    {0, 0, D, D},
    {h, h, h, h}
};

int signalFrequencies[K] = {1000, 2000, 4000, 8000};

//Sampling variables ***************************************************************************************************
double signal[3][SAMPLING_BUFFER_SIZE]; //Buffer for the signal
IntervalTimer timerSampling; //Timer for sampling

//Correlation variables ************************************************************************************************
double *maskSignal[K];
double correlations[K][3];
unsigned long last_millis_correlation = 0;


//Kalman filter ********************************************************************************************************
KalmanFilter kF = KalmanFilter(p0, theta0, P, al, be, ka, m, h, p);
unsigned long last_millis_ekf = 0;

//Loop timing **********************************************************************************************************
int EKF_period = 10; //Period of the EKF in ms
int correlation_period = 10; //Period of the correlation in ms

void setup() {
    Serial.begin(115200);

    //Generate correlation masks
    for(int i = 0; i < K; i++){
        maskSignal[i] = generateMask(CYCLES_IN_MASK, SAMPLING_FREQUENCY/signalFrequencies[i]);
    }

    //Setting up the sampling timer
    timerSampling.begin(readInputs, 1000000/SAMPLING_FREQUENCY);
    timerSampling.priority(0);

}

void loop() {

    //Correlation
    if (millis() - last_millis_correlation > correlation_period) {
        //Update last time
        last_millis_correlation = millis();

        //Copy signal to avoid interrupr changing it during computation
        double signalCopy[3][SAMPLING_BUFFER_SIZE];   
        for(int i = 0; i<3 ; i++){
            for (int j = 0; j < SAMPLING_BUFFER_SIZE; j++){
                signalCopy[i][j] = signal[i][j];
            }
        }
        
        //Compute correlation (i is the signal, j is the sensor)
        for(int i = 0; i < K; i++){
            for(int j = 0; j < 3; j++){
                correlations[i][j] = correlationShift(maskSignal[i], signalCopy[j], CYCLES_IN_MASK * SAMPLING_FREQUENCY/signalFrequencies[i], SAMPLING_FREQUENCY/signalFrequencies[i]);
            }
        }

        //printCorrelations(); //Print correlations
    }

    //EKF
    if (millis() - last_millis_ekf > EKF_period) {

        last_millis_ekf = millis(); //Update last time

        //Update measurements
        double r[3*K] = {correlations[0][0], correlations[0][1], correlations[0][2], correlations[1][0], correlations[1][1], correlations[1][2], correlations[2][0], correlations[2][1], correlations[2][2], correlations[3][0], correlations[3][1], correlations[3][2]};
        kF.updateMeasurements(r);

        kF.compute(); //Compute EKF

        //Forcing y and theta to 0
        //kF.xp[2] = 0;
        //kF.xp[4] = 0;

        kF.print();
    }
}

void readInputs(){
    //Shift the buffer
    for(int i = 0; i < SAMPLING_BUFFER_SIZE - 1; i++){
        signal[0][i] = signal[0][i+1];
        signal[1][i] = signal[1][i+1];
        signal[2][i] = signal[2][i+1];
    }
    //Read the new values
    signal[0][SAMPLING_BUFFER_SIZE - 1] = analogRead(A0);
    signal[1][SAMPLING_BUFFER_SIZE - 1] = analogRead(A1);
    signal[2][SAMPLING_BUFFER_SIZE - 1] = analogRead(A2);
}

void printCorrelations(){
    for(int i = 0; i < K; i++){
        for(int j = 0; j < 3; j++){
            Serial.print(correlations[i][j]);
            Serial.print(" ");
        }
    }
    Serial.println();
}