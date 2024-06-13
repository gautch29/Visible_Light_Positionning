#include <math.h>
#include <Arduino.h>
#include "kalman.h"
#include "signal.h"

void readInputs();

// System parameters ***************************************************************************************************
double h = 1;    // Height of the receiver
double D = 2;    // Distance between reference points
double be = 60 * PI / 180; // Beam elevation angle
double al[3] = {0, 120 * PI / 180, -120 * PI / 180}; // Beam azimuth angles
double si = sqrt(2.5e-40); // Measurement noise standard deviation
int nIter = 100; // Number of iterations
double T = 0.01; // Time step
double ka = 1; // Placeholder for constant 'ka'
double m = 1; // Placeholder for constant 'm'
double P = 1; // Placeholder for constant 'P'

// Kalman filter parameters *******************************************************************************************
// State transition matrix
double A[6][6] = {
    {1, T, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 0},
    {0, 0, 1, T, 0, 0},
    {0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 1, T},
    {0, 0, 0, 0, 0, 1}
};
// Process noise covariance
double Q[6][6] = {
    {0.01, 0, 0, 0, 0, 0},
    {0, 0.01, 0, 0, 0, 0},
    {0, 0, 0.01, 0, 0, 0},
    {0, 0, 0, 0.01, 0, 0},
    {0, 0, 0, 0, 0.01, 0},
    {0, 0, 0, 0, 0, 0.01}
};
// Measurement noise covariance, initialized later
double R[3*K][3*K];
// Reference point positions
double p[3][K] = {
    {0, D},
    {0, 0},
    {h, h}
};
// Placeholder for beam vectors
double v[3][3];
// Simulated measurement
double r[3*K];
//State prediction
double xm[6];
//Covariance prediction
double Pm[6][6]; 
//State update
double xp[6] = {0, 0, 0, 0, 0, 0};
//Covariance update
double Pp[6][6] = {
    {0.01, 0, 0,    0, 0,    0},
    {0,    1, 0,    0, 0,    0},
    {0,    0, 0.01, 0, 0,    0},
    {0,    0, 0,    1, 0,    0},
    {0,    0, 0,    0, 0.01, 0},
    {0,    0, 0,    0, 0,    1}
};
double H0[3][K];
double r0[3*K];

//Timing parameters
IntervalTimer timerSampling;
unsigned long last_millis_correlation = 0;

//Sampling variables
double* signal1;

//Variables for correlation
double* maskSignal1;
double c1;
double* maskSignal2;
double c2;

// Main EKF loop
void setup() {
    Serial.begin(115200);

    signal1 = new double[SAMPLING_BUFFER_SIZE];
    // Initialize measurement noise covariance R
    for (int i = 0; i < 3*K; i++) {
        for (int j = 0; j < 3*K; j++) {
            R[i][j] = (i == j) ? 0.01 : 0;
        }
    }

    // Generate beam vectors
    beamVectors(v, al, be);
    
    // Initialize state
    double p0[3] = {0.5,0,0};
    double theta0 = 0;

    // Initial measurement matrix H0 and initial range measurements r0
    comp_H0(ka, m, v, p, p0, theta0, H0);

    //R0 calculation
    r0Comp(P, H0, r0);

    for (int n = 0; n < nIter; n++) {
        // Simulated measurement noise
        for (int i = 0; i < 3*K; i++) {
            r[i] = r0[i] + si * random(1000)/1000.0;
        }

        //Kalman Filter
        kalmanFilter(A, xp, Pp, Q, r0, R, p, v, h, ka, m, xm, Pm, r, H0);
 
        // Output state estimate
        printKalmanResults(xp, n);
    }

    maskSignal1 = generateMask(CYCLES_IN_MASK, SAMPLING_FREQUENCY/SIGNAL1_FREQUENCY);
    maskSignal2 = generateMask(CYCLES_IN_MASK, SAMPLING_FREQUENCY/SIGNAL2_FREQUENCY);

    //Setting up the sampling timer
    timerSampling.begin(readInputs, 62);
    timerSampling.priority(0);

}

void loop() {

    //Correlation
    if (millis() - last_millis_correlation > 100) {
        //Correlation
        last_millis_correlation = millis();
        //Copy signal
        double* signalCopy = new double[SAMPLING_BUFFER_SIZE];
        for (int i = 0; i < SAMPLING_BUFFER_SIZE; i++){
            signalCopy[i] = signal1[i];
        }

        c1 = correlationShift(maskSignal1, signalCopy, CYCLES_IN_MASK * SAMPLING_FREQUENCY/SIGNAL1_FREQUENCY, SAMPLING_FREQUENCY/SIGNAL1_FREQUENCY);
        c2 = correlationShift(maskSignal2, signalCopy, CYCLES_IN_MASK * SAMPLING_FREQUENCY/SIGNAL2_FREQUENCY, SAMPLING_FREQUENCY/SIGNAL2_FREQUENCY);

        Serial.print(" Correlation 1 :");
        Serial.print(c1);
        Serial.print(" Correlation 2 :");
        Serial.print (c2);
        Serial.print(" Signal 1 :");
        Serial.println(signalCopy[SAMPLING_BUFFER_SIZE - 1]);

        delete[] signalCopy;
    }
}

void readInputs(){
    for(int i = 0; i < SAMPLING_BUFFER_SIZE - 1; i++){
        signal1[i] = signal1[i+1];
    }
    signal1[SAMPLING_BUFFER_SIZE - 1] = analogRead(A0);
}