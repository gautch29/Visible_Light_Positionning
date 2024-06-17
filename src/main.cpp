#include <math.h>
#include <Arduino.h>
#include "kalman.h"
#include "signal.h"

void readInputs();

// System parameters ***************************************************************************************************
double h = 0.35;    // Height of the receiver
double D = 0.5;    // Distance between reference points
double be = 60 * PI / 180; // Beam elevation angle
double al[3] = {0, 120 * PI / 180, -120 * PI / 180}; // Beam azimuth angles
double si = sqrt(2.5e-40); // Measurement noise standard deviation
int nIter = 100; // Number of iterations
double T = 0.01; // Time step
double ka = 1; // Placeholder for constant 'ka'
double m = 2; // Placeholder for constant 'm'
double P = 1; // Placeholder for constant 'P'

int signalFrequencies[2] = {1000, 2000};

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
double signal[3][SAMPLING_BUFFER_SIZE];

//Variables for correlation
double *maskSignal[K];
double correlations[K][3];

// Main EKF loop
void setup() {
    Serial.begin(115200);

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
        //for (int i = 0; i < 3*K; i++) {
        //    r[i] = r0[i] + si * random(1000)/1000.0;
        //}

        //Kalman Filter
        //kalmanFilter(A, xp, Pp, Q, r0, R, p, v, h, ka, m, xm, Pm, r, H0);
 
        // Output state estimate
        //printKalmanResults(xp, n);
    }

    //Generate masks
    for(int i = 0; i < K; i++){
        maskSignal[i] = generateMask(CYCLES_IN_MASK, SAMPLING_FREQUENCY/signalFrequencies[i]);
    }

    //Setting up the sampling timer
    timerSampling.begin(readInputs, 62);
    timerSampling.priority(0);

}

void loop() {

    //Correlation
    if (millis() - last_millis_correlation > 10) {
        //Update last time
        last_millis_correlation = millis();

        //Copy signal to avoid it changing because of interrupt during computation
        double signalCopy[3][SAMPLING_BUFFER_SIZE];   
        for(int i = 0; i<3 ; i++){
            for (int j = 0; j < SAMPLING_BUFFER_SIZE; j++){
                signalCopy[i][j] = signal[i][j];
            }
        }
        
        //Compute correlation
        //i is the signal, j is the sensor
        for(int i = 0; i < K; i++){
            for(int j = 0; j < 3; j++){
                correlations[i][j] = correlationShift(maskSignal[i], signalCopy[j], CYCLES_IN_MASK * SAMPLING_FREQUENCY/signalFrequencies[i], SAMPLING_FREQUENCY/signalFrequencies[i]);
            }
        }

        //Scaling correlation
        for(int i = 0; i < K; i++){
            for(int j = 0; j < 3; j++){
                correlations[i][j] = correlations[i][j]/3;
            }
        }

/*
        Serial.print("Correlation 1: ");
        Serial.print(correlations[0][0]);
        Serial.print(" Correlation 2: ");
        Serial.print(correlations[0][1]);
        Serial.print(" Correlation 3: ");
        Serial.print(correlations[0][2]);
        Serial.print(" Correlation 4: ");
        Serial.print(correlations[1][0]);
        Serial.print(" Correlation 5: ");
        Serial.print(correlations[1][1]);
        Serial.print(" Correlation 6: ");
        Serial.print(correlations[1][2]);
*/
        /*
        Serial.print(" Signal 1: ");
        Serial.print(signalCopy[0][SAMPLING_BUFFER_SIZE - 1]);
        Serial.print(" Signal 2: ");
        Serial.print(signalCopy[1][SAMPLING_BUFFER_SIZE - 1]);
        Serial.print(" Signal 3: ");
        Serial.print(signalCopy[2][SAMPLING_BUFFER_SIZE - 1]);
        */
        Serial.println();


        //Lampe 1
        //Capteur 1
            r[0] = correlations[0][0];
        //Capteur 2
            r[1] = correlations[0][1];
        //Capteur 3
            r[2] = correlations[0][2];

        //Lampe 2
        //Capteur 1
            r[3] = correlations[1][0];
        //Capteur 2
            r[4] = correlations[1][1];
        //Capteur 3
            r[5] = correlations[1][2];

        //Kalman Filter
        kalmanFilter(A, xp, Pp, Q, R, p, v, h, ka, m, xm, Pm, r, H0);
 
        //Forcing y and theta to 0
        xp[2] = 0;
        xp[4] = 0;

        // Output state estimate
        printKalmanResults(xp);

    }
}

void readInputs(){
    for(int i = 0; i < SAMPLING_BUFFER_SIZE - 1; i++){
        signal[0][i] = signal[0][i+1];
        signal[1][i] = signal[1][i+1];
        signal[2][i] = signal[2][i+1];
    }
    signal[0][SAMPLING_BUFFER_SIZE - 1] = analogRead(A0);
    signal[1][SAMPLING_BUFFER_SIZE - 1] = analogRead(A1);
    signal[2][SAMPLING_BUFFER_SIZE - 1] = analogRead(A2);
}