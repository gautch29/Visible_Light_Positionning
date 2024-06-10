#include <math.h>
#include <Arduino.h>
#include "mathK.h"

#define DEBUG

// System parameters
double h = 2;    // Height of the receiver
double D = 3;    // Distance between reference points
double be = 60 * PI / 180; // Beam elevation angle
double al[3] = {0, 120 * PI / 180, -120 * PI / 180}; // Beam azimuth angles
double si = sqrt(2.5e-40); // Measurement noise standard deviation
int nIter = 100; // Number of iterations
double T = 0.01; // Time step
double ka = 1; // Placeholder for constant 'ka'
double m = 1; // Placeholder for constant 'm'
double P = 1; // Placeholder for constant 'P'

// State transition matrix
double A[6][6] = {
    {1, T, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 0},
    {0, 0, 1, T, 0, 0},
    {0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 1, T},
    {0, 0, 0, 0, 0, 1}
};

// Process and measurement noise covariances
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
    {-D, 0, D, -D, 0, D, -D, 0, D},
    {D, D, D, 0, 0, 0, -D, -D, -D},
    {h, h, h, h, h, h, h, h, h}
};

// Placeholder for beam vectors
double v[3][3];

// Main EKF loop
void setup() {
    Serial.begin(115200);
    Serial.println("Starting EKF loop");

    // Initialize measurement noise covariance R
    for (int i = 0; i < 3*K; i++) {
        for (int j = 0; j < 3*K; j++) {
            R[i][j] = (i == j) ? 0.01 : 0;
        }
    }

    // Generate beam vectors
    for (int l = 0; l < 3; l++) {
        double tmp[3][3] = {{cos(al[l]), sin(al[l]), 0},
                           {-sin(al[l]), cos(al[l]), 0},
                           {0, 0, 1}};

        double tmp2[3][3] = {{cos(be), 0, -sin(be)},
                            {0, 1, 0},
                            {sin(be), 0, cos(be)}};
        double tmp_v[3] = {1, 0, 0};

        double res_v[3];
        double res_v2[3];

        matrix_vector_mult(tmp2, tmp_v, res_v);

        matrix_vector_mult(tmp, res_v, res_v2);

        for (int i = 0; i < 3; i++) {
            v[i][l] = res_v2[i];
        }
    }
    // Initialize state and covariance
    double p0[3];
    double theta0;

    p0[0] = -3.6287;
    p0[1] = 0.2866;
    p0[2] = 0;
    theta0 = 2.7332;

    // Initial measurement matrix H0 and initial range measurements r0
    double H0[3][K];
    comp_H0(ka, m, v, p, p0, theta0, H0);

    double r0[3*K] = {0.057357, 0.080822, 0.012166, 0.0040265, 0.026574, 0.0039414, 0, 0.0068034, 0.001001, 0.27668, 0.3789, 0.38121, 0.00013447, 0.041717, 0.028002, 0, 0.0074713, 0.0039121, 0.013699, 0.0174, 0.058898, 0, 0.011294, 0.01827, 0, 0.0039565, 0.0043419};
    /*
    for (int k = 0; k < K; k++) {
        double eyeK[9] = {0};
        eyeK[k] = 1;

        //H0 * P
        double H0P[3][K];
        for(int i = 0; i < 3; i++) {
            for (int j = 0; j < K; j++) {
                H0P[i][j] = 0;
                for (int l = 0; l < 3; l++) {
                    H0P[i][j] += H0[i][l] * P;
                }
            }
        }

        //H0 * P * eyeK(:, k) = H0 * P * eyeK[k]
        double H0PeyeK[3];
        for (int i = 0; i < 3; i++) {
            H0PeyeK[i] = 0;
            for (int j = 0; j < K; j++) {
                H0PeyeK[i] += H0P[i][j] * eyeK[j];
            }
        }

        for (int i = 0; i < 3; i++) {
            r0[3*k + i] = H0PeyeK[i];
        }
    }
    */

    double xp[6] = {0, 0, 0, 0, 0, 0};
    double Pp[6][6] = {
        {0.01, 0, 0,    0, 0,    0},
        {0,    1, 0,    0, 0,    0},
        {0,    0, 0.01, 0, 0,    0},
        {0,    0, 0,    1, 0,    0},
        {0,    0, 0,    0, 0.01, 0},
        {0,    0, 0,    0, 0,    1}
    };

    double r[3*K]; // Simulated measurement
    double Hl[3*K][6]; //Measurement matrix linearization
    double Kg[6][3*K]; //Kalman gain
    double xm[6]; //State prediction
    double Pm[6][6]; //Covariance prediction
    double rH[3*K];

    for (int n = 0; n < nIter; n++) {
        // Simulated measurement noise
        for (int i = 0; i < 3*K; i++) {
            r[i] = r0[i] + si * random(1000)/1000.0;
        }
        // State prediction
        state_prediction(A, xp, xm);

        // Covariance prediction
        cov_prediction(A, Pp, Q, Pm);

        // Measurement matrix linearization
        linear_H0(xm[0], xm[2], xm[4], Hl, p, v, h, ka, m);

        kalmanGain(Hl, Pm, R, Kg); // Kalman gain

        // Measurement prediction****************************
        double a[3] = {xm[0], xm[2], 0};

        comp_H0(ka, m, v, p, a, xm[4], H0);

        for (int k = 0; k < K; k++) {
            for (int l = 0; l < 3; l++) {
                rH[3*k + l] = H0[l][k];
            }
        }

        state_update(xm, Kg, r, rH, xp); // State update

        cov_update(Kg, Hl, Pm, Pp, R); // Covariance update
 

        // Output state estimate
        Serial.print("State estimate at iteration ");
        Serial.print(n);
        Serial.println(": ");
        for (int i = 0; i < 6; i++) {
            Serial.print(xp[i]);
            Serial.print(" ");
        }
        Serial.println();
    }
}

void loop() {
    // Empty loop
}