#include <Arduino.h>
#include <math.h>
#include "mathK.h"

// System parameters
float h = 2;    // Height of the receiver
float D = 3;    // Distance between reference points
float be = 60 * PI / 180; // Beam elevation angle
float al[3] = {0, 120 * PI / 180, -120 * PI / 180}; // Beam azimuth angles
float si = sqrt(2.5e-40); // Measurement noise standard deviation
int nIter = 3; // Number of iterations
float T = 0.01; // Time step
float ka = 1; // Placeholder for constant 'ka'
float m = 1; // Placeholder for constant 'm'
float P = 1; // Placeholder for constant 'P'

// State transition matrix
float A[6][6] = {
    {1, T, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 0},
    {0, 0, 1, T, 0, 0},
    {0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 1, T},
    {0, 0, 0, 0, 0, 1}
};

// Process and measurement noise covariances
float Q[6][6] = {
    {0.01, 0, 0, 0, 0, 0},
    {0, 0.01, 0, 0, 0, 0},
    {0, 0, 0.01, 0, 0, 0},
    {0, 0, 0, 0.01, 0, 0},
    {0, 0, 0, 0, 0.01, 0},
    {0, 0, 0, 0, 0, 0.01}
};

// Measurement noise covariance, initialized later
float R[3*K][3*K];

// Reference point positions
float p[3][K] = {
    {-D, 0, D, -D, 0, D, -D, 0, D},
    {D, D, D, 0, 0, 0, -D, -D, -D},
    {h, h, h, h, h, h, h, h, h}
};

// Placeholder for beam vectors
float v[3][3];

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
    for (int l = 0; l < 3; l++) {
        float tmp[3][3] = {{cos(al[l]), sin(al[l]), 0},
                           {-sin(al[l]), cos(al[l]), 0},
                           {0, 0, 1}};

        float tmp2[3][3] = {{cos(be), 0, -sin(be)},
                            {0, 1, 0},
                            {sin(be), 0, cos(be)}};
        float tmp_v[3] = {1, 0, 0};

        float res_v[3];
        float res_v2[3];

        matrix_vector_mult(tmp2, tmp_v, res_v);

        matrix_vector_mult(tmp, res_v, res_v2);

        for (int i = 0; i < 3; i++) {
            v[i][l] = res_v2[i];
        }
    }
    // Initialize state and covariance
    float p0[3];
    float theta0;

    p0[0] = -3.6287;
    p0[1] = 0.2866;
    p0[2] = 0;
    theta0 = 2.7332;

    // Initial measurement matrix H0 and initial range measurements r0
    float H0[3][K];
    comp_H0(ka, m, v, p, p0, theta0, H0);
    //Print out the initial measurement matrix
    Serial.println("Initial measurement matrix H0: ");
    for (int k = 0; k < K; k++) {
        for (int l = 0; l < 3; l++) {
            Serial.print(H0[l][k]);
            Serial.print(" ");
        }
        Serial.println();
    }

    float r0[3*K];

    for (int k = 0; k < K; k++) {
        for (int l = 0; l < 3; l++) {
            r0[3*k + l] = 0;
            for (int i = 0; i < 3; i++) {
                r0[3*k + l] += H0[l][k] * P;
            }
        }
    }

    float xp[6] = {0, 0, 0, 0, 0, 0};
    float Pp[6][6] = {
        {0.01, 0, 0,    0, 0,    0},
        {0,    1, 0,    0, 0,    0},
        {0,    0, 0.01, 0, 0,    0},
        {0,    0, 0,    1, 0,    0},
        {0,    0, 0,    0, 0.01, 0},
        {0,    0, 0,    0, 0,    1}
    };

    float r[3*K]; // Simulated measurement
    float Hl[3*K][6]; //Measurement matrix linearization
    float Kg[6][3*K]; //Kalman gain
    float xm[6]; //State prediction
    float Pm[6][6]; //Covariance prediction
    float rH[3*K];

    for (int n = 0; n < nIter; n++) {
        // Simulated measurement noise
        for (int i = 0; i < 3*K; i++) {
            r[i] = r0[i] + si * random(1000)/1000.0;
        }

        //Print out the simulated measurement
        Serial.print("Simulated measurement at iteration ");
        Serial.print(n);
        Serial.println(": ");
        for(int i = 0;i < 3; i++){
            for(int j = 0; j < K; j++){
                Serial.print(r[3*j + i]);
                Serial.print(" ");
            }
            Serial.println();
        }
        Serial.println();

        // State prediction
        state_prediction(A, xp, xm);

        //Print out the state prediction
        Serial.print("State prediction at iteration ");
        Serial.print(n);
        Serial.println(": ");
        for (int i = 0; i < 6; i++) {
            Serial.print(xm[i]);
            Serial.print(" ");
        }
        Serial.println();

        // Covariance prediction
        cov_prediction(A, Pp, Q, Pm);

        //Print out the covariance prediction
        Serial.print("Covariance prediction at iteration ");
        Serial.print(n);
        Serial.println(": ");
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                Serial.print(Pm[i][j]);
                Serial.print(" ");
            }
            Serial.println();
        }

        // Measurement matrix linearization
        linear_H0(xm[0], xm[2], xm[4], Hl, p, v, h, ka, m);

        kalmanGain(Hl, Pm, R, Kg); // Kalman gain

        // Measurement prediction****************************
        float a[3] = {xm[0], xm[2], 0};

        comp_H0(ka, m, v, p, a, xm[4], H0);

        for (int k = 0; k < K; k++) {
            for (int l = 0; l < 3; l++) {
                rH[3*k + l] = H0[l][k];
            }
        }

        //Print out the measurement prediction
        Serial.print("Measurement prediction at iteration ");
        Serial.print(n);
        Serial.println(": ");
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < K; j++) {
                Serial.print(rH[3*j + i]);
                Serial.print(" ");
            }
            Serial.println();
        }

        state_update(xm, Kg, r, rH, xp); // State update

        cov_update(Kg, Hl, Pm, Pp); // Covariance update

        // Output state estimate
        Serial.print("State estimate at iteration ");
        Serial.print(n);
        Serial.print(": ");
        for (int i = 0; i < 6; i++) {
            Serial.print(xp[i]);
            Serial.print(" ");
        }
        Serial.println();
    }
}

void loop() {
    // EKF updates in setup(), so no code needed here
}
