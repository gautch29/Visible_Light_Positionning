#include <math.h>
#include "mathK.h"
#include <MatrixMath.h>


// Rotation matrix
void rot_mat(float theta, float res[3][3]) {
    res[0][0] = cos(theta);
    res[0][1] = -sin(theta);
    res[0][2] = 0;
    res[1][0] = sin(theta);
    res[1][1] = cos(theta);
    res[1][2] = 0;
    res[2][0] = 0;
    res[2][1] = 0;
    res[2][2] = 1;
}
// Vector dot product
float dot_product(float a[3], float b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
// Vector norm
float norm(float a[3]) {
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}
// Cosine of angle between two vectors
float cos_ang(float a[3], float b[3]) {
    return dot_product(a, b) / (norm(a) * norm(b));
}

// Matrix-vector multiplication
void matrix_vector_mult(float mat[3][3], float vec[3], float res[3]) {
    for (int i = 0; i < 3; i++) {
        res[i] = 0;
        for (int j = 0; j < 3; j++) {
            res[i] += mat[i][j] * vec[j];
        }
    }
}

// Sign function
int sign(float value) {
    if (value > 0) {
        return 1;
    } else if (value < 0) {
        return -1;
    } else {
        return 0;
    }
}

// Kalman gain
void kalmanGain(float Hl[3*K][6], float Pm[6][6], float R[3*K][3*K], float res[6][3*K]){
    //Transposition of Hl to HlT
    float HlT[6][3*K];
    for(int i = 0;i < 3*K;i++){
        for(int j = 0;j < 6;j++){
            HlT[j][i] = Hl[i][j];
        }
    }

    //Matrix multiplication of Pm and HlT
    float Kg_tmp[6][3*K];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 3*K;j++){
            Kg_tmp[i][j] = 0;
            for(int k = 0;k < 6;k++){
                Kg_tmp[i][j] += Pm[i][k] * HlT[k][j];
            }
        }
    }

    //Matrix multiplication of Hl and Kg_tmp
    float Kg_tmp2[6][6];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 6;j++){
            Kg_tmp2[i][j] = 0;
            for(int k = 0;k < 3*K;k++){
                Kg_tmp2[i][j] += Hl[k][i] * Kg_tmp[j][k];
            }
        }
    }

    //Matrix addition of Kg_tmp2 and R
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 6;j++){
            Kg_tmp2[i][j] += R[i][j];
        }
    }

    //Matrix inversion of Kg_tmp2
    Matrix.Invert((double*)Kg_tmp2, 6);

    //Matrix multiplication of Kg_tmp2 and HlT
    float Kg_tmp3[6][3*K];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 3*K;j++){
            Kg_tmp3[i][j] = 0;
            for(int k = 0;k < 6;k++){
                Kg_tmp3[i][j] += Kg_tmp2[i][k] * HlT[k][j];
            }
        }
    }

    //Matrix multiplication of Kg_tmp3 and Pm
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 3*K;j++){
            res[i][j] = 0;
            for(int k = 0;k < 6;k++){
                res[i][j] += Kg_tmp3[i][k] * Pm[k][j];
            }
        }
    }
}

// Linearize H0 matrix
void linear_H0(float x, float y, float theta, float Hl[3*K][6], float p[3][K], float v[3][3], float h, float ka, float m) {
    int idx;
    for (int k = 0; k < K; k++) {
        for (int l = 0; l < 3; l++) {
            idx = 3*k + l;
            Hl[idx][0] = drdx(k, l, x, y, theta, p, v, h, ka, m);
            Hl[idx][1] = 0;
            Hl[idx][2] = drdy(k, l, x, y, theta, p, v, h, ka, m);
            Hl[idx][3] = 0;
            Hl[idx][4] = drdtheta(k, l, x, y, theta, p, v, h, ka, m);
            Hl[idx][5] = 0;

        }
    }
}

// Partial derivatives of range measurements
float drdx(int k, int l, float x, float y, float theta, float p[3][K], float v[3][3], float h, float ka, float m) {
    float xk = p[0][k], yk = p[1][k];
    float st = sin(theta), ct = cos(theta);
    float A = v[0][l] * ct - v[1][l] * st;
    float B = v[0][l] * st + v[1][l] * ct;
    float C = pow(xk - x, 2) + pow(yk - y, 2) + pow(h, 2);
    float D = A * (xk - x) + B * (yk - y) + h * v[2][l];
    if (is_out_fov(k, l, x, y, p, v, theta)) {
        return 0;
    } else {
        return -ka * (m + 1) * pow(h, m) * sign(D) * A / pow(C, (m + 3) / 2) + ka * (m + 1) * (m + 3) * pow(h, m) * fabs(D) * (xk - x) / pow(C, (m + 5) / 2);
    }


}

// Partial derivatives of range measurements
float drdy(int k, int l, float x, float y, float theta, float p[3][K], float v[3][3], float h, float ka, float m) {
    float xk = p[0][k], yk = p[1][k];
    float st = sin(theta), ct = cos(theta);
    float A = v[0][l] * ct - v[1][l] * st;
    float B = v[0][l] * st + v[1][l] * ct;
    float C = pow(xk - x, 2) + pow(yk - y, 2) + pow(h, 2);
    float D = A * (xk - x) + B * (yk - y) + h * v[2][l];
    if (is_out_fov(k, l, x, y, p, v, theta)) {
        return 0;
    } else {
        return -ka * (m + 1) * pow(h, m) * sign(D) * B / pow(C, (m + 3) / 2)
            + ka * (m + 1) * (m + 3) * pow(h, m) * fabs(D) * (yk - y) / pow(C, (m + 5) / 2);
    }
}
// Partial derivatives of range measurements
float drdtheta(int k, int l, float x, float y, float theta, float p[3][K], float v[3][3], float h, float ka, float m) {
    float xk = p[0][k], yk = p[1][k];
    float st = sin(theta), ct = cos(theta);
    float A = v[0][l] * ct - v[1][l] * st;
    float B = v[0][l] * st + v[1][l] * ct;
    float C = pow(xk - x, 2) + pow(yk - y, 2) + pow(h, 2);
    float D = A * (xk - x) + B * (yk - y) + h * v[2][l];
    if (is_out_fov(k, l, x, y, p, v, theta)) {
        return 0;
    } else {
        return ka * (m + 1) * (m + 3) * pow(h, m) * fabs(D) * ((xk - x) * (-v[0][l] * st - v[1][l] * ct)
            + (yk - y) * (v[0][l] * ct - v[1][l] * st)) / pow(C, (m + 5) / 2);
    }
}

// Check if reference point is outside field of view
bool is_out_fov(int k, int l, float x, float y, float p[3][K], float v[3][3], float theta) {
    float u[3], cosPsi;
    float rot[3][3];
    float tmp_v[3] = {v[0][l], v[1][l], v[2][l]};

    rot_mat(theta, rot);
    matrix_vector_mult(rot, tmp_v, u);

    float ref[3] = {p[0][k] - x, p[1][k] - y, p[2][k]};
    cosPsi = cos_ang(ref, u);
    return cosPsi < cos(FOV);
}

// Compute H0 matrix
void comp_H0(float ka, float m, float v[3][3], float p[3][K], float p0[3], float theta, float H0[3][K]) {
    float u[3], tmp[3], tmpn[3], d, cosPhi, cosPsi;
    float rot[3][3];
    rot_mat(theta, rot);

    for (int k = 0; k < K; k++) {
        for (int l = 0; l < 3; l++) {
            float tmp_v[3] = {v[0][l], v[1][l], v[2][l]};
            matrix_vector_mult(rot, tmp_v, u);

            tmp[0] = p0[0] - p[0][k];
            tmp[1] = p0[1] - p[1][k];
            tmp[2] = p0[2] - p[2][k];
            d = norm(tmp);
            float ref[3] = {0, 0, -1};
            cosPhi = cos_ang(tmp, ref);

            tmpn[0] = -tmp[0];
            tmpn[1] = -tmp[1];
            tmpn[2] = -tmp[2];

            cosPsi = cos_ang(tmpn, u);

            if (cosPsi >= cos(FOV)) {
                H0[l][k] = fabs(ka * (m + 1) / pow(d, 2) * pow(cosPhi, m) * cosPsi);
            } else {
                H0[l][k] = 0;
            }
        }
    }
}


// Random initial position and orientation
void rand_p0th0(float &x, float &y, float &theta) {
    x = A_BOUND * (2 * (random(1000)/1000.0) - 1);
    y = A_BOUND * (2 * (random(1000)/1000.0) - 1);
    theta = PI * (2 * (random(1000)/1000.0) - 1);
}

//State prediction
void state_prediction(float A[6][6], float xp[6], float xm[6]){
    for (int i = 0; i < 6; i++) {
        xm[i] = 0;
        for (int j = 0; j < 6; j++) {
            xm[i] += A[i][j] * xp[j];
        }
    }
}

//Covariance prediction
void cov_prediction(float A[6][6], float Pp[6][6], float Q[6][6], float Pm[6][6]){
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            Pm[i][j] = Q[i][j];
            for (int k = 0; k < 6; k++) {
                Pm[i][j] += A[i][k] * Pp[k][j];
            }
        }
    }
}

//State update
void state_update(float xm[6], float Kg[6][3*K], float r[3*K], float rH[3*K], float xp[6]){
    for (int i = 0; i < 6; i++) {
        xp[i] = xm[i];
        for (int j = 0; j < 3*K; j++) {
            xp[i] += Kg[i][j] * (r[j] - rH[j]);
        }
    }
}

//Covaariance update
void cov_update(float Kg[6][3*K], float Hl[3*K][6], float Pm[6][6], float Pp[6][6]){
    float eye6[6][6] = { {1, 0, 0, 0, 0, 0},
                        {0, 1, 0, 0, 0, 0},
                        {0, 0, 1, 0, 0, 0},
                        {0, 0, 0, 1, 0, 0},
                        {0, 0, 0, 0, 1, 0},
                        {0, 0, 0, 0, 0, 1}};

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            Pp[i][j] = 0;
            for (int k = 0; k < 6; k++) {
                Pp[i][j] += Kg[i][k] * Hl[k][j];
            }
            Pp[i][j] = eye6[i][j] - Pp[i][j];
            for (int k = 0; k < 6; k++) {
                Pp[i][j] += Pp[i][k] * Pm[k][j];
            }
        }
    }
}