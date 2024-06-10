#include <math.h>
#include "mathK.h"
#include "MatrixMath.h"


// Rotation matrix
void rot_mat(double theta, double res[3][3]) {
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
double dot_product(double a[3], double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
// Vector norm
double norm(double a[3]) {
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}
// Cosine of angle between two vectors
double cos_ang(double a[3], double b[3]) {
    return dot_product(a, b) / (norm(a) * norm(b));
}

// Matrix-vector multiplication
void matrix_vector_mult(double mat[3][3], double vec[3], double res[3]) {
    for (int i = 0; i < 3; i++) {
        res[i] = 0;
        for (int j = 0; j < 3; j++) {
            res[i] += mat[i][j] * vec[j];
        }
    }
}

// Sign function
int sign(double value) {
    if (value > 0) {
        return 1;
    } else if (value < 0) {
        return -1;
    } else {
        return 0;
    }
}

// Kalman gain
void kalmanGain(double Hl[3*K][6], double Pm[6][6], double R[3*K][3*K], double res[6][3*K]){
    //Transposition of Hl to HlT
    double HlT[6][3*K];
    for(int i = 0;i < 3*K;i++){
        for(int j = 0;j < 6;j++){
            HlT[j][i] = Hl[i][j];
        }
    }

    //Matrix multiplication of Pm and HlT
    double Kg_tmp[6][3*K];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 3*K;j++){
            Kg_tmp[i][j] = 0;
            for(int k = 0;k < 6;k++){
                Kg_tmp[i][j] += Pm[i][k] * HlT[k][j];
            }
        }
    }

    //Matrix multiplication of Hl and Pm
    double Kg_tmp2[3*K][6];
    for(int i = 0;i < 3*K;i++){
        for(int j = 0;j < 6;j++){
            Kg_tmp2[i][j] = 0;
            for(int k = 0;k < 6;k++){
                Kg_tmp2[i][j] += Hl[i][k] * Pm[k][j];
            }
        }
    }

    //Matrix multiplication of Kg_tmp2 and HlT
    double Kg_tmp3[3*K][3*K];
    for(int i = 0;i < 3*K;i++){
        for(int j = 0;j < 3*K;j++){
            Kg_tmp3[i][j] = 0;
            for(int k = 0;k < 6;k++){
                Kg_tmp3[i][j] += Kg_tmp2[i][k] * HlT[k][j];
            }
        }
    }

    //Matrix addition of R and Kg_tmp3
    double Kg_tmp4[3*K][3*K];
    for(int i = 0;i < 3*K;i++){
        for(int j = 0;j < 3*K;j++){
            Kg_tmp4[i][j] = R[i][j] + Kg_tmp3[i][j];
        }
    }

    //Matrix inversion of Kg_tmp4
    Matrix.Invert((double*)Kg_tmp4, 3*K);

    //Matrix multiplication of Kg_tmp and Kg_tmp4
    double kg[6][3*K];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 3*K;j++){
            kg[i][j] = 0;
            for(int k = 0;k < 3*K;k++){
                kg[i][j] += Kg_tmp[i][k] * Kg_tmp4[k][j];
            }
        }
    }

    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 3*K;j++){
            res[i][j] = kg[i][j];
        }
    }

}

// Linearize H0 matrix
void linear_H0(double x, double y, double theta, double Hl[3*K][6], double p[3][K], double v[3][3], double h, double ka, double m) {
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
double drdx(int k, int l, double x, double y, double theta, double p[3][K], double v[3][3], double h, double ka, double m) {
    double xk = p[0][k], yk = p[1][k];
    double st = sin(theta), ct = cos(theta);
    double A = v[0][l] * ct - v[1][l] * st;
    double B = v[0][l] * st + v[1][l] * ct;
    double C = pow(xk - x, 2) + pow(yk - y, 2) + pow(h, 2);
    double D = A * (xk - x) + B * (yk - y) + h * v[2][l];
    if (is_out_fov(k, l, x, y, p, v, theta)) {
        return 0;
    } else {
        return -ka * (m + 1) * pow(h, m) * sign(D) * A / pow(C, (m + 3) / 2) + ka * (m + 1) * (m + 3) * pow(h, m) * fabs(D) * (xk - x) / pow(C, (m + 5) / 2);
    }
}

// Partial derivatives of range measurements
double drdy(int k, int l, double x, double y, double theta, double p[3][K], double v[3][3], double h, double ka, double m) {
    double xk = p[0][k], yk = p[1][k];
    double st = sin(theta), ct = cos(theta);
    double A = v[0][l] * ct - v[1][l] * st;
    double B = v[0][l] * st + v[1][l] * ct;
    double C = pow(xk - x, 2) + pow(yk - y, 2) + pow(h, 2);
    double D = A * (xk - x) + B * (yk - y) + h * v[2][l];
    if (is_out_fov(k, l, x, y, p, v, theta)) {
        return 0;
    } else {
        return -ka * (m + 1) * pow(h, m) * sign(D) * B / pow(C, (m + 3) / 2)
            + ka * (m + 1) * (m + 3) * pow(h, m) * fabs(D) * (yk - y) / pow(C, (m + 5) / 2);
    }
}

// Partial derivatives of range measurements
double drdtheta(int k, int l, double x, double y, double theta, double p[3][K], double v[3][3], double h, double ka, double m) {
    double xk = p[0][k], yk = p[1][k];
    double st = sin(theta), ct = cos(theta);
    double A = v[0][l] * ct - v[1][l] * st;
    double B = v[0][l] * st + v[1][l] * ct;
    double C = pow(xk - x, 2) + pow(yk - y, 2) + pow(h, 2);
    double D = A * (xk - x) + B * (yk - y) + h * v[2][l];
    if (is_out_fov(k, l, x, y, p, v, theta)) {
        return 0;
    } else {
        return ka * (m + 1) * pow(h, m) * sign(D) / pow(C, (m + 3) / 2) * (-(xk-x)*(v[0][l]*st + v[1][l]*ct) + (yk -y)*(v[0][l]*ct - v[1][l]*st));
    }
}

// Check if reference point is outside field of view
bool is_out_fov(int k, int l, double x, double y, double p[3][K], double v[3][3], double theta) {
    double u[3], cosPsi;
    double rot[3][3];
    double tmp_v[3] = {v[0][l], v[1][l], v[2][l]};

    rot_mat(theta, rot);
    matrix_vector_mult(rot, tmp_v, u);

    double ref[3] = {p[0][k] - x, p[1][k] - y, p[2][k]};
    cosPsi = cos_ang(ref, u);
    return cosPsi < cos(FOV);
}

// Compute H0 matrix
void comp_H0(double ka, double m, double v[3][3], double p[3][K], double p0[3], double theta, double H0[3][K]) {
    double u[3], tmp[3], tmpn[3], d, cosPhi, cosPsi;
    double rot[3][3];
    rot_mat(theta, rot);

    for (int k = 0; k < K; k++) {
        for (int l = 0; l < 3; l++) {
            double tmp_v[3] = {v[0][l], v[1][l], v[2][l]};
            matrix_vector_mult(rot, tmp_v, u);

            tmp[0] = p0[0] - p[0][k];
            tmp[1] = p0[1] - p[1][k];
            tmp[2] = p0[2] - p[2][k];
            d = norm(tmp);
            double ref[3] = {0, 0, -1};
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
void rand_p0th0(double &x, double &y, double &theta) {
    x = A_BOUND * (2 * (random(1000)/1000.0) - 1);
    y = A_BOUND * (2 * (random(1000)/1000.0) - 1);
    theta = PI * (2 * (random(1000)/1000.0) - 1);
}

//State prediction
void state_prediction(double A[6][6], double xp[6], double xm[6]){
    for (int i = 0; i < 6; i++) {
        xm[i] = 0;
        for (int j = 0; j < 6; j++) {
            xm[i] += A[i][j] * xp[j];
        }
    }
}

//Covariance prediction
void cov_prediction(double A[6][6], double Pp[6][6], double Q[6][6], double Pm[6][6]){
    //Transposition of A to AT
    double AT[6][6];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 6;j++){
            AT[j][i] = A[i][j];
        }
    }

    //Matrix multiplication of A and Pp
    double Pm_tmp[6][6];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 6;j++){
            Pm_tmp[i][j] = 0;
            for(int k = 0;k < 6;k++){
                Pm_tmp[i][j] += A[i][k] * Pp[k][j];
            }
        }
    }

    //Matrix multiplication of Pm_tmp and AT and addition of Q
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 6;j++){
            Pm[i][j] = 0;
            for(int k = 0;k < 6;k++){
                Pm[i][j] += Pm_tmp[i][k] * AT[k][j];
            }
            Pm[i][j] += Q[i][j];
        }
    }
}

//State update
void state_update(double xm[6], double Kg[6][3*K], double r[3*K], double rH[3*K], double xp[6]){
    double rMinusRh[3*K];

    for(int i = 0;i < 3*K;i++){
        rMinusRh[i] = r[i] - rH[i];
    }

    //Kg * (r - rH)
    for(int i = 0;i < 6;i++){
        xp[i] = 0;
        for(int j = 0;j < 3*K;j++){
            xp[i] += Kg[i][j] * rMinusRh[j];
        }
    }

    //xp = xm + Kg * (r - rH)
    for(int i = 0;i < 6;i++){
        xp[i] += xm[i];
    }
}

//Covaariance update
void cov_update(double Kg[6][3*K], double Hl[3*K][6], double Pm[6][6], double Pp[6][6], double R[3*K][3*K]){
    double eye6[6][6] = { {1, 0, 0, 0, 0, 0},
                        {0, 1, 0, 0, 0, 0},
                        {0, 0, 1, 0, 0, 0},
                        {0, 0, 0, 1, 0, 0},
                        {0, 0, 0, 0, 1, 0},
                        {0, 0, 0, 0, 0, 1}};

    //Kg * Hl
    double KgHl[6][6];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 6;j++){
            KgHl[i][j] = 0;
            for(int k = 0;k < 3*K;k++){
                KgHl[i][j] += Kg[i][k] * Hl[k][j];
            }
        }
    }

    //eye6 - Kg * Hl
    double eye6MinusKgHl[6][6];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 6;j++){
            eye6MinusKgHl[i][j] = eye6[i][j] - KgHl[i][j];
        }
    }

    //eye6MinusKgHl * Pm
    double eye6MinusKgHlPm[6][6];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 6;j++){
            eye6MinusKgHlPm[i][j] = 0;
            for(int k = 0;k < 6;k++){
                eye6MinusKgHlPm[i][j] += eye6MinusKgHl[i][k] * Pm[k][j];
            }
        }
    }

    //transpose of eye6MinusKgHl
    double eye6MinusKgHlT[6][6];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 6;j++){
            eye6MinusKgHlT[j][i] = eye6MinusKgHl[i][j];
        }
    }

    //eye6minusKgHlPm * eye6MinusKgHlT
    double eye6MinusKgHlPmEye6MinusKgHlT[6][6];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 6;j++){
            eye6MinusKgHlPmEye6MinusKgHlT[i][j] = 0;
            for(int k = 0;k < 6;k++){
                eye6MinusKgHlPmEye6MinusKgHlT[i][j] += eye6MinusKgHlPm[i][k] * eye6MinusKgHlT[k][j];
            }
        }
    }

    //KgT
    double KgT[3*K][6];
    for(int i = 0;i < 3*K;i++){
        for(int j = 0;j < 6;j++){
            KgT[i][j] = Kg[j][i];
        }
    }

    //Kg * R
    double KgR[6][3*K];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 3*K;j++){
            KgR[i][j] = 0;
            for(int k = 0;k < 3*K;k++){
                KgR[i][j] += Kg[i][k] * R[k][j];
            }
        }
    }

    //Kg * R * KgT
    double KgRKgT[6][6];
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 6;j++){
            KgRKgT[i][j] = 0;
            for(int k = 0;k < 3*K;k++){
                KgRKgT[i][j] += KgR[i][k] * KgT[k][j];
            }
        }
    }

    //eye6MinusKgHlPm * eye6MinusKgHlT + KgRKgT
    for(int i = 0;i < 6;i++){
        for(int j = 0;j < 6;j++){
            Pp[i][j] = eye6MinusKgHlPmEye6MinusKgHlT[i][j] + KgRKgT[i][j];
        }
    }
}

void printMatrix(double* mat, int rows, int cols){
    for(int i = 0;i < rows;i++){
        for(int j = 0;j < cols;j++){
            Serial.print(mat[i*cols + j]);
            Serial.print(" ");
        }
        Serial.println();
    }
    Serial.println();
}