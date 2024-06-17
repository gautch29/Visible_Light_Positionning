#define K 2 // Number of reference points
#define FOV 1.5708 // 90 degrees in radians
#define A_BOUND 4 // RX pos in [-A, A, -A, A]

// Rotation matrix
void rot_mat(double theta, double res[3][3]);

// Vector dot product
double dot_product(double a[3], double b[3]);

// Vector norm
double norm(double a[3]);

// Cosine of angle between two vectors
double cos_ang(double a[3], double b[3]);

// Matrix-vector multiplication
void matrix_vector_mult(double mat[3][3], double vec[3], double res[3]);

// Sign function
int sign(double value);

// Kalman gain calculation
void kalmanGain(double Hl[3*K][6], double Pm[6][6], double R[3*K][3*K], double res[6][3*K]);

// Linearization of the measurement model
void linear_H0(double x, double y, double theta, double Hl[3*K][6], double p[3][K], double v[3][3], double h, double ka, double m);

// Derivative of the measurement model with respect to x
double drdx(int k, int l, double x, double y, double theta, double p[3][K], double v[3][3], double h, double ka, double m);

// Derivative of the measurement model with respect to y
double drdy(int k, int l, double x, double y, double theta, double p[3][K], double v[3][3], double h, double ka, double m);

// Derivative of the measurement model with respect to theta
double drdtheta(int k, int l, double x, double y, double theta, double p[3][K], double v[3][3], double h, double ka, double m);

// Check if a reference point is out of the field of view
bool is_out_fov(int k, int l, double x, double y, double p[3][K], double v[3][3], double theta);

// Measurement model
void comp_H0(double ka, double m, double v[3][3], double p[3][K], double p0[3], double theta, double H0[3][K]);

// State prediction
void state_prediction(double A[6][6], double xp[6], double xm[6]);

// Covariance prediction
void cov_prediction(double A[6][6], double Pp[6][6], double Q[6][6], double Pm[6][6]);

// State update
void state_update(double xm[6], double Kg[6][3*K], double r[3*K], double rH[3*K], double xp[6]);

// Covariance update
void cov_update(double Kg[6][3*K], double Hl[3*K][6], double Pm[6][6], double Pp[6][6], double R[3*K][3*K]);

// Print matrix
void printMatrix(double* mat, int rows, int cols);

// Kalman filter
void kalmanFilter(double A[6][6], double xp[6], double Pp[6][6], double Q[6][6], double R[3*K][3*K], double p[3][K], double v[3][3], double h, double ka, double m, double xm[6], double Pm[6][6], double r[3*K], double H0[3][K]);

// Print Kalman results
void printKalmanResults(double xp[6]);

// Compute r0
void r0Comp(double P, double H0[3][K], double r0[3*K]);

// Generate beam vectors
void beamVectors(double v[3][3], double al[3], double be);

class kalmanFilter{
    private:
        double T = 0.01; // Time step
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
        double R[3*K][3*K];
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
        double p[3][K];
        double v[3][3];
        double h;
        double ka;
        double m;
        double xm[6];
        double Pm[6][6];
        double r[3*K];
        double Hl[3*K][6];
        double Kg[6][3*K];
        double rH[3*K];


    public:
        double xp[6] = {0, 0, 0, 0, 0, 0};
        void compute(double h, double ka, double m, double p[3][K], double v[3][3]);
};