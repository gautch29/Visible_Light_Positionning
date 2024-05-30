#define K 9 // Number of reference points
#define FOV 1.5708 // 90 degrees in radians
#define A_BOUND 4 // RX pos in [-A, A, -A, A]

// Rotation matrix
void rot_mat(float theta, float res[3][3]);

// Vector dot product
float dot_product(float a[3], float b[3]);

// Vector norm
float norm(float a[3]);

// Cosine of angle between two vectors
float cos_ang(float a[3], float b[3]);

// Matrix-vector multiplication
void matrix_vector_mult(float mat[3][3], float vec[3], float res[3]);

// Sign function
int sign(float value);

// Kalman gain calculation
void kalmanGain(float Hl[3*K][6], float Pm[6][6], float R[3*K][3*K], float res[6][3*K]);

// Linearization of the measurement model
void linear_H0(float x, float y, float theta, float Hl[3*K][6], float p[3][K], float v[3][3], float h, float ka, float m);

// Derivative of the measurement model with respect to x
float drdx(int k, int l, float x, float y, float theta, float p[3][K], float v[3][3], float h, float ka, float m);

// Derivative of the measurement model with respect to y
float drdy(int k, int l, float x, float y, float theta, float p[3][K], float v[3][3], float h, float ka, float m);

// Derivative of the measurement model with respect to theta
float drdtheta(int k, int l, float x, float y, float theta, float p[3][K], float v[3][3], float h, float ka, float m);

// Check if a reference point is out of the field of view
bool is_out_fov(int k, int l, float x, float y, float p[3][K], float v[3][3], float theta);

// Measurement model
void comp_H0(float ka, float m, float v[3][3], float p[3][K], float p0[3], float theta, float H0[3][K]);

// State prediction
void state_prediction(float A[6][6], float xp[6], float xm[6]);

// Covariance prediction
void cov_prediction(float A[6][6], float Pp[6][6], float Q[6][6], float Pm[6][6]);

// State update
void state_update(float xm[6], float Kg[6][3*K], float r[3*K], float rH[3*K], float xp[6]);

// Covariance update
void cov_update(float Kg[6][3*K], float Hl[3*K][6], float Pm[6][6], float Pp[6][6]);