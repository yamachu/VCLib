void SPTK_mgcep(
    #ifdef DOUBLE
    double *spectrum,
    #else
    float *spectrum,
    #endif
    int sp_length, double alpha, double gamma,
    int order, int fft_length, int itype, int otype,
    int min_iter, int max_iter, int recursions,
    double eps, double end_cond, int etype, double min_det,
    #ifdef DOUBLE
    double *mcep
    #else
    float *mcep
    #endif
);
int SPTK_mlsadf(
    #ifdef DOUBLE
    double *wavform, int wavform_length, double *mcep, int mcep_length,
    #else
    float *wavform, int wavform_length, float *mcep, int mcep_length,
    #endif
    int order, double alpha, int frame_period, int i_period, int pade,
    int is_transpose, int is_invrese, int is_coef_b, int is_without_gain, 
    #ifdef DOUBLE
    double *y
    #else
    float *y
    #endif
);
void GetUserMcep(
    #ifdef DOUBLE
    double *sp, int length, double *result
    #else
    float *sp, int length, float *result
    #endif
);
int GetCompensationWavForm(
    #ifdef DOUBLE
    double *x,
    int x_length,
    double *userMcep,
    double *targetMcep,
    int mcep_length,
    double *y
    #else
    float *x,
    int x_length,
    float *userMcep,
    float *targetMcep,
    int mcep_length,
    float *y
    #endif
);
void Standardization1DArray(float *source, int length, int dim, float *result);
void UnStandardization1DArray(float *source, int length, int dim, float *means, float *sds, float *result);
void VarianceCompensation(float *source, int length, int dim, float *coef, float *result);
int GetInversedWavForm(
    #ifdef DOUBLE
    double *x,
    int x_length,
    double *userMcep,
    int mcep_length,
    double *y
    #else
    float *x,
    int x_length,
    float *userMcep,
    int mcep_length,
    float *y
    #endif
);
int GetFilteredWavForm(
    #ifdef DOUBLE
    double *x,
    int x_length,
    double *targetMcep,
    int mcep_length,
    double *y
    #else
    float *x,
    int x_length,
    float *targetMcep,
    int mcep_length,
    float *y
    #endif
);
void SPTK_mgc2sp(
    #ifdef DOUBLE
    double *fp /* mcep */,
    #else
    float *fp /* mcep */,
    #endif
    int mcep_length,
    double alpha,
    double gamma,
    int m,
    int norm,
    int mulg,
    int l /* FFT length */,
    int phase,
    int otype,
    #ifdef DOUBLE
    double *output
    #else
    float *output
    #endif
);