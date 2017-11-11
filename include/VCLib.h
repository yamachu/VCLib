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
void SPTK_mlsadf(
    #ifdef DOUBLE
    double *wavform, int wavform_length, double *mcep, int mcep_length,
    #else
    float *wavform, int wavform_length, float *mcep, int mcep_length,
    #endif
    int order, double alpha, int frame_period, int i_period, int pade,
    int is_tranpose, int is_invrese, int is_coef_b, int is_without_gain, 
    #ifdef DOUBLE
    double *y
    #else
    float *y
    #endif
);
void DifferentialMelCepstrumCompensation(
    #ifdef DOUBLE
    double *rawform, int rawform_length,
    double *sp, int sp_length,
    double *diff_mcep, int diff_mcep_length,
    #else
    float *rawform, int rawform_length,
    float *sp, int sp_length,
    float *diff_mcep, int diff_mcep_length,
    #endif
    double alpha, double gamma,
    int mcep_order,
    int fft_length,
    int itype,
    int otype,
    int min_iter,
    int max_iter,
    int recursions,
    double eps,
    double end_cond,
    int etype,
    double min_det,
    int frame_period,
    int interpolate_period,
    #ifdef DOUBLE
    double *out
    #else
    float *out
    #endif
);