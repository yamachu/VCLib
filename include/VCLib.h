void SPTK_mgcep(double *spectrum, int sp_length, double alpha, double gamma,
    int order, int fft_length, int itype, int otype,
    int min_iter, int max_iter, int recursions,
    double eps, double end_cond, int etype, double min_det, double *mcep);
void SPTK_mlsadf(double *wavform, int wavform_length, double *mcep, int mcep_length,
    int order, double alpha, int frame_period, int i_period, int pade,
    int is_tranpose, int is_invrese, int is_coef_b, int is_without_gain, double *y);