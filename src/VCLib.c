#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(WIN32)
#include "SPTK.h"
#else
#include <SPTK.h>
#endif

#include "VCLib.h"

void SPTK_mgcep(
#ifdef DOUBLE
    double *spectrum, int sp_length /* input */,
#else
    float *spectrum, int sp_length /* input */,
#endif
    double alpha /* a */,
    double gamma /* g */,
    int order /* m */,
    int fft_length /* l */,
    int itype /* q */,
    int otype /* o */,
    int min_iter /* i */,
    int max_iter /* j */,
    int recursions /* p */,
    double eps /* e, E */,
    double end_cond /* d */,
    int etype /* e, E */,
    double min_det /* f */,
#ifdef DOUBLE
    double *mcep /* output */)
#else
    float *mcep /* output */)
#endif
{
    int ilng, n, i, j, total_frame;
    double *x, *b;

    n = recursions;

    if (recursions == -1)
        n = fft_length - 1;

    if (itype == 0)
        ilng = fft_length;
    else
        ilng = fft_length / 2 + 1;

    total_frame = sp_length / (fft_length / 2 + 1);
    x = dgetmem(fft_length + order + order + 2);
    b = x + fft_length;

    for (i = 0; i < total_frame; i++)
    {
        for (j = 0; j < ilng; j++)
        {
            x[j] = spectrum[i * ilng + j];
        }

        mgcep(x, fft_length, b, order, alpha, gamma, n, min_iter, max_iter, end_cond, etype, eps, min_det, itype);
        
        if (otype == 0 || otype == 1 || otype == 2 || otype == 4)
            ignorm(b, b, order, gamma);

        if (otype == 0 || otype == 2 || otype == 4)
            if (alpha != 0.0)
                b2mc(b, b, order, alpha);

        if (otype == 2 || otype ==4)
            gnorm(b, b, order, gamma);

        if (otype == 4 || otype == 5)
            for (j = order; j >= 1; j--)
                b[j] *= gamma;

        for (j = 0; j < order + 1; j++)
            #ifdef DOUBLE
            mcep[i * (order + 1) + j] = b[j];
            #else
            mcep[i * (order + 1) + j] = (float)b[j];
            #endif
    }

    free(x);
}

void SPTK_mlsadf(
    #ifdef DOUBLE
    double *wavform, int wavform_length, /* input */
    double *mcep, int mcep_length, /* input */
    #else
    float *wavform, int wavform_length, /* input */
    float *mcep, int mcep_length, /* input */
    #endif
    int order, /* m */
    double alpha, /* a */
    int frame_period, /* p */
    int i_period, /* i */ 
    int pade, /* P */
    int is_tranpose, /* t */
    int is_invrese, /* v */
    int is_coef_b, /* b */
    int is_without_gain, /* k */
    #ifdef DOUBLE
    double *y /* output */)
    #else
    float *y /* output */)
    #endif
{
    int i, j, k, y_count, f0_length;
    double *c, *inc, *cc, *d, x;

    f0_length = mcep_length / (order + 1);

    c = dgetmem(3 * (order + 1) + 3 * (pade + 1) + pade * (order + 2));
    cc = c + order + 1;
    inc = cc + order + 1;
    d = inc + order + 1;

    for (i = 0; i < order + 1; i++)
        c[i] = mcep[i];

    if (is_coef_b == 0)
        mc2b(c, c, order, alpha);

    if (is_invrese)
    {
        for (i = 0; i <= order; i++)
            c[i] *= -1;
        c[0] = is_without_gain ? 0 : c[0];
    }

    for (i = 1, y_count = 0; i < f0_length; i++)
    {
        for (j = 0; j < order + 1; j++)
            cc[j] = mcep[i * (order + 1) + j];

        if (is_coef_b == 0)
            mc2b(cc, cc, order, alpha);

        if (is_invrese)
        {
            for (j = 0; j <= order; j++)
                cc[i] *= -1;
            cc[0] = is_without_gain ? 0 : cc[0];
        }

        for (j = 0; j <= order; j++)
            inc[j] = (cc[j] - c[j]) * (double)i_period / (double)frame_period;

        for (k = frame_period, j = (i_period + 1) / 2; k--;)
        {
            x = wavform[y_count];

            if (is_without_gain == 0)
                x *= exp(c[0]);
            if (is_tranpose)
                x = mlsadft(x, c, order, alpha, pade, d);
            else
                x = mlsadf(x, c, order, alpha, pade, d);

            #ifdef DOUBLE
            y[y_count++] = x;
            #else
            y[y_count++] = (float)x;
            #endif

            if (!--j)
            {
                for (j = 0; j <= order; j++)
                    c[j] += inc[j];
                j = i_period;
            }
        }

        movem(cc, c, sizeof(*cc), order + 1);
    }
}

