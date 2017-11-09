#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(WIN32)
#include "SPTK.h"
#else
#include <SPTK.h>
#endif

void SPTK_mgcep(double *spectrum, int f0_length /* input */,
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
    double *mcep /* output */)
{
    int ilng, n, i, j;

    n = recursions;

    if (recursions == -1)
        n = fft_length - 1;

    if (itype == 0)
        ilng = fft_length;
    else
        ilng = fft_length / 2 + 1;

    for (i = 0; i < f0_length; i++)
    {
        mgcep(&spectrum[i * ilng], fft_length, &mcep[i * (order + 1)], order, alpha, gamma, n, min_iter, max_iter, end_cond, etype, eps, min_det, itype);

        if (otype == 0 || otype == 1 || otype == 2 || otype == 4)
            ignorm(&mcep[i * (order + 1)], &mcep[i * (order + 1)], order, gamma);

        if (otype == 0 || otype == 2 || otype == 4)
            if (alpha != 0.0)
                b2mc(&mcep[i * (order + 1)], &mcep[i * (order + 1)], order, alpha);

        if (otype == 2 || otype ==4)
            gnorm(&mcep[i * (order + 1)], &mcep[i * (order + 1)], order, gamma);

        if (otype == 4 || otype == 5)
            for (j = order; j >= 1; j--)
                mcep[i * (order + 1) + j] *= gamma;
    }
}

void SPTK_mlsadf(double *wavform, int wavform_length, /* input */
    double *mcep, int f0_length, /* input */
    int order, /* m */
    double alpha, /* a */
    int frame_period, /* p */
    int i_period, /* i */ 
    int pade, /* P */
    int is_tranpose, /* t */
    int is_invrese, /* v */
    int is_coef_b, /* b */
    int is_without_gain, /* g */
    double *y /* output */)
{
    int i, j, k, y_count;
    double *c, *inc, *cc, *d, x;

    c = dgetmem(3 * (order + 1) + 3 * (pade + 1) + pade * (order + 2));
    cc = c + order + 1;
    inc = cc + order + 1;
    d = inc + order + 1;

    movem(mcep, c, sizeof(*mcep), order + 1);

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
        movem(&mcep[i * (order + 1)], cc, sizeof(*mcep), order + 1);

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

            y[y_count++] = x;

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

