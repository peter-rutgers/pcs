#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include "compute.h"

#ifdef GEN_PICTURES
static void do_draw(const struct parameters *p,
             size_t key, size_t h, size_t w,
             double (*restrict g)[h][w])
{
    begin_picture(key, w-2, h-2, p->io_tmin, p->io_tmax);
    size_t i, j;
    for (i = 1; i < h-1; ++i)
        for (j = 1; j < w-1; ++j)
            draw_point(j-1, i-1, (*g)[i][j]);
    end_picture();
}
#endif

static void do_copy(size_t h, size_t w,
             double (*restrict g)[h][w])
{
    size_t i;

    /* copy left and right column to opposite border */
    for (i = 0; i < h; ++i) {
        (*g)[i][w-1] = (*g)[i][1];
        (*g)[i][0] = (*g)[i][w-2];
    }
}

static void fill_report(const struct parameters *p, struct results *r,
                        size_t h, size_t w, 
                        double (*a)[h][w],
                        double maxdiff, 
                        double iter,
                        struct timeval *before)
{
    /* compute min/max/avg */
    double tmin = INFINITY, tmax = -INFINITY;
    double sum = 0.0;
    struct timeval after;

    for (size_t i = 1; i < h - 1; ++i)
        for (size_t j = 1; j < w - 1; ++j) 
        {
            double v = (*a)[i][j];
            sum += v;
            if (tmin > v) tmin = v;
            if (tmax < v) tmax = v;
        }

    r->niter = iter;
    r->maxdiff = maxdiff;
    r->tmin = tmin;
    r->tmax = tmax;
    r->tavg = sum / (p->N * p->M);

    gettimeofday(&after, NULL);
    r->time = (double)(after.tv_sec - before->tv_sec) + 
        (double)(after.tv_usec - before->tv_usec) / 1e6;
}

void do_compute(const struct parameters* p, struct results *r)
{

    size_t i, j;

    /* alias input parameters */
    const double (*restrict tinit)[p->N][p->M] = (const double (*)[p->N][p->M])p->tinit;
    const double (*restrict cinit)[p->N][p->M] = (const double (*)[p->N][p->M])p->conductivity;

    /* allocate grid data */
    const size_t h = p->N + 2;
    const size_t w = p->M + 2;
    double (*restrict g1)[h][w] = malloc(h * w * sizeof(double));
    double (*restrict g2)[h][w] = malloc(h * w * sizeof(double));

    /* allocate halo for conductivities */
    double (*restrict c)[h][w] = malloc(h * w * sizeof(double));

    struct timeval before;
    gettimeofday(&before, NULL);

    static const double c_cdir = 0.25 * M_SQRT2 / (M_SQRT2 + 1.0);
    static const double c_cdiag = 0.25 / (M_SQRT2 + 1.0);

    /* set initial temperatures and conductivities */
    for (i = 1; i < h - 1; ++i)
        for (j = 1; j < w - 1; ++j) 
        {
            (*g1)[i][j] = (*tinit)[i-1][j-1];
            (*c)[i][j] = (*cinit)[i-1][j-1];
        }

    /* smear outermost row to border */
    for (j = 1; j < w-1; ++j) {
        (*g1)[0][j] = (*g2)[0][j] = (*g1)[1][j];
        (*g1)[h-1][j] = (*g2)[h-1][j] = (*g1)[h-2][j];
    }



    /* compute */
    size_t iter;
    double maxdiff = 0.0;
    double (*restrict src)[h][w] = g2;
    double (*restrict dst)[h][w] = g1;

    for (iter = 1; iter <= p->maxiter; ++iter)
    {
#ifdef GEN_PICTURES
        do_draw(p, iter, h, w, src);
#endif
        /* swap source and destination */
        { void *tmp = src; src = dst; dst = tmp; }

        /* initialize halo on source */
        do_copy(h, w, src);

        /* compute */
        maxdiff = 0.0;
        for (i = 1; i < h - 1; ++i)
            for (j = 1; j < w - 1; ++j)
            {
                double w = (*c)[i][j];
                double restw = 1.0 - w;

                (*dst)[i][j] = w * (*src)[i][j] + 

                    ((*src)[i+1][j  ] + (*src)[i-1][j  ] + 
                     (*src)[i  ][j+1] + (*src)[i  ][j-1]) * (restw * c_cdir) +

                    ((*src)[i-1][j-1] + (*src)[i-1][j+1] + 
                     (*src)[i+1][j-1] + (*src)[i+1][j+1]) * (restw * c_cdiag);

                double diff = fabs((*src)[i][j] - (*dst)[i][j]);
                if (diff > maxdiff)
                    maxdiff = diff;
            }
     
        /* check for convergence */
        if (maxdiff < p->threshold) 
            { ++iter; break; }

        /* conditional reporting */
        if (iter % p->period == 0) {
            fill_report(p, r, h, w, dst, maxdiff, iter, &before);
            report_results(p, r);
        }
    }

    /* report at end in all cases */
    fill_report(p, r, h, w, dst, maxdiff, iter-1, &before);

    free(c);
    free(g2);
    free(g1);
}
