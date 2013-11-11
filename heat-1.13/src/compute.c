#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include "compute.h"
#include "fail.h"

void calculate_results(const struct parameters* p, struct results *result, double *tnew, int i, double maxdiff, double t1) {
    result->tmin = tnew[0];
    result->tmax = tnew[0];
    result->niter = i;
    result->maxdiff = maxdiff;
    double tsum = 0;
    double temp, t2;
    struct timeval now;

    for (int r = 0; r < p->N; r++) {
         for (int c = 0; c < p->M; c++) {
              temp = tnew[r * p->M + c];
              if (temp < result->tmin) {result->tmin = temp;}
              if (temp > result->tmax) {result->tmax = temp;}
              tsum += temp;
         }
    }

    result->tavg = tsum / (p->N * p->M);

    gettimeofday(&now,NULL);
    t2 = now.tv_sec + now.tv_usec / 1e6;
    result->time = t2 - t1;

    report_results(p, result);
}

void do_compute(const struct parameters* p, struct results *result)
{
    double *told, *tnew, *tmp;
    double *upper, *row, *lower;
    double cond, diff, t1;
    int i, r, c, left, right;

    struct timeval now;

    gettimeofday(&now,NULL);
    t1 = now.tv_sec + now.tv_usec / 1e6;

    if (!(told = calloc(p->N * p->M, sizeof(double)))) die("calloc");
    if (!(tnew = calloc(p->N * p->M, sizeof(double)))) die("calloc");

    for (r = 0; r < p->N; r++) {
         for (c = 0; c < p->M; c++) {
              tnew[r * p->M + c] = p->tinit[r * p->M + c];
         }
    }

    double direct_weight = sqrt(2) / (sqrt(2) + 1) / 4;
    double diagonal_weight = 1 / (sqrt(2) + 1) / 4;
    double direct_neighbours, diagonal_neighbours;

    double maxdiff = 0; // highest difference for this iteration, to compare with convergence threshold

    for (i = 0; i < p->maxiter; i++) {

        if (i > 0 && maxdiff < p->threshold) {break;}   // Check if convergence threshold is reached

        // swap told and tnew matrices
        tmp = told;
        told = tnew;
        tnew = tmp;

        maxdiff = 0;

        for (r = 0; r < p->N; r++) {

            upper = (double *) told + (r - 1) * p->M;
            row   = (double *) told +  r      * p->M;
            lower = (double *) told + (r + 1) * p->M;

            // Use fixed row containing initial values for boundary of cylinder
            if (r == 0) {
                upper = (double *) p->tinit;
            } else if (r == p->N - 1) {
                lower = (double *) p->tinit + r * p->M;
            }

            for (c = 0; c < p->M; c++) {

                // Calculate new temperature value for each point
                left  = (c - 1 + p->M) % p->M;
                right = (c + 1) % p->M;
                cond = p->conductivity[r * p->M + c];

                diagonal_neighbours = upper[left] + upper[right] + lower[left] + lower[right];
                direct_neighbours = upper[c] + lower[c] + row[left] + row[right];

                tnew[r * p->M + c] = cond * row[c]
                                     + (1-cond) * direct_weight * direct_neighbours
                                     + (1-cond) * diagonal_weight * diagonal_neighbours;

                // Update largest difference for convergence threshold
                diff = fabs(tnew[r * p->M + c] - told[r * p->M + c]);
                if (diff > maxdiff) {maxdiff = diff;}
            }
        }

        if ((i + 1) % p->period == 0) {
            calculate_results(p, result, tnew, i+1, maxdiff, t1);
        }

    }

    // Calculate results
    calculate_results(p, result, tnew, i, maxdiff, t1);

}
