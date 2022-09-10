/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */


#include <math.h>
#include <float.h>
#include <stdio.h>
#include <iostream>

namespace R {

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2005-6 Morten Welinder <terra@gnome.org>
 *  Copyright (C) 2005-10 The R Foundation
 *  Copyright (C) 2006-2015 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *
 *    double pgamma (double x, double alph, double scale,
 *               int lower_tail, int log_p)
 *
 *    double log1pmx    (double x)
 *    double lgamma1p (double a)
 *
 *    double logspace_add (double logx, double logy)
 *    double logspace_sub (double logx, double logy)
 *    double logspace_sum (double* logx, int n)
 *
 *
 *  DESCRIPTION
 *
 *    This function computes the distribution function for the
 *    gamma distribution with shape parameter alph and scale parameter
 *    scale.    This is also known as the incomplete gamma function.
 *    See Abramowitz and Stegun (6.5.1) for example.
 *
 *  NOTES
 *
 *    Complete redesign by Morten Welinder, originally for Gnumeric.
 *    Improvements (e.g. "while NEEDED_SCALE") by Martin Maechler
 *
 *  REFERENCES
 *
 */

/*----------- DEBUGGING -------------
 * make CFLAGS='-DDEBUG_p -g'
 * (cd `R-devel RHOME`/src/nmath; gcc -I. -I../../src/include -I../../../R/src/include  -DHAVE_CONFIG_H -fopenmp -DDEBUG_p -g -c ../../../R/src/nmath/pgamma.c -o pgamma.o)
 */

#ifndef NDEBUG
    //#define DEBUG_p
    #define REprintf printf
#endif

/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
#define M_LN_SQRT_2PI   0.918938533204672741780329736406    /* log(sqrt(2*pi)) */
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))
#define give_log log_p
#define R_D__0    (log_p ? -INFINITY : 0.)
#define R_D__1    (log_p ? 0. : 1.)            /* 1 */
#define R_D_exp(x)    (log_p    ?  (x)     : exp(x))    /* exp(x) */
#define SQR(x) ((x)*(x))
#define R_D_fexp(f,x) (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR

#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#else
#  define LDOUBLE double
#endif

#ifdef HAVE_NEARYINT
# define R_forceint(x)   nearbyint()
#else
# define R_forceint(x)   round(x)
#endif

#define R_P_bounds_01(x, x_min, x_max)    \
    if(x <= x_min) return R_DT_0;        \

/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]     =~=  -x */
static const double M_cutoff = M_LN2 * DBL_MAX_EXP / DBL_EPSILON;/*=3.196577e18*/

/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 *
 * auxilary in log1pmx() and lgamma1p()
 */

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif
#define M_LN_SQRT_PId2    0.225791352644727432363097614947
#define     M_2PI   6.28318530717958647692528676655900576
#define M_1_SQRT_2PI    0.398942280401432677939946059934    /* 1/sqrt(2pi) */
#define M_SQRT_32    5.656854249492380195206754896838
#define R_DT_0    (lower_tail ? R_D__0 : R_D__1)        /* 0 */
#define R_DT_1    (lower_tail ? R_D__1 : R_D__0)        /* 1 */

#ifdef HAVE_SINPI
#elif defined HAVE___SINPI
double sinpi(double x);
#else
// sin(pi * x)  -- exact when x = k/2  for all integer k
double sinpi(double x);
#endif

double dnorm(double x, double mu, double sigma, int give_log);

double attribute_hidden chebyshev_eval(double x, const double *a, const int n) ;

double attribute_hidden lgammacor(double x);

double attribute_hidden stirlerr(double n);

double gammafn(double x);

double lgammafn_sign(double x, int *sgn);

double lgammafn(double x);

// static double
// logcf (double x, double i, double d,
//        double eps /* ~ relative tolerance */);

/* Accurate calculation of log(1+x)-x, particularly for small x.  */
double log1pmx (double x);

/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
double lgamma1p (double a);
/* lgamma1p */


double fmax2(double x, double y);

/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
double logspace_add (double logx, double logy);

/*
 * Compute the log of a difference from logs of terms, i.e.,
 *
 *     log (exp (logx) - exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
double logspace_sub (double logx, double logy);

/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (sum_i  exp (logx[i]) ) =
 *     log (e^M * sum_i  e^(logx[i] - M) ) =
 *     M + log( sum_i  e^(logx[i] - M)
 *
 * without causing overflows or throwing much accuracy.
 */
#ifdef HAVE_LONG_DOUBLE
# define EXP expl
# define LOG logl
#else
# define EXP exp
# define LOG log
#endif
double logspace_sum (const double* logx, int n);

double attribute_hidden bd0(double x, double np);

double dpois_raw(double x, double lambda, int give_log);

/* dpois_wrap (x__1, lambda) := dpois(x__1 - 1, lambda);  where
 * dpois(k, L) := exp(-L) L^k / gamma(k+1)  {the usual Poisson probabilities}
 *
 * and  dpois*(.., give_log = TRUE) :=  log( dpois*(..) )
*/
// static double
// dpois_wrap (double x_plus_1, double lambda, int give_log);

/*
 * Abramowitz and Stegun 6.5.29 [right]
 */
// static double
// pgamma_smallx (double x, double alph, int lower_tail, int log_p);

// static double
// pd_upper_series (double x, double y, int log_p);

/* Continued fraction for calculation of
 *    scaled upper-tail F_{gamma}
 *  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
 */
// static double
// pd_lower_cf (double y, double d);
#undef NEEDED_SCALE


// static double
// pd_lower_series (double lambda, double y);

/*
 * Compute the following ratio with higher accuracy that would be had
 * from doing it directly.
 *
 *         dnorm (x, 0, 1, FALSE)
 *       ----------------------------------
 *       pnorm (x, 0, 1, lower_tail, FALSE)
 *
 * Abramowitz & Stegun 26.2.12
 */
#define SIXTEN    16 /* Cutoff allowing exact "*" and "/" */

void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);

// static double
// dpnorm (double x, int lower_tail, double lp);

/*
 * Asymptotic expansion to calculate the probability that Poisson variate
 * has value <= x.
 * Various assertions about this are made (without proof) at
 * http://members.aol.com/iandjmsmith/PoissonApprox.htm
 */

double pnorm(double x, double mu, double sigma, int lower_tail, int log_p);


// static double
// ppois_asymp (double x, double lambda, int lower_tail, int log_p);

double pgamma_raw (double x, double alph, int lower_tail, int log_p);

double pgamma(double x, double alph, double scale, int lower_tail, int log_p);

/* From: terra@gnome.org (Morten Welinder)
 * To: R-bugs@biostat.ku.dk
 * Cc: maechler@stat.math.ethz.ch
 * Subject: Re: [Rd] pgamma discontinuity (PR#7307)
 * Date: Tue, 11 Jan 2005 13:57:26 -0500 (EST)

 * this version of pgamma appears to be quite good and certainly a vast
 * improvement over current R code.  (I last looked at 2.0.1)  Apart from
 * type naming, this is what I have been using for Gnumeric 1.4.1.

 * This could be included into R as-is, but you might want to benefit from
 * making logcf, log1pmx, lgamma1p, and possibly logspace_add/logspace_sub
 * available to other parts of R.

 * MM: I've not (yet?) taken  logcf(), but the other four
 */

/*
int main(int argc, const char * argv[]) {
    double x=10;
    double df=135;
    int lower_tail=0;
    int log_p=1;
    std::cout << pgamma(x, df/2., 2., lower_tail, log_p) << "\n";
}
*/

} // R
