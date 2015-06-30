/* 
   size.c -- sample size determination functions. 

   Copyright (C) 2003 Philippe Glaziou

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include <stdio.h>
#include <math.h>

#include "cmdline.h"
#include "protos.h"
#include "core.h"

#define PC 100
#define MAXSIZE 1.0e8


	/* Returns exact confidence interval limits of a proportion */
Bounds
cbounds (double k, double n, double level)
{
  Bounds ci;
  if (k < n)
    ci.ul = bdtri (k, n, (1 - level) / 2);
  else if (k >= n)
    ci.ul = 1;
  if (k > 0)
    ci.ll = bdtri (k - 1, n, 1 - (1 - level) / 2);
  else
    ci.ll = 0;
  return (ci);
}


	/* Returns specified confidence interval */
Bounds
spbounds (const Cmdline * const par)
{
  Bounds ci;
  ci.ll = (par->prevalence - par->precision) * PC;
  if (ci.ll < 0)
    ci.ll = 0;

  ci.ul = (par->prevalence + par->precision) * PC;
  if (ci.ul > 100)
    ci.ul = 100;
  return (ci);
}


	/* returns za */
double
getza (const Cmdline * const par)
{
  return (par->onesidedP) ?
    ndtri (1 - par->alpha) : ndtri (1 - (par->alpha / 2));
}


	/* Returns an exact binomial estimate */
double
small_sampsi (const Cmdline * const par)
{
  double cubinom = 0.0;
  double top, mid;
  double bot = 0;

  mid = par->observed;

  while (cubinom < par->level)
    {
      cubinom =
	incbet (par->observed, (mid - par->observed + 1), par->prevalence);
      mid *= 2;
    }

  do
    {
      if (cubinom >= par->level)
	top = mid;
      else
	bot = mid;
      mid = (bot + top) / 2;
      cubinom =
	incbet (par->observed, (mid - par->observed + 1), par->prevalence);
    }
  while (fabs (bot - mid) > 0.5);

  return ceil (top);
}


	/* Display small results */
void
display_small (const Cmdline * const par, double nb)
{
  printf ("\nAssumptions:\n\n");
  printf ("   prevalence                  = %8g%%\n", par->prevalence * PC);
  printf
    ("   minimum number of successes =       %d  in the sample\n",
     par->observed);
  printf ("   with probability            = %8g%%\n\n", par->level * PC);
  printf ("Estimated sample size:\n");
  printf ("                             n = %8g\n\n", nb);

  exit (0);
}


	/* Returns exact sample size if prevalence null
	 * and approximate if prevalence > 0
	 */
double
sampsi (const Cmdline * const par)
{
  double q, nb, alpha, pq, s2;

  q = 1 - par->prevalence;
  pq = par->prevalence * (1 - par->prevalence);
  s2 = par->precision * par->precision;

  alpha = ndtri (1 - (1 - par->level) / 2);
  alpha *= alpha;

  if (par->prevalence == 0)
    {
      /* This formula is a simplification of the
       * cumulative binomial when k = 0,
       * it will return an exact sample size
       */
      nb = log (1 - (1 - (1 - par->level) / 2)) / log (1 - par->precision);
    }
  else
    {
      if (par->pop > 0)
	/* The following formula includes a finite
	 * population correction to the standard formula
	 * based on the normal distribution */
	nb =
	  (par->pop * (alpha / s2) * pq) / (par->pop - 1 + (alpha / s2) * pq);
      else
	nb = alpha * pq / s2;
    }

  return ceil (nb);
}


	/* Display survey results */
void
display_surv (const Cmdline * par, double nb, int infpop)
{
  double level, k;
  double maxsize = MAXSIZE;
  Bounds cexact, cspec;

  level = par->level;

  printf ("\nAssumptions:\n\n");
  printf ("     Precision          = %3.2f  %%\n", par->precision * PC);
  printf ("     Prevalence         = %3.2f %%\n", par->prevalence * PC);
  (!infpop) ? printf ("     Population size    = %.0f\n\n", par->pop)
    : printf ("     Population size    = infinite\n\n");

  if (par->prevalence == 0)
    {
      level = 1 - (1 - level) / 2;

      printf
	("A null prevalence was entered, we assume here that the observed \n");
      printf
	("prevalence in the sample will be null, and compute the sample\n");
      printf
	("size for an upper limit of the confidence interval equal to the\n");
      printf
	("entered precision. The population size is not taken into account\n");
      printf
	("even if a finite size was entered. Sample sizes are calculated\n");
      printf ("using the binomial exact method.\n\n");
      printf ("%10g%% one-sided Confidence Interval upper limit = %3.2f\n\n",
	      level * PC, par->precision * PC);
      printf ("Estimated sample size:\n");
      printf ("                      n = %10g\n\n", nb);
      exit (0);
    }

  k = nb * (par->prevalence);
  if (k - floor (k) >= .5)
    k = floor (k) + 1;
  else
    k = floor (k);
  cspec = spbounds (par);	/* specified confidence bounds */

  printf
    ("%10g%% Confidence Interval specified limits [%4g%%   -- %4g%% ]\n",
     level * PC, cspec.ll, cspec.ul);
  printf ("\t(these limits equal prevalence plus or minus precision)\n\n");
  printf ("Estimated sample size:\n");
  printf ("                      n = %8g\n\n", nb);

  if (infpop && nb < maxsize)
    {
      cexact = cbounds (k, nb, level);	/* exact confidence bounds */
      printf
	("%10g%% Binomial Exact Confidence Interval with n = %8g\n",
	 level * PC, nb);
      printf ("\tand n * prevalence = %6g observed events:\n", k);
      printf ("\t[%6g%%   -- %6g%% ]\n\n", cexact.ll * PC, cexact.ul * PC);
    }
  if (k < 5 || (nb - k) < 5)
    {
      printf ("\nWarning: #succ or (sample size - #succ) < 5\n");
      printf ("You may try different combinations of #obs and #succ\n");
      printf ("and assess binomial exact confidence intervals.\n\n");
    }
  exit (0);
}


	/* Display binomial exact confidence bounds */
void
binom_ci (const Cmdline * const par)
{
  Bounds ci;
  int n, k;
  double level;

  level = par->level;
  n = par->binomial[0];
  k = par->binomial[1];

  if (k > n)
    sperror ("#obs < #succ");

  ci = cbounds (k, n, level);

  printf ("\nAssumptions:\n\n");
  printf ("   number of observations = %d\n", n);
  printf ("   number of successes    = %d\n\n", k);
  printf ("Binomial Exact ");
  if (ci.ll == 0 || ci.ul == 1)
    printf ("one-sided %4g%% Confidence Interval: \n   [%6g%% -- %6g%%]\n\n",
	    PC * (1 - (1 - level) / 2), ci.ll * PC, ci.ul * PC);
  else
    printf ("%4g%% Confidence Interval: \n   [%6g%% -- %6g%%]\n\n",
	    level * PC, ci.ll * PC, ci.ul * PC);
  exit (0);
}


	/* Comparison of percentages */
void
comp (const Cmdline * const par)
{
  double p, p0, p1, pb, diff, za, zb;
  double w0, w1, n0, n1, n2, r1, fraction;

  p0 = par->comp[0] / PC;
  p1 = par->comp[1] / PC;
  if (p0 == p1)
    sperror ("-cp parameters should not be equal to each other");

  diff = fabs (p0 - p1);
  w0 = p0 * (1 - p0);
  w1 = p1 * (1 - p1);
  za = getza (par);
  zb = ndtri (par->power);

  if (par->onesampleP)
    {
      n1 = (za * sqrt (w0) + zb * sqrt (w1)) / diff;
      n1 *= n1;
      n1 = ceil (n1);
      printf
	("\nEstimated sample size for one-sample comparison of percentage to\n");
      printf ("hypothesized value\n\n");
      printf ("Test Ho:     p = %10g%%, ", par->comp[0]);
      printf ("where p is the percentage in the population\n\n");
      printf ("Assumptions:\n");
      printf ("         alpha = %10g%%", par->alpha * PC);
      (par->onesidedP) ?
	printf (" (one-sided)\n") : printf (" (two-sided)\n");
      printf ("         power = %10g%%\n", par->power * PC);
      printf (" alternative p = %10g%%\n\n", par->comp[1]);
      printf ("Estimated sample size:\n");
      printf ("             n = %10g\n\n", n1);
    }
  else
    {
      if (par->correctionP)
	{
	  r1 = par->ratio + 1;
	  pb = (p0 + par->ratio * p1) / r1;
	  n0 =
	    (za * sqrt (r1 * pb * (1 - pb)) +
	     zb * sqrt (par->ratio * w0 + w1));
	  n0 *= n0;
	  n0 /= (par->ratio * diff * diff);
	  n1 = 1 + sqrt (1 + 2 * r1 / (n0 * par->ratio * diff));
	  n1 *= n1;
	  n1 *= n0 / 4;
	  n1 = ceil (n1);
	  n2 = ceil (par->ratio * n1);
	}
      else			/* no continuity correction */
	{
	  fraction = 1 / (par->ratio + 1);
	  p = fraction * p0 + (1 - fraction) * p1;
	  n1 =
	    za * sqrt ((par->ratio + 1) * p * (1 - p)) +
	    zb * sqrt (par->ratio * p0 * (1 - p0) + p1 * (1 - p1));
	  n1 *= n1 / par->ratio / ((p0 - p1) * (p0 - p1));
	  n2 = n1 * par->ratio;
	  n1 = ceil (n1);
	  n2 = ceil (n2);
	}

      printf
	("\nEstimated sample size for two-sample comparison of percentages\n");
      printf
	("\nTest H:    p1 = p2, where p1 is the percentage in population 1\n");
      printf ("and p2 is the percentage in population 2\n\n");
      printf ("Assumptions:\n");
      printf ("        alpha = %10g%%", par->alpha * PC);
      (par->
       onesidedP) ? printf (" (one-sided)\n") : printf (" (two-sided)\n");
      printf ("        power = %10g%%\n", par->power * PC);
      printf ("           p1 = %10g%%\n", par->comp[0]);
      printf ("           p2 = %10g%%\n\n", par->comp[1]);
      printf ("Estimated sample size");
      if (par->correctionP)
	printf (" (with continuity correction)");
      printf (":\n");
      printf ("           n1 = %10g\n", n1);
      printf ("           n2 = %10g\n\n", n2);
    }

  if (par->obsclusP || par->numclusP)
    cluster (par, n1);
  exit (0);
}


	/* Comparison of means */
void
means (const Cmdline * const par)
{
  double m1, m2, sd1, sd2, diff, w;
  double za, zb, n1, n2;

  m1 = par->means[0];
  m2 = par->means[1];
  sd1 = par->means[2];
  sd2 = par->means[3];

  if (sd1 < 0)
    sperror ("sd1 out of range");
  if (sd2 < 0)
    sperror ("sd2 out of range");
  if (par->meansC == 4 && par->onesampleP)
    sperror
      ("only one sd(#) must be specified for one-sample comparison of means");
  if (par->meansC == 3 && !par->onesampleP)
    sperror
      ("two sd(#) must be specified for two-sample comparison of means");

  za = getza (par);
  zb = ndtri (par->power);
  diff = fabs (m1 - m2);
  w = (za + zb) / diff;
  w *= w;

  if (sd2 == 0)
    sd2 = sd1;
  if (par->onesampleP)
    {
      if (sd1 == 0)
	sperror ("sd1 should not be null in a one-sample comparison");
      n1 = w * sd1 * sd1;
      n1 = ceil (n1);
      printf
	("\nEstimated sample size for one-sample comparison of mean to\n");
      printf ("hypothesized value\n\n");
      printf
	("Test Ho:    m = %10g, where m is the mean in the population\n\n",
	 m1);
      printf ("Assumptions:\n");
      printf ("        alpha = %10g", par->alpha * PC);
      (par->onesidedP) ?
	printf (" (one-sided)\n") : printf (" (two-sided)\n");
      printf ("        power = %10g\n", par->power * PC);
      printf ("alternative m = %10g\n", m2);
      printf ("           sd = %10g\n", sd1);
      printf ("\nEstimated sample size:\n");
      printf ("            n = %10g\n\n", n1);
    }
  else
    {
      if (sd1 == 0)
	sd1 = sd2;
      n1 = w * (sd1 * sd1 + sd2 * sd2 / par->ratio);
      n1 = ceil (n1);
      n2 = n1 * par->ratio;
      n2 = ceil (n2);

      printf
	("\nEstimated sample size for two-sample comparison of means\n\n");
      printf ("Test Ho:   m1 = m2, where m1 is the mean in population 1\n");
      printf ("and m2 is the mean in population 2\n\n");
      printf ("Assumptions: \n\n");
      printf ("        alpha = %10g", par->alpha * PC);
      (par->onesidedP) ?
	printf (" (one-sided)\n") : printf (" (two-sided)\n");
      printf ("        power = %10g\n", par->power * PC);
      printf ("           m1 = %10g\n", m1);
      printf ("           m2 = %10g\n", m2);
      printf ("          sd1 = %10g\n", sd1);
      printf ("          sd2 = %10g\n", sd2);
      printf ("        n2/n1 = %10g\n\n", par->ratio);
      printf ("Estimated sample size:\n\n");
      printf ("           n1 = %10g\n", n1);
      printf ("           n2 = %10g\n\n", n2);
    }
  if (par->obsclusP || par->numclusP)
    cluster (par, n1);

  exit (0);
}

void
cluster (const Cmdline * const par, double n)
{
  double obsclus, numclus;
  double n1, n2, ratio, obs, rh;
  double newn1, newn2, totaln, numcl, newobs;

  n1 = n;
  ratio = par->ratio;
  n2 = n * ratio;

  obs = par->obsclus;
  numcl = par->numclus;
  rh = par->rho;

  if (par->obsclusP)
    {
      newn1 = ceil (n1 * (1 + rh * (obs - 1)));
      newn2 = ceil (newn1 * ratio);
      totaln = newn1 + newn2;
      numcl = ceil (totaln / obs);
    }
  if (par->numclusP)
    {
      if (numcl <= (n1 * rh) + (n2 * rh))
	{
	  numcl = ceil ((n1 * rh) + (n2 * rh));
	  printf ("\nSample size adjusted for cluster design:\n\n");
	  printf
	    ("Error: for rho = %10g, the minimum number of clusters possible is:\n",
	     rh);
	  printf ("       numclus = %10g\n\n", numcl);
	  exit (1);
	}
      obs = (n1 - n1 * rh + n2 - n2 * rh) / (numcl - (n1 * rh) - (n2 * rh));
      obs = rint (obs);
      newn1 = ceil (n1 * (1 + rh * (obs - 1)));
      newn2 = ceil (newn1 * ratio);
      totaln = newn1 + newn2;
      newobs = ceil (totaln / numcl);
      if (newobs != obs)
	{
	  obs = newobs;
	  newn1 = ceil (n1 * (1 + rh * (obs - 1)));
	  newn2 = ceil (ratio * newn1);
	  totaln = newn1 + newn2;
	  numcl = ceil (totaln / obs);
	}
    }
  printf ("\nSample size adjusted for cluster design:\n\n");
  printf ("  Intraclass correlation     = %10g\n", rh);
  printf ("  Average obs. per cluster   = %10g\n", obs);
  printf ("  Minimum number of clusters = %10g\n\n", numcl);
  printf ("Estimated sample size per group:\n\n");
  printf ("           n1 (corrected)    = %10g\n", newn1);
  (!par->onesampleP) ? printf ("           n2 (corrected)    = %10g\n\n",
			       newn2) : printf ("\n");


  exit (0);
}

void
nequivp (const Cmdline * const par)
{
  double n, za, zb, ps, pt, delta, den;

  ps = par->comp[0] / PC;
  pt = par->comp[1] / PC;
  delta = par->delta / PC;
  za = ndtri (1 - par->alpha);
  zb = ndtri (par->power);
  n = (za + zb) * (za + zb) * (ps * (1 - ps) + pt * (1 - pt));
  den = ps - pt - delta;
  den *= den;
  n /= den;
  n = ceil (n);

  printf
    ("\nEstimated sample size for two-sample non-inferiority comparison of percentages\n");
  printf
    ("\nTest H0:    p2 < p1 - delta, where p1 is the percentage in population 1,\n");
  printf
    ("p2 is the percentage in population 2, and delta is the critical value\n\n");
  printf ("Assumptions:\n");
  printf ("        alpha = %10g%%", par->alpha * PC);
  (par->onesidedP) ? printf (" (one-sided)\n") : printf (" (two-sided)\n");
  printf ("        power = %10g%%\n", par->power * PC);
  printf ("           p1 = %10g%%\n", par->comp[0]);
  printf ("           p2 = %10g%%\n", par->comp[1]);
  printf ("        delta = %10g%%\n\n", par->delta);
  printf ("Estimated sample size per group:\n");
  printf ("           n = %10g\n\n", n);

  exit (0);
}



void
nequivm (const Cmdline * const par)
{
  double n, za, zb, delta, den;
  double ms, me, vs, ve;

  ms = par->means[0];
  me = par->means[1];
  vs = par->means[2];

  if (par->means[3] > 0)
    ve = par->means[3];
  else
    ve = vs;
  vs *= vs;
  ve *= ve;

  delta = par->delta;
  za = ndtri (1 - par->alpha);
  zb = ndtri (par->power);
  den = ms - me - delta;
  den *= den;
  n = (za + zb) * (za + zb) * (vs + ve) / den;
  n = ceil (n);

  printf
    ("\nEstimated sample size for two-sample non-inferiority comparison of means\n");
  printf
    ("\nTest H0:    diff < delta, where diff is the difference in means\n");
  printf ("Assumptions:\n");
  printf ("         alpha = %10g%%", par->alpha * PC);
  (par->onesidedP) ? printf (" (one-sided)\n") : printf (" (two-sided)\n");
  printf ("         power = %10g%%\n", par->power * PC);
  printf (" mean standard = %10g\n", ms);
  printf ("mean treatment = %10g\n", me);
  printf ("     std dev s = %10g\n", sqrt (vs));
  printf ("     std dev t = %10g\n", sqrt (ve));
  printf ("         delta = %10g\n\n", delta);
  printf ("Estimated sample size per group:\n");
  printf ("           n = %10g\n\n", n);

  exit (0);
}
