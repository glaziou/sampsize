/* 
   power.c -- power determination functions

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
#define ERFC_LIM 37.519379

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))


const char onesided[] = " (one-sided)";
const char twosided[] = " (two-sided)";


	/* Power determination -- percentages */
void
ppower (const Cmdline * const par)
{
  double p0, p1, pb, diff, za, zb;
  double w0, w1, n0, n1, n2, r1;
  double tmp, ratio, power;
  double pm, ds, ex, sd;

  p0 = par->comp[0] / PC;
  p1 = par->comp[1] / PC;

  if (p0 == p1)
    sperror ("-cp parameters should not be equal to each other");

  diff = fabs (p0 - p1);
  w0 = p0 * (1 - p0);
  w1 = p1 * (1 - p1);
  za = getza (par);
  n1 = par->sample;
  n2 = ceil (n1 * par->ratio);

  if (par->onesampleP)		/* one sample */
    {
      tmp = (diff * sqrt (n1) - za * sqrt (w0)) / sqrt (w1);
      power = (fabs (tmp) < ERFC_LIM) ? ndtr (tmp) : 1;
      printf ("\nEstimated power for one-sample comparison of percentage\n");
      printf ("to hypothesized value\n\n");
      printf ("Test Ho:    p = %10g%%", p0 * PC);
      printf (", where p is the percentage in the population\n\n");
      printf ("Assumptions:\n\n");
      printf ("        alpha = %10g%%", par->alpha * PC);
      if (par->onesidedP)
	printf ("%s\n", onesided);
      else
	{
	  printf ("%s\n", twosided);
	  tmp = (-diff * sqrt (n1) - za * sqrt (w0)) / sqrt (w1);
	  power = (fabs (tmp) < ERFC_LIM) ? power + ndtr (tmp) : 1;
	  if (power > 1)
	    power = 1;
	}
      printf ("alternative p = %10g%%\n", p1 * PC);
      printf ("  sample size = %10g\n\n", n1);
      printf ("Estimated power:\n\n");
      printf ("        power = %10g%%\n\n", power * PC);
    }
  else				/* two samples */
    {
      ratio = n2 / n1;
      if (par->correctionP)
	{
	  r1 = ratio + 1;
	  pb = (p0 + ratio * p1) / r1;
	  n0 = (n1 - r1 / (2 * ratio * diff));
	  n0 *= n0;
	  n0 /= n1;
	  zb =
	    (diff * sqrt (ratio * n0) -
	     za * sqrt (r1 * pb * (1 - pb))) / sqrt (ratio * w0 + w1);
	  power = (fabs (zb) < ERFC_LIM) ? ndtr (zb) : 1;

	}
      else			/* no continuity correction */
	{
	  pm = (n1 * p0 + n2 * p1) / (n1 + n2);
	  ds = za * sqrt ((1 / n1 + 1 / n2) * pm * (1 - pm));
	  ex = fabs (p0 - p1);
	  sd = sqrt (w0 / n1 + w1 / n2);
	  power = (fabs (((-ds - ex) / sd)) < ERFC_LIM) ?
	    1 - ndtr ((ds - ex) / sd) + ndtr ((-ds - ex) / sd) : 1;
	}
      printf
	("\nEstimated power for two-sample comparison of percentages\n\n");
      printf ("Test Ho:    p1 = p2");
      printf (", where p1 is the percentage in the population 1\n");
      printf ("and p2 is the percentage in the population 2\n\n");
      printf ("Assumptions:\n\n");
      printf ("         alpha = %10g", par->alpha * PC);
      if (par->onesidedP)
	printf ("%s\n", onesided);
      else
	printf ("%s\n", twosided);
      printf ("            p1 = %10g%%\n", p0 * PC);
      printf ("            p2 = %10g%%\n", p1 * PC);
      printf ("sample size n1 = %10g\n", n1);
      printf ("sample size n2 = %10g\n", n2);
      printf ("         n2/n1 = %10g\n", ratio);
      printf ("\nEstimated power:\n\n");
      printf ("        power = %10g%%", power * PC);
      if (par->correctionP)
	printf ("  (with continuity correction)\n\n");
      else
	printf ("\n\n");
    }
  if (par->obsclusP || par->numclusP)
    cluster (par, n1);
  exit (0);
}

	/* Power determination when comparing means */
void
mpower (const Cmdline * const par)
{
  double diff, za, zb;
  double w, n0, n1, n2, r1;
  double m1, m2, sd1, sd2;
  double tmp, ratio, power;

  n1 = par->sample;
  n2 = n1 * par->ratio;
  n2 = ceil (n2);
  m1 = par->means[0];
  m2 = par->means[1];
  sd1 = par->means[2];
  sd2 = par->means[3];

  if (sd1 < 0)
    sperror ("sd1 out of range");
  if (par->meansC == 4 && par->onesampleP)
    sperror
      ("only one sd(#) must be specified for one-sample comparison of means");
  if (par->meansC == 3 && !par->onesampleP)
    sperror
      ("two sd(#) must be specified for two-sample comparison of means");

  za = getza (par);
  diff = fabs (m1 - m2);

  if (par->onesampleP)		/* one sample */
    {
      tmp = diff * sqrt (n1) / sd1 - za;
      power = (fabs (tmp) < ERFC_LIM) ? ndtr (tmp) : 1;
      printf ("\nEstimated power for one-sample comparison of mean\n");
      printf ("to hypothesized value\n\n");
      printf ("Test Ho:    m = %10g%%", m1);
      printf (", where m is the mean in the population\n\n");
      printf ("Assumptions:\n\n");
      printf ("        alpha = %10g%%", par->alpha * PC);
      if (par->onesidedP)
	printf ("%s\n", onesided);
      else
	{
	  printf ("%s\n", twosided);
	  tmp = -diff * sqrt (n1) / sd1 - za;
	  power = (fabs (tmp) < ERFC_LIM) ? power + ndtr (tmp) : 1;
	  if (power > 1)
	    power = 1;
	}
      printf ("alternative m = %10g\n", m2);
      printf ("           sd = %10g\n", sd1);
      printf ("  sample size = %10g\n\n", n1);
      printf ("\nEstimated power:\n\n");
      printf ("        power = %10g%%\n\n", power * PC);
    }
  else				/* two samples */
    {
      ratio = n2 / n1;
      if (sd2 == 0)
	sd2 = sd1;
      w = sqrt ((sd1 * sd1) / n1 + (sd2 * sd2) / n2);
      tmp = diff / w - za;
      power = (fabs (tmp) < ERFC_LIM) ? ndtr (tmp) : 1;
      printf ("\nEstimated power for two-sample comparison of means\n\n");
      printf ("Test Ho:    m1 = m2");
      printf (", where m1 is the mean in the population 1\n");
      printf ("and m2 is the mean in the population 2\n\n");
      printf ("Assumptions:\n\n");
      printf ("         alpha = %10g%%", par->alpha * PC);
      if (par->onesidedP)
	printf ("%s\n", onesided);
      else
	{
	  printf ("%s\n", twosided);
	  tmp = -diff / w - za;
	  power = (fabs (tmp) < ERFC_LIM) ? power + ndtr (tmp) : 1;
	  if (power > 1)
	    power = 1;
	}
      printf ("            m1 = %10g\n", m1);
      printf ("            m2 = %10g\n", m2);
      printf ("           sd1 = %10g\n", sd1);
      printf ("           sd2 = %10g\n", sd2);
      printf ("sample size n1 = %10g\n", n1);
      printf ("sample size n2 = %10g\n", n2);
      printf ("         n2/n1 = %10g\n", ratio);
      printf ("\nEstimated power:\n\n");
      printf ("         power = %10g%%\n\n", power * PC);
    }
  if (par->obsclusP || par->numclusP)
    cluster (par, n1);
  exit (0);
}

	/* Power of a case-control study */
void
ccpower (const Cmdline * const par)
{
  double alpha, power;
  double za, zb, c, n;
  double p0, q0, p1, q1, p, q;
  double controls;

  c = par->ratio;
  p0 = par->exposed;
  n = par->sample;

  if (!par->exposed)
    {
      printf ("\nMinimum options and minimum values of parameters\n");
      printf ("for a case-control study power computation:\n\n");
      printf ("\tsampsize -or xx -exp yy -n zz\n");
      printf ("\t(yy is the percentage of exposed controls)\n\n");
      printf ("\txx > 0    and    0 < yy < 100\n\n");
      exit (1);
    }

  alpha = (par->onesidedP) ? 2 * par->alpha : par->alpha;
  za = getza (par);
  p1 = (p0 * par->odds_ratio) / (1 + p0 * (par->odds_ratio - 1));
  p = (p1 + c * p0) / (1 + c);
  q = 1 - p;
  q1 = 1 - p1;
  q0 = 1 - p0;

  zb = sqrt (n * (p1 - p0) * (p1 - p0)) - za * sqrt ((1 + 1 / c) * p * q);
  zb /= sqrt (p1 * q1 + p0 * q0 / c);
  power = (fabs (zb) < ERFC_LIM) ? ndtr (zb) : 1;

  controls = rint (n * c);

  printf ("\nAssumptions:\n\n");
  printf ("     Odds ratio             = %10g\n", par->odds_ratio);
  printf ("     Exposed controls       = %10g%%\n", p0 * PC);
  if (par->onesidedP)
    printf ("     One-sided alpha risk   = %10g%%\n", par->alpha * PC);
  else
    printf ("     Alpha risk             = %10g%%\n", par->alpha * PC);
  printf ("     Controls / Case ratio  = %10g\n\n", c);
  printf ("     Total exposed          = %10g%%\n\n", p * PC);
  printf ("     Number of cases        = %10g\n", n);
  printf ("     Number of controls     = %10g\n\n", controls);
  printf ("     Total                  = %10g\n\n", n + controls);
  printf ("Estimated power:\n\n");
  printf ("                      power = %10g%%\n\n", power * PC);

  exit (0);
}


	/* Minimum detectable odds-ratio */
void
ccmin (const Cmdline * const par)
{
  double or1, or2, aa, bb, cc;
  double n, p0, za, zb, alpha;

  p0 = par->exposed;
  n = par->sample;

  alpha = (par->onesidedP) ? 2 * par->alpha : par->alpha;
  za = getza (par);
  zb = ndtri (par->power);

  aa = za + zb, aa *= aa;
  bb = 1 + 2 * p0;
  cc = 2 * p0 * (n * (1 - p0) - aa * p0);
  or1 = 1 + sqrt (aa) * (bb * sqrt (aa) + sqrt (aa * bb * bb + 4 * cc)) / cc;
  or2 = 1 + sqrt (aa) * (bb * sqrt (aa) - sqrt (aa * bb * bb + 4 * cc)) / cc;
  if (or2 < 0)
    or2 = 0;

  printf ("\nAssumptions:\n\n");
  printf ("     Exposed controls       = %10g%%\n", p0 * PC);
  if (par->onesidedP)
    printf ("     One-sided alpha risk   = %10g%%\n", par->alpha * PC);
  else
    printf ("     Alpha risk             = %10g%%\n", par->alpha * PC);
  printf ("     Power                  = %10g\n", par->power * PC);
  printf ("     Number of cases        = %10g\n", n);
  printf ("     Number of controls     = %10g\n\n", n);
  printf ("     Total                  = %10g\n\n", 2 * n);
  printf ("Estimated smallest and highest detectable odds-ratio:\n\n");
  printf ("   Highest  odds-ratio < 1  = %10g\n", min (or1, or2));
  printf ("   Smallest odds-ratio > 1  = %10g\n\n", max (or1, or2));

  exit (0);
}

void
mccpower (const Cmdline * const par)
{
  double m, za, zb, pp, power;
  double alpha, p0, q0, p1, q1, pe, mm;

  p0 = par->exposed;
  mm = par->sample;
  alpha = (par->onesidedP) ? 2 * par->alpha : par->alpha;
  za = getza (par);
  p1 =
    (par->exposed * par->odds_ratio) / (1 +
					par->exposed * (par->odds_ratio - 1));
  q1 = 1 - p1;
  q0 = 1 - p0;
  pp = par->odds_ratio / (1 + par->odds_ratio);
  pe = (p0 * q1 + p1 * q0);
  m = mm * pe;
  zb = -za / 2 + sqrt (m * (pp - 0.5) * (pp - 0.5));
  zb /= sqrt (pp * (1 - pp));

  power = (fabs (zb) < ERFC_LIM) ? ndtr (zb) : 1;

  printf ("\nAssumptions:\n\n");
  printf ("     Odds ratio                 = %10g\n", par->odds_ratio);
  printf ("     Exposed controls           = %10g%%\n", p0 * PC);
  if (par->onesidedP)
    printf ("     One-sided alpha risk       = %10g%%\n", par->alpha * PC);
  else
    printf ("     Alpha risk                 = %10g%%\n", par->alpha * PC);
  printf ("     Number of pairs            = %10g\n\n", mm);
  printf ("     Probability of an\n");
  printf ("    exposure-discordant pair    = %10g%%\n\n", pe * PC);
  printf ("Estimated power:\n\n");
  printf ("                          power = %10g%%\n\n", power * PC);

  exit (0);
}
