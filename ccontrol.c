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


void
case_control (const Cmdline * const par)
{
  double p1, q, p, nb, alpha;
  double za, zb;
  double p0, q0, q1, c;
  double pe, pp, m, mm;

  c = par->ratio;
  p0 = par->exposed;

  if (!par->exposed)
    {
      printf ("\nMinimum options and minimum values of parameters\n");
      printf ("for a case-control study sample size computation:\n\n");
      printf ("\tsampsize -or xx -exp yy\n");
      printf ("\t(yy is the percentage of exposed controls)\n\n");
      printf ("\txx > 0    and    0 < yy < 100\n\n");
      exit (1);
    }

  za = getza (par);
  zb = ndtri (1 - (1 - par->power));
  p1 =
    (par->exposed * par->odds_ratio) / (1 +
					par->exposed * (par->odds_ratio - 1));
  p = (p1 + c * par->exposed) / (1 + c);
  q = 1 - p;
  q1 = 1 - p1;
  q0 = 1 - p0;

  nb = za * sqrt ((1 + 1 / c) * p * q) + zb * sqrt (p1 * q1 + p0 * q0 / c);
  nb *= nb / ((p1 - p0) * (p1 - p0));
  nb = ceil (nb);

  if (par->matchedP)		/* matched */
    {
      pp = par->odds_ratio / (1 + par->odds_ratio);
      pe = (p0 * q1 + p1 * q0);
      m = za / 2 + zb * sqrt (pp * (1 - pp));
      m *= m / ((pp - 0.5) * (pp - 0.5));
      mm = m / pe;
      m = ceil (m);
      mm = ceil (mm);
      display_mcc (par, m, mm, pe);
    }

  display_cc (par, p, nb);
  exit (0);
}

	/* Display case-control results */
void
display_cc (const Cmdline * const par, double tte, double nb)
{
  printf ("\nAssumptions:\n\n");
  printf ("     Odds ratio             = %10g\n", par->odds_ratio);
  printf ("     Exposed controls       = %10g%%\n", par->exposed * PC);
  if (par->onesidedP)
    printf ("     One-sided alpha risk   = %10g%%\n", par->alpha * PC);
  else
    printf ("     Alpha risk             = %10g%%\n", par->alpha * PC);
  printf ("     Power                  = %10g%%\n", par->power * PC);
  printf ("     Controls / Case ratio  = %10g\n\n", par->ratio);
  printf ("     Total exposed          = %10g%%\n\n", tte * PC);
  printf ("Estimated sample size:\n\n");
  printf ("     Number of cases        = %10g\n", nb);
  printf ("     Number of controls     = %10g\n\n", nb * par->ratio);
  printf ("     Total                  = %10g\n\n", nb * (1 + par->ratio));
  exit (0);
}

	/* Display matched case-control results */
void
display_mcc (const Cmdline * const par, double m, double mm, double pe)
{
  printf ("\nAssumptions:\n\n");
  printf ("     Odds ratio             = %10g\n", par->odds_ratio);
  printf ("     Exposed controls       = %10g%%\n", par->exposed * PC);
  if (par->onesidedP)
    printf ("     One-sided alpha risk   = %10g%%\n", par->alpha * PC);
  else
    printf ("     Alpha risk             = %10g%%\n", par->alpha * PC);
  printf ("     Power                  = %10g%%\n", par->power * PC);
  printf ("     Probability of an\n");
  printf ("  exposure-discordant pair  = %10g%%\n\n", pe * PC);
  printf ("Estimated sample size (number of pairs):\n\n");
  printf ("     Number of exposure discordant-pairs  = %10g\n", m);
  printf ("     Number of pairs                      = %10g\n", mm);
  printf ("     Total sample size                    = %10g\n\n", 2 * mm);
  exit (0);
}
