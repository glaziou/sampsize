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

#include "cmdline.h"
#include "core.h"


	/* Check for wrong combination of options and wrong parameters */
void
check_param (const Cmdline * const par)
{
  int flag = 0;
  if (par->compP && par->meansP)
    ++flag;
  if (par->compP && par->binomialP)
    ++flag;
  if (par->binomialP && par->meansP)
    ++flag;
  if (par->precisionP && (par->meansP || par->compP || par->binomialP))
    ++flag;
  if ((par->odds_ratioP || par->exposedP)
      && (par->meansP || par->compP || par->binomialP || par->precisionP
	  || par->onesampleP))
    ++flag;
  if (par->sampleP && par->binomialP)
    ++flag;

  if (flag)
    sperror ("wrong combination of options");

  if (par->precisionP && par->onesidedP)
    sperror ("option -onesided does not make sense with option -e");

  if (par->precisionP && par->onesampleP)
    sperror ("option -onesample does not make sense with option -e");

  if (par->observedP && par->onesidedP)
    sperror ("option -onesided does not make sense with option -nob");

  if (par->observedP && par->onesampleP)
    sperror ("option -onesample does not make sense with option -nob");

  if (par->observedP && par->prevalence == 0)
    sperror
      ("an event cannot occur if its probability of occurence equals zero");

  if (par->alphaP && par->alpha == 0 || par->alpha == 100)
    sperror ("-alpha parameter out of range");

  if (par->powerP && par->power == 100)
    sperror ("-power parameter out of range");

  if (par->levelP && par->level == 100)
    sperror ("-level parameter out of range");

  if (par->ratioP && par->ratio == 0)
    sperror ("-c parameter out of range");

  if ((par->ratio != 1) && par->matchedP)
    sperror ("-c may not be used with -matched");

  if (par->numclusP && par->obsclusP)
    sperror ("must choose either -numclus or -obsclus but not both");

  if ((par->numclusP || par->obsclusP) && !par->compP && !par->means)
    sperror ("option -cp or -means missing");

  if ((par->numclusP || par->obsclusP) && par->precisionP)
    sperror ("cluster options not allowed with option -e");

  if ((par->numclusP || par->obsclusP) && par->observedP)
    sperror ("cluster options not allowed with option -nob");

  if ((par->numclusP || par->obsclusP) && par->deltaP)
    sperror ("cluster options not allowed with option -d");

  if (par->deltaP && par->matchedP)
    sperror ("delta option not allowed with option -matched");

  if (par->deltaP && par->onesidedP)
    sperror ("delta option not allowed with option -onesided");

  if (par->deltaP && par->correctionP)
    sperror ("delta option not allowed with option -cc");

  if (par->odds_ratioP && par->correctionP)
    sperror ("or option not allowed with option -cc");

  if (par->sampleP && par->exposedP && par->matchedP && !par->odds_ratioP)
    sperror
      ("minimum/maximum detectable odds-ratio not available with\nmatched case-control designs");

  if (par->sampleP && (par->ratio != 1) && par->exposedP && !par->odds_ratioP)
    sperror
      ("minimum/maximum detectable odds-ratio not available with option -c");

  if (par->delta > 100 && par->compP)
    sperror ("delta between 0 and 100 if comparing two percentages");

  if (par->delta == 0 && (par->compP || par->means) && par->deltaP)
    sperror ("delta should be > 0");

  return;
}


	/* Display error messages and exit 1 */
int
sperror (const char error_text[])
{
  fprintf (stderr, "\nError: ");
  fprintf (stderr, "%s\n\n", error_text);
  exit (1);
}
