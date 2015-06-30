/* 
   sampsize -- sample size and power determination. 
   v. 0.6.0

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "system.h"
#include "cmdline.h"
#include "protos.h"
#include "mconf.h"
#include "core.h"

#define PC 100


int
main (int argc, char **argv)
{
  char *xmalloc ();
  char *xrealloc ();
  char *xstrdup ();
  int infpop;
  double nb;

  /* Check parameters */
  Cmdline *cmd = parseCmdline (argc, argv);

  if ((cmd->show_helpP) | (argc == 1))
    usage ();

  if (cmd->show_versionP)
    {
      printf ("%s %s\n", argv[0], VERSION);
      exit (0);
    }

  check_param (cmd);

  infpop = (cmd->pop == 0) ? 1 : 0;

  cmd->precision /= PC;
  cmd->prevalence /= PC;
  cmd->level /= PC;
  cmd->alpha /= PC;
  cmd->power /= PC;
  cmd->exposed /= PC;

  if (cmd->observedP)
    {
      nb = small_sampsi (cmd);
      display_small (cmd, nb);
    }

  else if (cmd->odds_ratioP && !cmd->sampleP)
    {
      cmd->ratio = floor (cmd->ratio);
      if (cmd->ratio < 1)
	sperror ("option -c should be >= 1");

      case_control (cmd);
    }

  /* Absolute precision then sample size equals population size */
  else if (cmd->precision == 0 && cmd->pop > 0)
    {
      nb = cmd->pop;
      display_surv (cmd, nb, infpop);
    }

  else if (cmd->precisionP)
    {
      nb = sampsi (cmd);
      display_surv (cmd, nb, infpop);
    }
  else if (cmd->binomialP)
    binom_ci (cmd);

  else if (cmd->compP && !cmd->sampleP && !cmd->deltaP)
    comp (cmd);

  else if (cmd->meansP && !cmd->sampleP && !cmd->deltaP)
    means (cmd);

  else if (cmd->sampleP && cmd->exposedP && !cmd->odds_ratioP && cmd->powerP
	   && !cmd->matchedP)
    ccmin (cmd);

  else if (cmd->sampleP && cmd->compP && !cmd->deltaP)
    ppower (cmd);

  else if (cmd->sampleP && cmd->meansP && !cmd->deltaP)
    mpower (cmd);

  else if (cmd->sampleP && cmd->odds_ratioP && cmd->exposedP
	   && !cmd->matchedP)
    ccpower (cmd);

  else if (cmd->sampleP && cmd->matchedP && cmd->odds_ratioP && cmd->exposedP)
    mccpower (cmd);

  else if (cmd->deltaP && cmd->compP && !cmd->sampleP)
    nequivp (cmd);

  else if (cmd->deltaP && cmd->meansP && !cmd->sampleP)
    nequivm (cmd);

  else
    sperror ("wrong combination of options, or missing options");

  exit (0);
}
