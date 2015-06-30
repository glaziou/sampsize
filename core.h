/* 
 * Primary functions
 * 
 * Sampsize -- sample size and power determination
 *
 * Copyright (C) 2003 Philippe Glaziou
 */

typedef struct s_Bounds
{ 
  double ll; 
  double ul;
}
Bounds; /* holds exact confidence limits */
    
	/* functions in errors.c */
extern void check_param (const Cmdline *);
extern int sperror (const char *);

	/* functions in size.c */
extern double sampsi (const Cmdline *);
extern Bounds cbounds (double, double, double);
extern Bounds spbounds (const Cmdline *);
extern double getza (const Cmdline *);
extern void display_surv (const Cmdline *, double, int);
extern double small_sampsi (const Cmdline *);
extern void display_small (const Cmdline *, double);
extern void binom_ci (const Cmdline *);
extern void comp (const Cmdline *);
extern void means (const Cmdline *);
extern void cluster (const Cmdline *, double);
extern void nequivp (const Cmdline *);
extern void nequivm (const Cmdline *);

	/* functions in ccontrol.c */
extern void case_control (const Cmdline *);
extern void display_cc (const Cmdline *, double, double);
extern void ccpower (const Cmdline *);
extern void ccmin (const Cmdline *);
extern void display_mcc (const Cmdline *, double, double, double);

	/* functions in power.c */
extern void ppower (const Cmdline *);
extern void mpower (const Cmdline *);
extern void mccpower (const Cmdline *);
