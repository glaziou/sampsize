#ifndef __cmdline__
#define __cmdline__
/*****
  command line parser interface -- generated by clig 
  (http://wsd.iitb.fhg.de/~kir/clighome/)

  The command line parser `clig':
  (C) 1995---2001 Harald Kirsch (kirschh@lionbioscience.com)
*****/

typedef struct s_Cmdline {
  /***** -h: show usage information */
  char show_helpP;
  /***** -v: show program version */
  char show_versionP;
  /***** -onesided: one-sided test for comparisons */
  char onesidedP;
  /***** -onesample: one-sample test */
  char onesampleP;
  /***** -matched: 1:1 matched case-control study */
  char matchedP;
  /***** -cc: continuity correction */
  char correctionP;
  /***** -pop: target population size */
  char popP;
  float pop;
  int popC;
  /***** -e: precision of estimate (%) -- prevalence study option */
  char precisionP;
  float precision;
  int precisionC;
  /***** -pr: prevalence (%) -- prevalence study option */
  char prevalenceP;
  float prevalence;
  int prevalenceC;
  /***** -level: level of confidence interval (%) */
  char levelP;
  float level;
  int levelC;
  /***** -alpha: Risk alpha (%) */
  char alphaP;
  float alpha;
  int alphaC;
  /***** -power: Power of the test (%) */
  char powerP;
  float power;
  int powerC;
  /***** -c: number of controls per case, to be used with case control options */
  char ratioP;
  float ratio;
  int ratioC;
  /***** -or: Odds Ratio -- case-control study option */
  char odds_ratioP;
  float odds_ratio;
  int odds_ratioC;
  /***** -exp: Exposed controls (%) -- case-control study option */
  char exposedP;
  float exposed;
  int exposedC;
  /***** -cp: #p1 (%) #p2 (%) -- two sample comparison of percentages */
  char compP;
  float *comp;
  int compC;
  /***** -means: #m1 #m2 #sd1 #sd2 -- two sample comparison of means */
  char meansP;
  float *means;
  int meansC;
  /***** -rho: intraclass correlation coefficient -- cluster option */
  char rhoP;
  float rho;
  int rhoC;
  /***** -d: delta critical value -- equivalence */
  char deltaP;
  float delta;
  int deltaC;
  /***** -nob: minimum number of events observed in the sample */
  char observedP;
  int observed;
  int observedC;
  /***** -bi: #obs #succ (binomial events, #obs >= #succ) */
  char binomialP;
  int *binomial;
  int binomialC;
  /***** -n: determine power from sample size */
  char sampleP;
  int sample;
  int sampleC;
  /***** -obsclus: number of observations per cluster */
  char obsclusP;
  int obsclus;
  int obsclusC;
  /***** -numclus: minimum number of clusters */
  char numclusP;
  int numclus;
  int numclusC;
  /***** uninterpreted command line parameters */
  int argc;
  /*@null*/char **argv;
  /***** the whole command line concatenated */
  char *tool;
} Cmdline;


extern char *Program;
extern void usage(void);
extern /*@shared*/Cmdline *parseCmdline(int argc, char **argv);

#endif

