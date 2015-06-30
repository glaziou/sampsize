/* double precision functions
 *   
 */
extern double bdtrc ( int k, int n, double p );
extern double bdtr ( int k, int n, double p );
extern double bdtri ( int k, int n, double y );
extern double btdtr ( double a, double b, double x );
extern double gamma ( double x );
extern double incbet ( double aa, double bb, double xx );
extern double incbi ( double aa, double bb, double yy0 );
extern int mtherr ( char *name, int code );
extern double ndtr ( double a );
extern double erfc ( double a );
extern double erf ( double x );
extern double ndtri ( double y0 );
extern double log1p ( double x );
extern double exp1m ( double x );
extern double polevl ( double x, void *P, int n );
extern double p1evl ( double x, void *P, int n );

/* These are presumed in libm. */
extern double exp ( double x );
extern double fabs ( double x );
extern double floor ( double x );
extern double log ( double x );
extern double powi ( double x, int n );
extern double pow ( double x, double y );
extern double sin ( double x );
extern double cos ( double x );
extern double sqrt ( double x );
