/* -*- mode: c; kept-new-versions: 30; kept-old-versions: 20 -*- */

/* file lmrob.c
 * was	roblm/src/roblm.c - version 0.6	 by Matias Salibian-Barreras

 * Includes the stable correct asymptotic variance estimators
 * of Croux, Dhaene, Hoorelbeke
 * Includes the fast-s algorithm
*/

/* Robust MM regression estimates *
 * ------------------------------ */

/* comment code *
 *
 * adapt other sampler <<<<<<<<<<  R's random number generator !!!!

 * replace abort for too many singular resamples by
 * returning the number of singular ones
 */

/* MM:
   -  Done:  fast_s[_large_n]() both had FIXED seed (= 37),
	     and effectively discarded the seed_rand argument below
   -  Done:  drop 'register' : today's compilers do optimize well!
   -  Done:  using Calloc() / Free() instead of malloc()/free()
*/

/* kollerma:
   Added alternative psi functions callable via psifun, chifun and
   wgtfun. ipsi is used to distinguish between the different types.
   The ipsi argument works for the S-estimator as well as for the
   MM-estimator.
*/

#include <R.h>
#include <Rmath.h>

#include <R_ext/BLAS.h>


#include "robustbase.h"

/* these  will also move to "lmrob.h" ---
 *  but first make many of these 'static' <<< FIXME!
 */
void fast_s_large_n(double *X, double *y,
		    int *nn, int *pp, int *nRes,
		    int *ggroups, int *nn_group,
		    int *K, int *max_k, double *rel_tol, int *converged,
		    int *best_r, double *bb, double *rrhoc, int *iipsi,
		    double *bbeta, double *sscale, int *trace_lev);

void fast_s(double *X, double *y,
	    int *nn, int *pp, int *nRes,
	    int *K, int *max_k, double *rel_tol, int *converged,
	    int *best_r, double *bb, double *rrhoc, int *iipsi,
	    double *bbeta, double *sscale, int *trace_lev);

int rwls(double **a, int n, int p,
	 double *estimate, double *i_estimate,
	 double *resid, double *loss,
	 double scale, double epsilon,
	 int *max_it, const double rho_c, const int ipsi, int trace_lev);

#ifdef _NO_LONGUER_USED_
void sample_n_outof_N(int n, int N, int *x);
#endif
static void sample_noreplace(int *x, int n, int k, int *ind_space);


int lu(double **a,int *P, double *x);

double norm (double *x, int n);
double norm1(double *x, int n);
double norm_diff (double *x, double *y, int n);
double norm1_diff(double *x, double *y, int n);
double normcnst(double c, int ipsi);

double rho(double x, double c, int ipsi);
double psi(double x, double c, int ipsi);
double psip(double x, double c, int ipsi);
double wgt(double x, double c, int ipsi);

double rho_biwgt(double x, double c);
double psi_biwgt(double x, double c);
double psip_biwgt(double x, double c);
double wgt_biwgt(double x, double c);

double rho_gwgt(double x, double c);
double psi_gwgt(double x, double c);
double psip_gwgt(double x, double c);
double wgt_gwgt(double x, double c);

double rho_opt(double x, double c);
double psi_opt(double x, double c);
double psip_opt(double x, double c);
double wgt_opt(double x, double c);

double rho_hmpl(double x, double c);
double psi_hmpl(double x, double c);
double psip_hmpl(double x, double c);
double wgt_hmpl(double x, double c);

double rho_gws(double x, double c);
double psi_gws(double x, double c);
double psip_gws(double x, double c);
double wgt_gws(double x, double c);

double sum_rho	(double *x, int n, double c, int ipsi);
double sum_rho_sc(double *r, double scale, int n, int p, double c, int ipsi);

void get_weights_rhop(double *r, double s, int n,
		      double rhoc, int ipsi, double *w);

int refine_fast_s(double **x, double *y, double *weights,
		  int n, int p, double *res,
		  double *tmp, double *tmp2,
		  double **tmp_mat, double **tmp_mat2,
		  double *beta_cand,
		  int kk, Rboolean *conv, int max_k, double rel_tol,
		  int trace_lev,
		  double b, double rhoc, int ipsi, double initial_scale,
		  double *beta_ref, double *scale);

#ifdef _UNUSED_
void fast_s_irwls(double **x, double *y,
		  double *weights, int n, int p, double *beta_ref,
		  double **tmp_mat, double *tmp, double *tmp2);
#endif

int fast_s_with_memory(double **x, double *y,
		       int *nn, int *pp, int *nRes,
		       int *K, int *max_k, double *rel_tol, int *trace_lev,
		       int *best_r, int *ind_space, double *bb, double *rrhoc, int *iipsi,
		       double **best_betas, double *best_scales);

/* for "tracing" only : */
void disp_mat(double **a, int n, int m);
void disp_vec(double *a, int n);

double kthplace(double *, int, int);

int find_max(double *a, int n);

double find_scale(double *r, double b, double rhoc, int ipsi,
		  double initial_scale, int n, int p);

void r_sum_w_x(double **x, double *w, int n, int p,
	       double *tmp,
	       double *sum);
void r_sum_w_x_xprime(double **x, double *w, int n, int p,
		      double **tmp, double **ans);

double median_abs(double *, int, double *);
double MAD(double *a, int n, double center, double *tmp, double *tmp2);

double vecprime_vec(double *a, double *b, int n);

void mat_prime_mat_w(double **a, double *w, double **c, int n, int m);
void mat_prime_vec(double **a, double *b, double *c, int n, int m);
void outer_vec_vec(double **a, double *v1, double *v2, int n);

void zero_mat(double **a, int n, int m);
void zero_vec(double *a, int n);
void scalar_mat(double **a, double b, double **c, int n, int m);
void scalar_vec(double *a, double b, double *c, int n);
void sum_mat(double **a, double **b, double **c, int n, int m);
void sum_vec(double *a, double *b, double *c, int n);

#define SETUP_FAST_S(_n_, _p_)					\
    ind_space =  (int *) R_alloc(n, sizeof(int));		\
    res	    = (double *) R_alloc(_n_, sizeof(double));		\
    weights = (double *) R_alloc(_n_, sizeof(double));		\
    tmp	    = (double *) R_alloc(_n_, sizeof(double));		\
    tmp2    = (double *) R_alloc(_n_, sizeof(double));		\
    /* 'Matrices' (pointers to pointers): use Calloc(), 	\
     *  so we can Free() in precise order: */			\
    tmp_mat  = (double **) Calloc(_p_, double *);		\
    tmp_mat2 = (double **) Calloc(_p_, double *);		\
    for(i=0; i < _p_; i++) {					\
	tmp_mat [i] = (double *) Calloc( _p_,	 double);	\
	tmp_mat2[i] = (double *) Calloc((_p_)+1, double);	\
    }

/* This assumes that 'p' is correctly defined, and 'j' can be used in caller: */
#define COPY_beta(BETA_FROM, BETA_TO) \
  for(j=0; j < p; j++) BETA_TO[j] = BETA_FROM[j];
/* In theory BLAS should be fast, but this seems slightly slower,
 * particularly for non-optimized BLAS :*/
/* static int one = 1; */
/* #define COPY_beta(BETA_FROM, BETA_TO) \ */
/*     F77_CALL(dcopy)(&p, BETA_FROM, &one, BETA_TO, &one);  */

#define ZERO 1e-10
#define INFI 1e+20
/*UNUSED #define MAX_ITER_REFINE_S 50 */
#define MAX_NO_TRY_SAMPLES 1000 /* was 500 */
#define MAX_ITER_FIND_SCALE 200
#define TOL_INVERSE ZERO

/* This function computes an S-regression estimator */
void R_lmrob_S(double *X, double *y, int *n, int *P,
	       int *nRes, double *scale, double *beta_s,
	       double *rrhoc, int *iipsi, double *bb,
	       int *best_r, int *Groups, int *N_group,
	       int *K_s, int *max_k, double *rel_tol, int *converged,
	       int *trace_lev)
{
  /* best_r = 't' of Salibian-Barrera_Yohai(2006),
   *	      = no. of best candidates to be iterated further
   *		("refined") */
  
  /* Rprintf("R_lmrob_s %d\n", *iipsi); */
  
  if ( *nRes > 0) { 
    if( *n > 2000 )
      fast_s_large_n(X, y, n, P, nRes,
		     Groups, N_group,
		     K_s, max_k, rel_tol, converged,
		     best_r, bb, rrhoc, iipsi, beta_s, scale, trace_lev);
    else
      fast_s(X, y, n, P, nRes,
	     K_s, max_k, rel_tol, converged,
	     best_r, bb, rrhoc, iipsi, beta_s, scale, trace_lev);
  } else {
    *scale = find_scale(y, *bb, *rrhoc, *iipsi, *scale, *n, *P);
  }
}

/* This function performs RWLS iterations starting from
 * an S-regression estimator (and associated residual scale).
 * So, in itself, this is ``just'' an M-estimator :
 *
 * NOTE: rel_tol now controls the *relative* changes in beta,
 *      instead of being hard-wired to EPS = 1e-7 and bounding the
 *	absolute || beta_1 - beta_2 ||
 */
void R_lmrob_MM(double *X, double *y, int *n, int *P,
		double *beta_initial, double *scale,
		double *beta_m, double *resid,
		int *max_it, double *rho_c, int *ipsi, double *loss,
		double *rel_tol, int *converged, int *trace_lev)
{
    double **x;
    int N = *n, p = *P, i, j;
    x = (double **) Calloc(N, double *);

    /* rearranges X into a matrix of n x p,
     * i.e. [, 0:(p-1)] and cbind()'s  y as column [,p] : */
    for(i=0; i < N; i++) {
	double *x_i = x[i] = (double *) Calloc( (p+1), double);
	for(j=0; j < p; j++)
	    x_i[j]=X[j*N+i];
	x_i[p]=y[i];
    }

/* starting from the S-estimate (beta_initial), use
 * irwls to compute the MM-estimate (beta_m)  */

    *converged = rwls(x,N,p,beta_m, beta_initial, resid, loss,
		      *scale, *rel_tol,
		      max_it, *rho_c, *ipsi, *trace_lev);
    if (!converged)
	COPY_beta(beta_initial, beta_m);

    for(i=0; i < N; i++)
	Free(x[i]);
    Free(x);
}


/* This function solves a linear system of equations.
 * It solves for "x"
 *	a[, 0:(p-1)] x = a[,p]
 * using the LU decomposition of the p x p matrix  a[0:(p-1), 0:(p-1)]
 */
int lu(double **a, int *P, double *x)
{
    int *pp; /* pp[] : vector storing the permutations */
    int i,j,k, p = *P;
    double s;

    if ((pp = (int *) Calloc(p, int))==NULL)
	return(1);

    for(j=0; j < p; j++) { /* cols */
	pp[j]=j;
	for(i=j; i < p; i++)   /* rows */
	    if ( fabs( a[i][j] ) > fabs( a[pp[j]][j] ) )
		pp[j]=i;
	if ( pp[j] != j ) { /* permute rows */
	    double *kk;
	    kk= a[j]; a[j] = a[pp[j]]; a[pp[j]] = kk;
	}
	/* return if singular (det=0)
	 * if pivot (j,j) is "zero"	 */
	if ( fabs(a[j][j]) < TOL_INVERSE ) {
	    Free(pp);
	    return(1);
	}
	for(k=(j+1); k < p; k++)
	    a[k][j] = a[k][j] / a[j][j];
	for(k=(j+1); k < p; k++)
	    for(i=(j+1); i < p; i++)
		a[k][i] = a[k][i] - a[k][j] * a[j][i];

    } /* end of j for loop*/
    for(i=0; i < p; i++) {
	s=0.;
	for(j=0; j < i; j++)
	    s += a[i][j] * x[j];
	x[i] = a[i][p] - s; /* y[i]=a[i][p] */
    }
    for(i=(p-1); i>=0; i--) {
	s=0;
	for(j=(i+1); j < p; j++)
	    s += a[i][j] * x[j];
	x[i] = (x[i] - s) / a[i][i];
    }
    Free(pp);
    return(0);
}


void R_psifun(double *xx, double *cc, int *iipsi, int *dderiv, int *llength) {
  /* 
   * Calculate psi for vectorized x with psip(0) = 1
   * deriv -1: rho (not normalized)
   * deriv 0: psi
   * deriv 1: psip (psip(0) = 1)
   */
  
  int i;
  double nc;
  nc = normcnst(*cc, *iipsi);

  for(i = 0; i < *llength; i++) {
    switch(*dderiv) {
    case -1: xx[i] = rho(xx[i], *cc, *iipsi) / nc; break;
    default:
    case 0: xx[i] = psi(xx[i], *cc, *iipsi); break;
    case 1: xx[i] = psip(xx[i], *cc, *iipsi); break;
    }
  }  
}

void R_chifun(double *xx, double *cc, int *iipsi, int *dderiv, int *llength) {
  /* 
   * Calculate chi for vectorized x with rho(inf) = 1
   * deriv 0: chi (normalized)
   * deriv 1: chi'
   * deriv 2: chi''
   */
  
  int i;
  double nc;
  nc = normcnst(*cc, *iipsi);

  for(i = 0; i < *llength; i++) {
    switch(*dderiv) {
    default:
    case 0: xx[i] = rho(xx[i], *cc, *iipsi); break;
    case 1: xx[i] = psi(xx[i], *cc, *iipsi) * nc; break;
    case 2: xx[i] = psip(xx[i], *cc, *iipsi) * nc; break;
    }
  }  
}

void R_wgtfun(double *xx, double *cc, int *iipsi, int *llength) {
  /* 
   * Calculate wgt for vectorized x 
   */
  
  int i;

  for(i = 0; i < *llength; i++) 
    xx[i] = wgt(xx[i], *cc, *iipsi);
}

double normcnst(double c, int ipsi) {
  /* 
   * return normalizing constant for psi functions
   */

  switch(ipsi) {
  default:
  case 1: return(6./c/c); // biweight
  case 2: return(1./c/c); // GaussWeight
  case 3: return(1./3.25/c/c); // Optimal
  case 4: return(1./7.5/c/c); // Hampel
  case 5: 
    switch((int)(c*100)) { // convention: b = 1.x, x.__ eff or bp
    default:
    case 95: return(1/5.309853); break;
    case 85: return(1/2.804693); break;
    case 50: return(1/0.3748076); break;
    case 595: return(1/4.779906); break;
    case 585: return(1/2.446574); break;
    case 550: return(1/0.4007054); break;
    }
  }
}


double rho(double x, double c, int ipsi) 
{
  /* 
   * return the correct rho according to ipsi
   * this is normalized to 1
   */
  switch(ipsi) {
  default:
  case 1: return(rho_biwgt(x, c)); // biweight
  case 2: return(rho_gwgt(x, c)); // GaussWeight
  case 3: return(rho_opt(x, c)); // Optimal
  case 4: return(rho_hmpl(x, c)); // Hampel
  case 5: return(rho_gws(x, c)); // GWS
  }
}

double psi(double x, double c, int ipsi) 
{
  /* 
   * return the correct psi according to ipsi
   * this is actually rho' and not psi
   */
  switch(ipsi) {
  default:
  case 1: return(psi_biwgt(x, c)); // biweight
  case 2: return(psi_gwgt(x, c)); // GaussWeight
  case 3: return(psi_opt(x, c)); // Optimal
  case 4: return(psi_hmpl(x, c)); // Hampel
  case 5: return(psi_gws(x, c)); // GWS
  }
}

double psip(double x, double c, int ipsi) 
{
  /* 
   * return the correct ppsi according to ipsi
   * this is actually rho'' and not psip
   */
  switch(ipsi) {
  default:
  case 1: return(psip_biwgt(x, c)); // biweight
  case 2: return(psip_gwgt(x, c)); // GaussWeight
  case 3: return(psip_opt(x, c)); // Optimal
  case 4: return(psip_hmpl(x, c)); // Hampel
  case 5: return(psip_gws(x, c)); // GWS
  }
}

double wgt(double x, double c, int ipsi) 
{
  /* 
   * return the correct wgt according to ipsi
   * wgt: rho'(x) / x
   */
  switch(ipsi) {
  default:
  case 1: return(wgt_biwgt(x, c)); // biweight
  case 2: return(wgt_gwgt(x, c)); // GaussWeight
  case 3: return(wgt_opt(x, c)); // Optimal
  case 4: return(wgt_hmpl(x, c)); // Hampel
  case 5: return(wgt_gws(x, c)); // GWS
  }
}

double rho_biwgt(double x, double c)
{
/*
 * Tukey's bisquare loss function  == R's  tukeyChi()
 */
    if (fabs(x) > c)
	return(1.);
    else {
	double t;
	t = x / c;
	t *= t; /* = t^2 */
	return( t*(3. + t*(-3. + t)) );
    }
}

double psi_biwgt(double x, double c)
{
/*
 * First derivative of Tukey's bisquare loss function
 */
    if (fabs(x) > c)
	return(0.);
    else {
      double a, u;
      a = x / c;
      u = 1. - a*a;
      return( x * u * u );
    }
}

double psip_biwgt(double x, double c)
{
/*
 * Second derivative of Tukey's bisquare loss function
 */
    if (fabs(x) > c)
	return(0.);
    else {
      double x2;
      x /= c;
      x2 = x*x;
      return( (1. - x2) * (1 - 5 * x2));
    }
}

double wgt_biwgt(double x, double c)
{
/*
 * Weights for Tukey's bisquare loss function
 */
  if( fabs(x) > c )
    return(0.);
  else {
    double a;
    a = x / c;
    a = (1. - a)*(1. + a);
    return( a * a );
  }
}

double rho_gwgt(double x, double c)
{
  /* 
   * Gauss Weight Loss function
   */
  double ac;
  ac = x / c;
  ac *= ac; // ac^2
  return(-1*exp(-1*ac/2) + 1);
}

double psi_gwgt(double x, double c)
{
  /* 
   * Gauss Weight Loss function
   */
  double a, ac;
  a = x / c;
  ac = a*a; 
  return(x*exp(-1*ac/2));
}

double psip_gwgt(double x, double c)
{
  /* 
   * Gauss Weight Loss function
   */
  double ac;
  x /= c;
  ac = x*x; 
  return(exp(-1*ac/2) * (1 - ac));
}

double wgt_gwgt(double x, double c)
{
  /* 
   * Gauss Weight Loss function
   */
  double ac;
  ac = x / c;
  ac *= ac; // ac^2
  return(exp(-1*ac/2));
}

double rho_opt(double x, double c) 
{
  /* 
   * Optimal psi Function, thank you robust package
   */
  double R1 = -1.944/2., R2 = 1.728/4., R3 = -0.312/6., R4 = 0.016/8.;
  double ax, ac;
  ac = x / c; // AX=S/XK
  ax = fabs(ac); // AX=ABST/XK
  if (ax > 3) // IF (AX .GT. 3.D0) THEN
    return(1); // rlRHOm2=3.25D0*XK*XK
  else if (ax > 2.) //    ELSE IF (AX .GT. 2.D0) THEN
    return((R1*pow(ax,2)+R2*pow(ax,4)+R3*pow(ax,6)+R4*pow(ax,8)+1.792)/3.25);
	   // rlRHOm2=XK*XK*(R1*AX**2+R2*AX**4+R3*AX**6+R4*AX**8+1.792D0)
  else //      ELSE
    return(ac*ac/6.5); // rlRHOm2=S2/2.D0
}

double psi_opt(double x, double c) 
{
  /* 
   * Optimal psi Function, thank you robust package
   */
  double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016;
  double ax, ac;
  ac = x / c; // AX=S/XK
  ax = fabs(ac); // AX=ABST/XK
  if (ax > 3.) //    IF (AX .GT. 3.D0) THEN
    return(0); // rlPSIm2=0.D0
  else if (ax > 2.) { //  ELSE IF(AX .GT. 2.D0) THEN
    if (ac > 0.) //     IF (AX .GT. 0.D0) THEN
      return(fmax2(0., c*(R4*pow(ac,7)+R3*pow(ac,5)+R2*pow(ac,3)+R1*ac))); 
    // rlPSIm2=DMAX1(0.D0,XK*(R4*AX**7+R3*AX**5+R2*AX**3+R1*AX))
    else // ELSE
      return(-1*fabs(c*(R4*pow(ac,7)+R3*pow(ac,5)+R2*pow(ac,3)+R1*ac))); 
    //  rlPSIm2=-DABS(XK*(R4*AX**7+R3*AX**5+R2*AX**3+R1*AX))
  } // ENDIF
  else // ELSE
    return(x); // rlPSIm2=S 
}

double psip_opt(double x, double c) 
{
  /* 
   * Optimal psip Function, thank you robust package
   */
  double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016;
  double ax, ac;
  ac = x / c; 
  ax = fabs(ac); 
  if (ax > 3.) 
    return(0); 
  else if (ax > 2.) 
    return(7*R4*pow(ax,6)+5*R3*pow(ax,4)+3*R2*pow(ax,2)+R1); 
  else
    return(1); 
}

double wgt_opt(double x, double c) 
{
  /* 
   * Optimal psi Function, thank you robust package
   */
  double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016;
  double ax, ac;
  ac = x / c; 
  ax = fabs(ac); 
  if (ax > 3.) 
    return(0.); 
  else if (ax > 2.)
    return(fmax2(0., R4*pow(ac,6)+R3*pow(ac,4)+R2*pow(ac,2)+R1)); // /c/3.25
  else 
    return(1.); //  /c/3.25
}

double rho_hmpl(double x, double k)
{ 
  /*
   * Hampel redescending rho function
   * constants a, b, c s.t. slope of psi
   * is 1 in the center and (-)1/2 outside
   * this function is normalized s.t. rho(inf) = 1
   */
  double a, b, c, nc, u;
  a = k * 1.5;
  b = k * 3.5;
  c = k * 8.;
  nc = k * k * 7.5; 
  u = fabs(x);

  if (u <= a) 
    return( x*x/2 / nc );
  else if (u <= b)
    return( ( u - a/2 ) * a / nc );
  else if (u <= c)
    return( ( b - a/2 + (u - b) * (1 - ( u - b ) / ( c - b ) / 2 )) * a / nc);
  else 
    return( 1 );
}

double psi_hmpl(double x, double k)
{ 
  /*
   * Hampel redescending psi function
   * constants a, b, c s.t. slope of psi
   * is 1 in the center and (-)1/2 outside
   */
  double a, b, c, u;
  a = k * 1.5;
  b = k * 3.5;
  c = k * 8.;
  u = fabs(x);

  if (u <= a) 
    return( x );
  else if (u <= b) {
    if (x > 0) 
      return( a );
    else 
      return( -1*a );
  } else if (u <= c) {
    if (x > 0)
      return( a * ( c - u ) / ( c - b ) );
    else 
      return( a * ( u - c ) / ( c - b ) );
  } else 
    return( 0 );
}

double psip_hmpl(double x, double k)
{ 
  /*
   * Hampel redescending psi function
   * constants a, b, c s.t. slope of psi
   * is 1 in the center and (-)1/2 outside
   */
  double a, b, c, u;
  a = k * 1.5;
  b = k * 3.5;
  c = k * 8.;
  u = fabs(x);

  if (u <= a) 
    return( 1 );
  else if (u <= b) 
      return( 0 );
  else if (u <= c) 
    return( a / ( b - c ) );
  else 
    return( 0 );
}

double wgt_hmpl(double x, double k)
{ 
  /*
   * Hampel redescending psi function
   * constants a, b, c s.t. slope of psi
   * is 1 in the center and (-)1/2 outside
   */
  double a, b, c, u;
  a = k * 1.5;
  b = k * 3.5;
  c = k * 8.;
  u = fabs(x);

  if (u <= a) 
    return( 1 );
  else if (u <= b)
    return( a / u );
  else if (u <= c) 
    return( a * ( c - u ) / ( c - b ) / u );
  else 
    return( 0 );
}

double rho_gws(double x, double k)
{
  /*
   * Gauss Weight with constant center
   */
    const double C[6][20] = { // 0: b = 1, 95% efficiency
	{0.094164571656733, -0.168937372816728, 0.00427612218326869,
	 0.336876420549802, -0.166472338873754, 0.0436904383670537,
	 -0.00732077121233756, 0.000792550423837942, -5.08385693557726e-05,
	 1.46908724988936e-06, -0.837547853001024, 0.876392734183528,
	 -0.184600387321924, 0.0219985685280105, -0.00156403138825785,
	 6.16243137719362e-05, -7.478979895101e-07, -3.99563057938975e-08,
	 1.78125589532002e-09, -2.22317669250326e-11}, 
	// 1: b = 1, 85% efficiency
	{0.178253087655439, -0.244138991750008, 0.0404898134806095,
	 0.714486444976252, -0.507472939028585, 0.186085935389365,
	 -0.0425613907399418, 0.00616368021086328, -0.000520587101220777,
	 1.95868060120858e-05, -1.03032276702423, 1.47740024324938,
	 -0.486412165591269, 0.094417423595595, -0.0118150586827857,
	 0.00097972595308303, -5.32531744603713e-05, 1.80226431504463e-06,
	 -3.36181507658105e-08, 2.5027378473363e-10},
	// 2: b = 1, bp 0.5
	{1.33401778864039, -0.244523256745066, 0.111362610980414,
	 5.34941025041789, -10.3962071901345, 10.4296258493464,
	 -6.52564058684854, 2.58502986114483, -0.597187179950296,
	 0.0614544169812606, -1.03503907657601, 4.05896865343778,
	 -3.66686851501837, 1.95625373441388, -0.674485893541406,
	 0.154697356974158, -0.0234024029265975, 0.00222780253107263,
	 -0.000119254676258370, 2.66438378741771e-06},
	// 3: b = 1.5, 95% efficiency
	{0.104604570079252, 0.0626649856211545, -0.220058184826331,
	 0.403388189975896, -0.213020713708997, 0.102623342948069,
	 -0.0392618698058543, 0.00937878752829234, -0.00122303709506374,
	 6.70669880352453e-05, 0.632651530179424, -1.14744323908043,
	 0.981941598165897, -0.341211275272191, 0.0671272892644464,
	 -0.00826237596187364, 0.0006529134641922, -3.23468516804340e-05,
	 9.17904701930209e-07, -1.14119059405971e-08},
	// 4: b = 1.5, 85% efficiency
	{0.204367401015115, 0.150893608928378, -0.577989001048967,
	 1.0294512751671, -0.575516566945057, 0.259880185675216,
	 -0.109860846048747, 0.0323047137569345, -0.00523755203714339,
	 0.000351158243420723, 1.16858313773729, -2.94521801333508,
	 3.31178336703781, -1.69454782208968, 0.505158881547528,
	 -0.0954552533042206, 0.0116753187335009, -0.00090089548565088,
	 4.00306861156273e-05, -7.83125565664021e-07},
	// 5: b = 1.5, bp 0.5
	{1.24779939864579, 0.150908420620023, -1.42828766163334,
	 6.28563903546494, -8.68253013847564, 9.68713466780148,
	 -10.1184959700703, 7.35195633088053, -2.94530819532254,
	 0.487942915516091, 1.17510828685948, -7.31091885661658,
	 20.2953910903282, -25.6614599466583, 18.9099061867153,
	 -8.83450648767371, 2.67203673126514, -0.509920849882194,
	 0.0560449261448248, -0.00271236230476857}};
    const double end[6] = {18.5527638190955, 12.7638190954774, 
			   .71356783919598, 
			   11.4974874371859, 7.42713567839196, 
			   2.99497487437186};
    int j;
    double c;
    switch((int)(k*100)) { // convention: b = 1.x, x.__ eff or bp
    default:
    case 95:  j = 0; c = 1.694;     break;
    case 85:  j = 1; c = 1.319;     break;
    case 50:  j = 2; c = 0.482284;  break;
    case 595: j = 3; c = 1.063;     break;
    case 585: j = 4; c = 0.941;     break;
    case 550: j = 5; c = 0.3808319; break;
    }
    x = fabs(x);
    if (x <= c)
	return(C[j][0]*x*x);
    else if (x <= 3*c)
	return(C[j][1] + C[j][2]*x + C[j][3]*x*x + C[j][4]*pow(x, 3) + 
	       C[j][5]*pow(x, 4) + C[j][6]*pow(x, 5) + C[j][7]*pow(x, 6) + 
	       C[j][8]*pow(x, 7) + C[j][9]*pow(x, 8));
    else if (x <= end[j])
	return(C[j][10] + C[j][11]*x + C[j][12]*x*x + C[j][13]*pow(x, 3) + 
	       C[j][14]*pow(x, 4) + C[j][15]*pow(x, 5) + 
	       C[j][16]*pow(x, 6) + C[j][17]*pow(x, 7) + 
	       C[j][18]*pow(x, 8) + C[j][19]*pow(x, 9));
    else return(1);
}

double psi_gws(double x, double k)
{
  /*
   * Gauss Weight with constant center
   */
  // set a,b,c
  double a, b, c;
  switch((int)(k*100)) { // convention: b = 1.x, x.__ eff or bp
  default:
  case 95: a = 0.648; b = 1; c = 1.694; break;
  case 85: a = 0.440; b = 1; c = 1.319; break;
  case 50: a = 0.1607910; b = 1; c = 0.482284; break;
  case 595: a = 1.387; b = 1.5; c = 1.063; break;
  case 585: a = 0.703; b = 1.5; c = 0.941; break;
  case 550: a = 0.1809862;  b = 1.5;  c = 0.3808319; break;
  }
  double ax;
  ax = fabs(x);
  if (ax < c) 
    return(x);
  else 
    return(x*exp(-1*pow(ax-c,b)/2/a));
}

double psip_gws(double x, double k)
{
  /*
   * Gauss Weight with constant center
   */
  // set a,b,c
  double a, b, c;
  switch((int)(k*100)) { // convention: b = 1.x, x.__ eff or bp
  default:
  case 95: a = 0.648; b = 1; c = 1.694; break;
  case 85: a = 0.440; b = 1; c = 1.319; break;
  case 50: a = 0.1607910; b = 1; c = 0.482284; break;
  case 595: a = 1.387; b = 1.5; c = 1.063; break;
  case 585: a = 0.703; b = 1.5; c = 0.941; break;
  case 550: a = 0.1809862;  b = 1.5;  c = 0.3808319; break;
  }
  double ax;
  ax = fabs(x);
  if (ax < c)
    return(1);
  else
    return(exp(-1*pow(ax-c,b)/2/a)*(1-b/2/a*ax*pow(ax-c,b-1)));
}

double wgt_gws(double x, double k)
{
  /*
   * Gauss Weight with constant center
   */
  // set a,b,c
  double a, b, c;
  switch((int)(k*100)) { // convention: b = 1.x, x.__ eff or bp
  default:
  case 95: a = 0.648; b = 1; c = 1.694; break;
  case 85: a = 0.440; b = 1; c = 1.319; break;
  case 50: a = 0.1607910; b = 1; c = 0.482284; break;
  case 595: a = 1.387; b = 1.5; c = 1.063; break;
  case 585: a = 0.703; b = 1.5; c = 0.941; break;
  case 550: a = 0.1809862;  b = 1.5;  c = 0.3808319; break;
  }
  double ax;
  ax = fabs(x);
  if (ax < c)
    return(1);
  else 
    return(exp(-1*pow(ax-c,b)/2/a));
}

/* this function finds the k-th place in the
 * vector a, in the process it permutes the
 * elements of a
 */
double kthplace(double *a, int n, int k)
{
    int jnc,j;
    int l,lr;
    double ax,w;
    k--;
    l=0;
    lr=n-1;
    while (l < lr) {
	ax=a[k];
	jnc=l;
	j=lr;
	while (jnc <= j) {
	    while (a[jnc] < ax) jnc++;
	    while (a[j] > ax) j--;
	    if (jnc <= j) {
		w=a[jnc];
		a[jnc]=a[j];
		a[j]=w;
		jnc++;
		j--;
	    }
	}
	if (j < k) l=jnc;
	if (k < jnc) lr=j;
    }
    return(a[k]);
}

#ifdef _UNUSED_
void sampler_i(int n, int *x)
{
/* function to get a random sample of
 * indices (0 to n-1)
 * *x receives the output
 * rand() returns an integer between 0 and RAND_MAX
 */
    int i;
    for(i=0; i < n; i++)
	x[i] = (int) ( (double) rand() / RAND_MAX * (double) (n-1) );
}
#endif

/* This is from VR's bundle, MASS package  VR/MASS/src/lqs.c : */
/*
   Sampling k from 0:n-1 without replacement.
 */
static void sample_noreplace(int *x, int n, int k, int *ind_space)
{
    int i, j, nn=n;
#define II ind_space

    for (i = 0; i < n; i++) II[i] = i;
    for (i = 0; i < k; i++) {
	j = nn * unif_rand();
	x[i] = II[j];
	II[j] = II[--nn];
    }
#undef II
}

#ifdef _NO_LONGUER_USED_
void sample_n_outof_N(int n, int N, int *x)
{
/* function to get a random sample of size n
 * of the indices (0 to N) WITHOUT replication
 * *x receives the output
 */
    int i;
    if( N < n ) {
	/* printf("\nCant get %d out of %d \ without replication\n", n, N); */
	for(i=0; i < n; i++) x[i] = i;
    } else {
	int j, is_previous, cand=0;
	for(i=0; i < n; i++) {
	    is_previous=1;
	    while (is_previous) {
		is_previous = 0;
		cand = N * unif_rand();
		/* cand = (int) ( (double) rand() / RAND_MAX * (double) N ); */
		for(j=0; j < i; j++) {
		    if(cand == x[j]) {
			is_previous = 1; break;
		    }
		}
	    }
	    x[i]=cand;
	}
    }
}
#endif

/* C = A + B */
void sum_mat(double **a, double **b, double **c, int n, int m)
{
    int i,j;
    for(i=0; i < n; i++)
	for(j=0; j < m; j++)
	    c[i][j] = a[i][j] + b[i][j];
}

/* A = v1 %*% t(v2) */
void outer_vec_vec(double **a, double *v1, double *v2, int n)
{
    int i,j;
    for(i=0; i < n; i++)
	for(j=0; j < n; j++)
	    a[i][j] = v1[i] * v2[j];
}

/* C = A * b */
void scalar_mat(double **a, double b, double **c, int n, int m)
{
    int i,j;
    for(i=0; i < n; i++)
	for(j=0; j < m; j++)
	    c[i][j]  = b * a[i][j];
}

/* c = a * b */
void scalar_vec(double *a, double b, double *c, int n)
{
    int i;
    for(i=0; i < n; i++)
	c[i]  = b * a[i];
}

/* returns the inner product of a and b, i.e. t(a) %*% b */
double vecprime_vec(double *a, double *b, int n)
{
    int i;
    double s = 0.;
    for(i=0; i < n; i++) s += a[i] * b[i];
    return(s);
}

/* c = a + b */
void sum_vec(double *a, double *b, double *c, int n)
{
    int i;
    for(i=0; i < n; i++) c[i] = a[i] + b[i];
}

/* c = a - b */
void dif_vec(double *a, double *b, double *c, int n)
{
    int i;
    for(i=0; i < n; i++) c[i] = a[i] - b[i];
}

/* c = a / b */
void div_vec(double *a, double *b, double *c, int n)
{
    int i;
    for(i=0; i < n; i++) c[i] = a[i] / b[i];
}


/* C = A - B */
void dif_mat(double **a, double **b, double **c, int n, int m)
{
    int i,j;
    for(i=0; i < n; i++)
	for(j=0; j < m; j++) c[i][j] = a[i][j] - b[i][j];
}

/* c = A %*% b */
void mat_vec(double **a, double *b, double *c, int n, int m)
{
    int i,j;
    for(i=0; i < n; i++)
	for(c[i]=0,j=0; j < m; j++) c[i] += a[i][j] * b[j];
}

/* C = A %*% B */
void mat_mat(double **a, double **b, double **c, int n,
		int m, int l)
{
    int i,j,k;
    for(i=0; i < n; i++)
	for(j=0; j < l; j++) {
	    c[i][j] = 0;
	    for(k=0; k < m; k++) c[i][j] += a[i][k] * b[k][j];
	}
}

/* RWLS iterations starting from i_estimate,
 * ---- the workhorse of the "lmrob_MM" algorithm;
 * in itself,  ``just'' an M-estimator :
 *
 * FIXME: rather than using C-matrices  **a, **b, and copying
 *        should use BLAS and Lapack
 */
int rwls(double **a, int n, int p,
	 double *estimate, double *i_estimate,
	 double *resid, double* loss,
	 double scale, double epsilon,
	 int *max_it, /* on Input:  maximal number of iterations;
			 on Output: number of iterations */
	 const double rho_c, const int ipsi, int trace_lev)
{
    double **b, *beta1, *beta2, *beta0, *weights;
    double d_beta = 0.;
    int i,j, iterations = 0;
    Rboolean converged = FALSE;

    if ( (b = (double **) Calloc(p, double *)) == NULL)
	return(1);
    for (i=0; i < p; i++)
	if ( (b[i] = (double *) Calloc((p+1), double)) == NULL)
	    return(1);/* FIXME: memory leaks if this happens! */
    beta1 = (double *) Calloc(p, double);
    beta2 = (double *) Calloc(p, double);
    beta0 = (double *) Calloc(p, double);
    weights = (double *) Calloc(n, double);

    COPY_beta(i_estimate, beta1);

    /* main loop */
    while(!converged &&	 ++iterations < *max_it) {

	double r,s;
	int k;
#ifdef LAMBDA_Iterations_remnant_no_longer_needed
	double loss2, lambda;
	int iter_lambda;
#endif

	R_CheckUserInterrupt();
	for(i=0; i < n; i++) {
	    s=0;
	    for(j=0; j < p; j++)
		s += a[i][j] * beta1[j];
	    r = a[i][p]- s;
	    weights[i] = wgt(r/scale, rho_c, ipsi);
	    /* if(fabs(r/scale) < 1e-7) */
	    /* 	weights[i] = 1.0 / scale / rho_c; */
	    /* else */
	    /* 	weights[i] = psi_biwgt(r/scale, rho_c) / r; */
	}

	COPY_beta(beta1, beta2);

	/* get the residuals and loss for beta2 */
	for(i=0; i < n; i++) {
	    s = 0;
	    for(j=0; j < p; j++)
		s += a[i][j] * beta2[j];
	    resid[i] = a[i][p] - s;
	}
#ifdef LAMBDA_Iterations_remnant_no_longer_needed
	loss2 = sum_rho(resid,n,rho_c,ipsi);
#endif
	if(trace_lev >= 2) {
#ifdef LAMBDA_Iterations_remnant_no_longer_needed
	    *loss = loss2;
#else
	    *loss = sum_rho(resid,n,rho_c,ipsi);
#endif
	    Rprintf(" it %4d: L(b2) = %12g ", iterations, *loss);
	}

	/* S version of the following code:
	 * A <- matrix(0, p, p)
	 * Y <- rep(0, p)
	 * for(i in 1:n) {
	 *   A <- A + w[i] * a[i,0:(p-1)] %*% t(a[i,0:(p-1)])
	 *   Y <- Y + w[i] * a[i,0:(p-1)] * a[i,p]
	 * }
	 * beta1 <- solve(A, Y)
	 */
	for(j=0; j < p; j++)
	    for(k=0; k <= p; k++) {
		b[j][k]=0.;
		for(i=0; i < n; i++)
		    b[j][k] += a[i][j] * a[i][k] * weights[i];
	    }

	if( lu(b,&p,beta1) ) { /* is system singular? */
	    /* converged = FALSE; */
	    if(trace_lev)
		Rprintf("rwls(.): stopping early because of singular LU decomposition\n");
	    break; /* out of while(.) */
	}

	/* else -- lu() was not singular : */

	/* get the residuals and loss for beta1 */
	for(i=0; i < n; i++) {
	    s = 0;
	    for(j=0; j < p; j++)
		s += a[i][j] * beta1[j];
	    resid[i] = a[i][p] - s;
	}
	*loss = sum_rho(resid,n,rho_c,ipsi);

#ifdef LAMBDA_Iterations_remnant_no_longer_needed
	/* Is beta1 good enough? --
	 * if not, find beta0 in [beta1, beta2] and use that. */

	COPY_beta(beta1, beta0);
	/* from now on:	 *loss = loss( beta0 ) ... */
	lambda = 1.;
	iter_lambda = 1;
	if(trace_lev >= 3)
	    Rprintf(" -- lambda iterations: L(b0) = %12g\n%s\n",
		    *loss, "  l.it |      lambda ,   L(<beta0>)");
	while( *loss > loss2 ) {
	    lambda /= 2.;
	    for(j=0; j < p; j++)
		beta0[j] = (1 - lambda) * beta2[j] + lambda * beta1[j];
	    /* get the residuals and loss for beta0 */
	    for(i=0; i < n; i++) {
		s = 0;
		for(j=0; j < p; j++)
		    s += a[i][j] * beta0[j];
		resid[i] = a[i][p] - s;
	    }
	    *loss = sum_rho(resid,n,rho_c,ipsi);

	    if(trace_lev >= 3 && iter_lambda % 20 == 1)
		Rprintf("  %4d | %11g , %12g\n", iter_lambda, lambda, *loss);

	    if( ++iter_lambda > 1000) {
		warning("rwls(): not converged in 1000 lambda iterations");
		COPY_beta(beta2, beta0);
		/* FIXME - should we signal non-convergence ?
		 * converged = FALSE; */
		break; /* out of while(.) lambda iterations */
	    }
	} /* end while (*loss > loss2) */
	if(trace_lev >= 2)
	    Rprintf(" used %4d lambda iter.;", iter_lambda);

#endif /*lambda iterations */

	d_beta = norm1_diff(beta1,beta2, p);
	if(trace_lev >= 2) {
	    if(trace_lev >= 3) {
		Rprintf("\n  b2 = (");
		for(j=0; j < p; j++)
		    Rprintf("%s%11g", (j > 0)? ", " : "", beta0[j]);
		Rprintf(");");
	    }
	    Rprintf(" ||b1 - b2||_1 = %g\n", d_beta);
	}


#ifdef LAMBDA_Iterations_remnant_no_longer_needed
	converged = d_beta <= epsilon * fmax2(epsilon, norm1(beta0, p));
#else
	converged = d_beta <= epsilon * fmax2(epsilon, norm1(beta1, p));
#endif

    } /* end while(!converged & iter <=...) */

#ifdef LAMBDA_Iterations_remnant_no_longer_needed
    COPY_beta(beta0, estimate);
#else
    COPY_beta(beta1, estimate);
#endif

    if(trace_lev)
	Rprintf(" rwls() used %d it.; last ||b1 - b2||_1 = %g;%sconvergence\n",
		iterations, d_beta, (converged ? " " : " NON-"));

    Free(weights); Free(beta1); Free(beta2); Free(beta0);
    for(i=0; i < p; i++) Free(b[i]);
    Free(b);

    *max_it = iterations;

    return (int)converged;

} /* rwls() */

/* sets the entries of a matrix to zero */
void zero_mat(double **a, int n, int m)
{
    int i,j;
    for(i=0; i < n; i++)
	for(j=0; j < m; j++)
	    a[i][j] = 0.;
}

/* sets the entries of a vector to zero */
void zero_vec(double *a, int n)
{
    int i;
    for(i=0; i < n; i++) a[i] = 0.;
}


/*
 *
 * 2004 / 5 -- Matias Salibian-Barrera & Victor Yohai
 * Department of Statistics, University of British Columbia
 * matias@stat.ubc.ca
 * Department of Mathematics, University of Buenos Aires
 * vyohai@uolsinectis.com.ar
 *
 *
 * Reference: A fast algorithm for S-regression estimates,
 * 2005, Salibian-Barrera and Yohai.
*/

/* This function implements the "large n" strategy
*/
void fast_s_large_n(double *X, double *y,
		    int *nn, int *pp, int *nRes,
		    int *ggroups, int *nn_group,
		    int *K, int *max_k, double *rel_tol, int *converged,
		    int *best_r, double *bb, double *rrhoc, int *iipsi,
		    double *bbeta, double *sscale,
		    int *trace_lev)
{
/* *X  = the n x p  design matrix (incl. intercept if appropriate),
 *	 in column order as sent by R)
 * *y  = the ( n ) response vector
 * *nn =: n = the length of y
 * *pp =: p = the number of columns in X
 * *nRes  = number of re-sampling candidates to be used in each partition
 * *K     = number of refining steps for each candidate (typically 1 or 2)
 * *max_k = number of refining steps for each candidate (typically 1 or 2)
             [used to be hard coded to MAX_ITER_REFINE_S = 50 ]
 * *rel_tol= convergence tolerance for iterative refinement iterations
             [used to be hard coded to EPS = 1e-7 ]
 * *converged: will become 0(FALSE)  iff at least one of the best_r iterations
 *             did not converge (in max_k steps to rel_tol precision)
 * *ggroups = number of groups in which to split the
 *	      random subsample
 * *nn_group = size of each of the (*ggroups) groups
 *	       to use in the random subsample
 * *best_r = no. of best candidates to be iterated further ("refined")
 * *bb	   = right-hand side of S-equation (typically 1/2)
 * *rrhoc  = tuning constant for loss function
 *	     (this should be associated with *bb)
 * *iipsi  = indicator for type of psi function to be used
 * *bbeta  = final estimator
 * *sscale = associated scale estimator (or -1 when problem)
*/
    int i,j,k, k2, it_k;
    int n = *nn, p = *pp, kk = *K;
    int groups = *ggroups, n_group = *nn_group, ipsi = *iipsi;
    double b = *bb, rhoc = *rrhoc, sc, best_sc, worst_sc;
    int pos_worst_scale;
    Rboolean conv;

    /* (Pointers to) Arrays - to be allocated */
    int *ind_space, *indices;
    double **best_betas, *best_scales, *weights;
    double *tmp, *tmp2, **tmp_mat, **tmp_mat2;
    double **x, **xsamp, *ysamp, *res, *beta_ref;
    double **final_best_betas, *final_best_scales;

    SETUP_FAST_S(n, p);

    beta_ref = (double *) Calloc(p, double);
    final_best_betas = (double **) Calloc(*best_r, double *);
    for(i=0; i < *best_r; i++)
	final_best_betas[i] = (double *) Calloc(p, double);
    final_best_scales = (double *) Calloc(*best_r, double);
    k = *best_r * groups;
    best_betas = (double **) Calloc(k,	double *);
    best_scales = (double *) Calloc(k,	double );
    for(i=0; i < k; i++)
	best_betas[i] = (double*) Calloc(p, double);
    x = (double**) Calloc(n, double *);
    for(i=0; i < n; i++)
	x[i] = (double*) Calloc(p, double);

    k = n_group * groups;
    indices = (int *) Calloc(k, int);
    xsamp = (double**) Calloc(k, double *);
    ysamp = (double*) Calloc(k, double);
    for(i=0; i < k; i++)
	xsamp[i] = (double*) Calloc(p, double);
    for(i=0; i < n; i++)
	for(j=0; j < p; j++)
	    x[i][j] = X[j*n+i];

    /* assume that n > 2000 */

    /*	set the seed */
    /* srand((long)*seed_rand); */
    GetRNGstate();

    /* get a sample of k indices */
    /* sample_n_outof_N(k, n-1, indices); */
    sample_noreplace(indices, n, k, ind_space);
    /* FIXME: Also look at lqs_setup(),
     * -----  and  xr[.,.] "fortran-like" matrix can be used from there!*/

    /* get the sampled design matrix and response */
    for(i=0; i < k; i++) {
	for(j=0; j < p; j++)
	    xsamp[i][j] = x[indices[i]][j];
	ysamp[i] = y[indices[i]];
    }

/* For each (of 'groups') group : get the *best_r best betas : */

    for(i=0; i < groups; i++) {
	if( !fast_s_with_memory(xsamp+i*n_group, ysamp+i*n_group,
				&n_group, pp, nRes, K, max_k, rel_tol,
				trace_lev, best_r, ind_space, bb, rrhoc,
				iipsi, best_betas + i* *best_r,
				best_scales+ i* *best_r)) {

	    *sscale = -1.; /* problem */
	    goto cleanup_and_return;
	}
    }
/* now	iterate (refine) these "best_r * groups"
 * best betas in the (xsamp,ysamp) sample
 * with kk C-steps and keep only the "best_r" best ones
*/
    conv = FALSE;
    pos_worst_scale = 0;
    for(i=0; i < *best_r; i++)
	final_best_scales[i] = INFI;
    worst_sc = INFI;
    /* set the matrix to zero */
    zero_mat(final_best_betas, *best_r, p);
    k = n_group * groups;
    for(i=0; i < (*best_r * groups); i++) {
	if(*trace_lev >= 3) {
	    Rprintf("sample[%3d]: before refine_(*, conv=FALSE):\n", i);
	    Rprintf("beta_cand : "); disp_vec(best_betas[i],p);
	    Rprintf("with scale %.15g\n", best_scales[i]);
	    
	}
	refine_fast_s(xsamp, ysamp, weights, k, p, res,
		      tmp, tmp2, tmp_mat, tmp_mat2, best_betas[i],
		      kk, &conv/* = FALSE*/, *max_k, *rel_tol, *trace_lev,
		      b, rhoc, ipsi, best_scales[i], /* -> */ beta_ref, &sc);
	if(*trace_lev >= 3) {
	  Rprintf("after refine: beta_ref : "); disp_vec(beta_ref,p);
	  Rprintf("with scale %.15g\n", sc);
	}
	if ( sum_rho_sc(res, worst_sc, k, p, rhoc, ipsi) < b ) {
	    /* scale will be better */
	  sc = find_scale(res, b, rhoc, ipsi, sc, k, p);
	    k2 = pos_worst_scale;
	    final_best_scales[ k2 ] = sc;
	    COPY_beta(beta_ref, final_best_betas[k2]);
	    pos_worst_scale = find_max(final_best_scales, *best_r);
	    worst_sc = final_best_scales[pos_worst_scale];
	}
    }

/* now iterate the best "best_r"
 * betas in the whole sample until convergence (max_k, rel_tol)
*/
    best_sc = INFI; *converged = 1;  k = 0;

    for(i=0; i < *best_r; i++) {
	conv = TRUE;
	it_k = refine_fast_s(x, y, weights, n, p, res, tmp, tmp2,
			     tmp_mat, tmp_mat2, final_best_betas[i], kk,
			     &conv/* = TRUE */, *max_k, *rel_tol, *trace_lev,
			     b, rhoc, ipsi, final_best_scales[i],
			     /* -> */ beta_ref, &sc);
	if(*trace_lev)
	    Rprintf(" best [%d]: %d iterations, i.e. %sconverged\n",
		    i, it_k, conv ? " " : "NOT ");
	if(best_sc > sc) {
	    best_sc = sc;
	    COPY_beta(beta_ref, bbeta);
	}
	if (!conv && *converged) *converged = 0;
	if (k < it_k) k = it_k;
    }
    *sscale = best_sc;
    *max_k = k;

/* Done. Now clean-up. */

cleanup_and_return:
    PutRNGstate();

    for(i=0; i < n; i++) Free(x[i]); Free(x);
    Free(best_scales);
    k = *best_r * groups;
    for(i=0; i < k; i++) Free( best_betas[i] );
    Free(best_betas); Free(indices); Free(ysamp);
    k = n_group * groups;
    for(i=0; i < k; i++) Free(xsamp[i]);
    Free(xsamp);
    for(i=0; i < p; i++) {
	Free(tmp_mat[i]);
	Free(tmp_mat2[i]);
    }
    Free(tmp_mat); Free(tmp_mat2);
    for(i=0; i < *best_r; i++)
	Free(final_best_betas[i]);
    Free(final_best_betas);
    Free(final_best_scales);
    Free(beta_ref);

} /* fast_s_large_n() */

int fast_s_with_memory(double **x, double *y,
		       int *nn, int *pp, int *nRes,
		       int *K, int *max_k, double *rel_tol, int *trace_lev,
		       int *best_r, int *ind_space, double *bb, double *rrhoc, int *iipsi,
		       double **best_betas, double *best_scales)
{
/*
 * Called from fast_s_large_n(), the adjustment for large "n",
 * same as fast_s, but it returns the best_r best betas,
 * and their associated scales.
 *
 * x an	 n x p design matrix (including intercept if appropriate)
 * y an	 n vector
 * *nn = n, *pp = p
 * *nRes   = number of re-sampling candidates to be taken
 * *K	   = number of refining steps for each candidate
 * *best_r = number of (refined) to be retained for full iteration
 * *bb	   = right-hand side of the S-equation (typically 1/2)
 * *rrhoc  = tuning constant for loss function
 *	     (this should be associated with *bb)
 * *iipsi  = indicator for type of loss function to be used
 * *best_betas	= returning the best ... coefficient vectors
 * *best_scales = returning their associated residual scales
*/
    int i,j,k, no_try_samples;
    int n = *nn, p = *pp, nResample = *nRes, *b_i;
    Rboolean lu_sing, conv = FALSE;
    double **x_samp, *beta_cand, *beta_ref, *res;
    int ipsi = *iipsi;
    double b = *bb, rhoc = *rrhoc, sc, worst_sc = INFI;
    double *weights;
    int pos_worst_scale, is_ok = 1;
    double *tmp, *tmp2, **tmp_mat2, **tmp_mat;

    res	      = (double *) Calloc(n, double);
    tmp	      = (double *) Calloc(n, double);
    tmp2      = (double *) Calloc(n, double);
    weights   = (double *) Calloc(n, double);
    beta_cand = (double *) Calloc(p, double);
    beta_ref  = (double *) Calloc(p, double);
    b_i	      = (int *) Calloc(n, int);
    x_samp    = (double **) Calloc(n, double *);
    tmp_mat   = (double **) Calloc(p, double *);
    tmp_mat2  = (double **) Calloc(p, double *);
    for(i=0; i < n; i++)
	x_samp[i]  = (double *) Calloc((p+1), double);
    for(i=0; i < p; i++) {
	tmp_mat [i] = (double *) Calloc(p, double);
	tmp_mat2[i] = (double *) Calloc((p+1), double);
    }

    for(i=0; i < *best_r; i++)
	best_scales[i] = INFI;
    pos_worst_scale = 0;

/* resampling approximation  */

    for(i=0; i < nResample; i++) {

	/* find a candidate */
	no_try_samples = 0;
	do {
	    R_CheckUserInterrupt();
	    if( (++no_try_samples) > MAX_NO_TRY_SAMPLES ) {
		REprintf("\nToo many singular resamples\n"
			 "Aborting fast_s_w_mem()\n\n");
		is_ok = 0;
		goto cleanup_and_return;
	    }
	    /* take a sample of the indices  */
	    /* sample_n_outof_N(p, n-1, b_i); */
	    sample_noreplace(b_i, n, p, ind_space);

	    /* build the submatrix */
	    for(j=0; j < p; j++) {
		for(k=0; k < p; k++)
		    x_samp[j][k]=x[b_i[j]][k];
		x_samp[j][p]=y[b_i[j]];
	    }
	    /* solve the system, lu() = TRUE means matrix is singular
	     */
	    lu_sing = lu(x_samp,pp,beta_cand);

	} while(lu_sing);

	/* disp_vec(beta_cand,p); */

	/* improve the re-sampling candidate */

	/* conv = FALSE : do *K refining steps */
	refine_fast_s(x, y, weights, n, p, res,
		      tmp, tmp2, tmp_mat, tmp_mat2, beta_cand,
		      *K, &conv/* = FALSE*/, *max_k, *rel_tol, *trace_lev,
		      b, rhoc, ipsi, -1., /* -> */ beta_ref, &sc);

	/* FIXME: if sc ~ 0 ---> return beta_cand and be done */

	if ( sum_rho_sc(res, worst_sc, n, p, rhoc, ipsi) < b )	{
	    /* scale will be better */
	  sc = find_scale(res, b, rhoc, ipsi, sc, n, p);
	    k = pos_worst_scale;
	    best_scales[ k ] = sc;
	    for(j=0; j < p; j++)
		best_betas[k][j] = beta_ref[j];
	    pos_worst_scale = find_max(best_scales, *best_r);
	    worst_sc = best_scales[pos_worst_scale];
	}

    } /* for(i ) */

cleanup_and_return:
    Free(tmp); Free(tmp2);
    Free(res); Free(weights); Free(beta_cand);
    Free(beta_ref); Free(b_i);
    for(i=0; i < n; i++) {
	Free(x_samp[i]);
    }
    for(i=0; i < p; i++) {
	Free(tmp_mat[i]);
	Free(tmp_mat2[i]);
    }
    Free(x_samp); Free(tmp_mat); Free(tmp_mat2);

    return is_ok;
} /* fast_s_with_memory() */

void fast_s(double *X, double *y,
	    int *nn, int *pp, int *nRes,
	    int *K, int *max_k, double *rel_tol, int *converged,
	    int *best_r, double *bb, double *rrhoc, int *iipsi,
	    double *bbeta, double *sscale, int *trace_lev)
{
/* *X  = the n x p  design matrix (incl. intercept if appropriate),
 *	 in column order as sent by R)
 * *y  = the ( n ) response vector
 * *nn =: n = the length of y
 * *pp =: p = the number of columns in X
 * *nRes   = number of re-sampling candidates to be taken
 * *K	   = number of refining steps for each candidate
 * *best_r = number of (refined) to be retained for full iteration
 * *converged: will become FALSE  iff at least one of the best_r iterations
 *	       did not converge (in max_k steps to rel_tol precision)
 * *bb	   = right-hand side of the S-equation (typically 1/2)
 * *rrhoc  = tuning constant for loss function
 *	     (this should be associated with *bb)
 * *iipsi  = indicator for type of loss function to be used
 * *bbeta  = final estimator
 * *sscale = associated scale estimator (or -1 when problem)
*/
    int i,j,k, it_k, no_try_samples;
    int n = *nn, p = *pp, nResample = *nRes, ipsi = *iipsi;
    Rboolean lu_sing, conv;
    double b = *bb, rhoc = *rrhoc;
    double sc, best_sc, worst_sc, aux;
    int pos_worst_scale;

    /* Rprintf("fast_s %d\n", ipsi); */

    /* (Pointers to) Arrays - to be allocated */
    int *ind_space, *b_i;
    double **x, **x_samp, *beta_cand, *beta_ref, *res;
    double **best_betas, *best_scales, *weights;
    double *tmp, *tmp2, **tmp_mat, **tmp_mat2;

    SETUP_FAST_S(n, p);

    best_betas = (double **) Calloc(*best_r, double *);
    best_scales = (double *) Calloc(*best_r, double);
    for(i=0; i < *best_r; i++) {
	best_betas[i] = (double*) Calloc(p, double);
	best_scales[i] = INFI;
    }
    beta_cand = (double *) Calloc(p, double);
    beta_ref  = (double *) Calloc(p, double);
    b_i	  = (int *) Calloc(n, int);
    x	  = (double **) Calloc(n, double *);
    x_samp= (double **) Calloc(n, double *);
    for(i=0; i < n; i++) {
	x[i]	   = (double *) Calloc(p, double);
	x_samp[i]  = (double *) Calloc((p+1), double);
    }
    for(i=0; i < n; i++)
	for(j=0; j < p; j++)
	    x[i][j] = X[j*n+i];

    /* disp_mat(x, n, p); */

    pos_worst_scale = 0; conv = FALSE; worst_sc = INFI;

    /* srand((long)*seed_rand); */
    GetRNGstate();

/* resampling approximation  */

    for(i=0; i < nResample; i++) {

	/* find a candidate */
	no_try_samples = 0;
	do {
	    R_CheckUserInterrupt();
	    if( (++no_try_samples) > MAX_NO_TRY_SAMPLES ) {
		REprintf("\nToo many singular resamples\n"
			 "Aborting fast_s()\n\n");
		*sscale = -1.;
		goto cleanup_and_return;
	    }
	    /* take a sample of the indices  */
	    /* sample_n_outof_N(p, n-1, b_i); */
	    sample_noreplace(b_i, n, p, ind_space);

	    /* build the submatrix */
	    for(j=0; j < p; j++) {
		for(k=0; k < p; k++)
		    x_samp[j][k] = x[b_i[j]][k];
		x_samp[j][p] = y[b_i[j]];
	    }
	    /* solve the system, lu() = TRUE means matrix is singular
	     */
	    lu_sing = lu(x_samp,pp,beta_cand);

	} while(lu_sing);

	/* disp_vec(beta_cand,p); */

	/* improve the re-sampling candidate */

	/* conv = FALSE : do *k refining steps */
	refine_fast_s(x, y, weights, n, p, res,
		      tmp, tmp2, tmp_mat, tmp_mat2, beta_cand,
		      *K, &conv/* = FALSE*/, *max_k, *rel_tol, *trace_lev,
		      b, rhoc, ipsi, -1., /* -> */ beta_ref, &sc);
	if(*trace_lev >= 2) {
	    double del = norm_diff(beta_cand, beta_ref, p);
	    Rprintf("sample[%3d]: after refine_(*, conv=FALSE):\n", i);
	    Rprintf("beta_ref : "); disp_vec(beta_ref,p);
	    Rprintf(" with ||beta_ref - beta_cand|| = %.12g, --> sc = %.15g\n",
		    del, sc);
	}
	if(fabs(sc) == 0.) { /* exact zero set by refine_*() */
	    if(*trace_lev >= 1)
		Rprintf("too many exact zeroes -> leaving refinement!\n");
	    *sscale = sc;
	    COPY_beta(beta_cand, bbeta);
	    goto cleanup_and_return;
	}
	if ( sum_rho_sc(res, worst_sc, n, p, rhoc, ipsi) < b )	{
	    /* scale will be better */
	  sc = find_scale(res, b, rhoc, ipsi, sc, n, p);
	    k = pos_worst_scale;
	    best_scales[ k ] = sc;
	    COPY_beta(beta_ref, best_betas[k]);
	    pos_worst_scale = find_max(best_scales, *best_r);
	    worst_sc = best_scales[pos_worst_scale];
	}

    } /* for(i ) */

/* now look for the very best */
    if(*trace_lev)
	Rprintf("now refine() to convergence for %d very best ones:\n",
		*best_r);

    best_sc = INFI; *converged = 1;  k = 0;
    for(i=0; i < *best_r; i++) {
	conv = TRUE;
	it_k = refine_fast_s(x, y, weights, n, p, res, tmp, tmp2,
			     tmp_mat, tmp_mat2, best_betas[i], *K,
			     &conv /* = TRUE */, *max_k, *rel_tol, *trace_lev,
			     b, rhoc, ipsi, best_scales[i], /* -> */ beta_ref, &aux);
	if(*trace_lev)
	    Rprintf("i=%2d: %sconvergence (%d iter.):",
		    i, (conv) ? "" : "NON ", it_k);
	if(aux < best_sc) {
	    if(*trace_lev)
		Rprintf(" -> improved scale to %.15g", aux);
	    best_sc = aux;
	    COPY_beta(beta_ref, bbeta);
	}
	if(*trace_lev) Rprintf("\n");
	if (!conv && *converged) *converged = 0;
	if (k < it_k) k = it_k;
    }
    *sscale = best_sc;
    *max_k = k;

cleanup_and_return:

    PutRNGstate();

    Free(best_scales);
    Free(beta_cand);
    Free(beta_ref); Free(b_i);
    for(i=0; i < *best_r; i++)
	Free(best_betas[i]);
    Free(best_betas);
    for(i=0; i < n; i++) {
	Free(x[i]);
	Free(x_samp[i]);
    }
    Free(x); Free(x_samp);
    for(i=0; i < p; i++) {
	Free(tmp_mat[i]);
	Free(tmp_mat2[i]);
    }
    Free(tmp_mat); Free(tmp_mat2);

    return;
} /* fast_s() */

int refine_fast_s(double **x, double *y, double *weights,
		  int n, int p, double *res,
		  double *tmp, double *tmp2,
		  double **tmp_mat, double **tmp_mat2, double *beta_cand,
		  int kk, Rboolean *conv, int max_k, double rel_tol,
		  int trace_lev,
		  double b, double rhoc, int ipsi, double initial_scale,
		  double *beta_ref, double *scale)
{
/*
 * x	   = matrix (n x p) of explanatory variables
 * y	   = vector ( n )   of responses
 * weights = robustness weights wt[] * y[]	(of length n)
 * res	   = residuals	y[] - x[,] * beta	(of length n)
 * conv:  FALSE means do kk refining steps
 *	  TRUE  means refine until convergence(rel_tol, max_k)
 *      In the latter case, 'conv' *returns* TRUE if refinements converged
 * beta_cand= candidate beta[] (of length p)	Input *and* Output
 * is	    = initial scale			input

 * beta_ref = resulting beta[] (of length p)	Output
 * scale    = final scale			Output

 * tmp	    = aux vector of length n
 * tmp2	    = aux vector of length n
 * tmp_mat  = aux matrix p x p
 * tmp_mat2 = aux matrix p x (p+1)
*/

    int i,j, zeroes=0, one = 1;
    Rboolean converged = FALSE;/* Wall */
    double s0;

    for(j=0; j < n; j++) {
	res[j] = y[j] - F77_CALL(ddot)(&p, x[j], &one, beta_cand, &one);
	if( fabs(res[j]) < ZERO )
	    zeroes++;
    }
/* if "perfect fit", return it with a 0 assoc. scale */
    if( zeroes > (((double)n + (double)p)/2.) ) /* <<- FIXME: depends on 'b' ! */
    {
	COPY_beta(beta_cand, beta_ref);
	*scale = 0.;
	return 0;
    }

    if( initial_scale < 0. )
	initial_scale = MAD(res, n, 0., tmp, tmp2);
    s0 = initial_scale;
    if( *conv )
	kk = max_k;

    for(i=0; i < kk; i++) {

	/* one step for the scale */
      s0 = s0 * sqrt( sum_rho_sc(res, s0, n, p, rhoc, ipsi) / b );
	/* compute weights for IRWLS */
      get_weights_rhop(res, s0, n, rhoc, ipsi, weights);
	/* compute the matrix for IRWLS */
	r_sum_w_x_xprime(x, weights, n, p, tmp_mat, tmp_mat2);
	/* compute the vector for IRWLS */
	for(j=0; j < n; j++)
	    weights[j] *= y[j];
	r_sum_w_x(x, weights, n, p, tmp, tmp2);
	for(j=0; j < p; j++)
	    tmp_mat2[j][p] = tmp2[j];
	/* solve the system for IRWLS */
	lu(tmp_mat2, &p, beta_ref);
	if(*conv) { /* check for convergence */
	    double del = norm_diff(beta_cand, beta_ref, p);
	    double nrmB= norm(beta_cand, p);
	    if(trace_lev >= 3)
		Rprintf(" i = %d, ||b[i]||= %.12g, ||b[i] - b[i-1]|| = %.15g\n",
			i, nrmB, del);
	    converged = (del < rel_tol * fmax2(rel_tol, nrmB));
	    if(converged)
		break;
	}
	for(j=0; j < n; j++)
	    res[j] = y[j] - F77_CALL(ddot)(&p, x[j], &one, beta_ref, &one);
	COPY_beta(beta_ref, beta_cand);
    } /* for(i = 0; i < kk ) */

    if(*conv) {
	/* was "if(0)",	 since default lead to 'NOT converged' */
	if(!converged) {
	    *conv = FALSE;
	    warning("S refinements did not converge (to tol=%g) in %d iterations",
		    rel_tol, i);
	}
	if(trace_lev >= 2)
	    Rprintf("refinements %sconverged in %d iterations\n",
		    converged ? " " : "NOT ", i);
    }
    *scale = s0;
    return i; /* number of refinement steps */
} /* refine_fast_s() */


#ifdef _UNUSED_
void fast_s_irwls(double **x, double *y,
		double *weights, int n, int p, double *beta_ref,
		double **tmp_mat, double *tmp, double *tmp2)
{
    int i;
    for(i=0; i < n; i++)
	tmp[i] = weights[i] * y[i];
    mat_prime_vec(x, tmp, tmp2, n, p);
    mat_prime_mat_w(x, weights, tmp_mat, n, p);
    for(i=0; i < p; i++)
	tmp_mat[i][p] = tmp2[i];
    lu(tmp_mat, &p, beta_ref);
}
#endif

void get_weights_rhop(double *r, double s, int n, double rhoc, int ipsi,
		      double *w)
{
    int i;
    /* double a; */
    for(i=0; i < n; i++) 
      w[i] = wgt(r[i] / s, rhoc, ipsi);
      /* { */
      /* 	a = r[i] / s / rhoc; */
      /* 	if( fabs(a) > 1 ) */
      /* 	  w[i] = 0; */
      /* 	else { */
      /* 	  a = (1 - a)*(1 + a); */
      /* 	  w[i] = a * a; */
      /* 	} */
      /* } */
}

double find_scale(double *r, double b, double rhoc, int ipsi,
		  double initial_scale, int n, int p)
{
    int max_it = MAX_ITER_FIND_SCALE, it = 0;
    double e = 1, scale;

    while( (++it < max_it) && (fabs(e) > ZERO) )
    {
	scale = initial_scale *
	  sqrt( sum_rho_sc(r, initial_scale, n, p, rhoc, ipsi) / b ) ;
	e = fabs( scale / initial_scale - 1);
	initial_scale = scale;
    }
    return(scale);
}

int find_max(double *a, int n)
{
    if(n==1)
	return(0);
    else {
	int i, k = 0;
	double tt = a[0];
	for(i=1; i < n; i++)
	    if(tt < a[i]) {
		tt = a[i];
		k = i;
	    }
	return(k);
    }
}

void r_sum_w_x(double **x, double *w, int n, int p,
	       double *tmp, double *sum)
{
/*
 * Given a matrix x (n x p) and a vector w of n weights,
 * computes the p-vector  sum[] := \sum_{i=1}^n w_i x_i[]  , i.e.,
 *                        sum_j := \sum_{i=1}^n x_{ij} w_i , i.e., sum = X %*% w
 * Need space for p doubles in *tmp
*/
    int i;
    zero_vec(sum, p);
    for(i=0; i < n; i++) {
	scalar_vec(x[i], w[i], tmp, p);
	sum_vec(sum, /* + */ tmp, /* -> */ sum, p);
    }
}


void r_sum_w_x_xprime(double **x, double *w, int n, int p,
		      double **tmp, double **ans)
{
/* Given a matrix x (n x p) and a vector w ( n ) , compute the matrix
 * ans := \sum_{i=1}^n  w_i x_i x_i'
 * Need space for p x p "doubles" in tmp
*/
    int i;
    zero_mat(ans, p, p);

    for(i=0; i < n; i++) {
	outer_vec_vec(tmp, x[i], x[i], p);
	scalar_mat(tmp, w[i], tmp, p, p);
	sum_mat(ans, tmp, ans, p, p);
    }
}


double sum_rho(double *x, int n, double c, int ipsi)
{
    int i;
    double s = 0;
    for(i=0; i < n; i++)
      s += rho(x[i], c, ipsi);
    return(s);
}

double sum_rho_sc(double *r, double scale, int n, int p, double c, int ipsi)
{
    int i;
    double s = 0;
    for(i=0; i < n; i++)
      s += rho(r[i]/scale, c, ipsi);
    return(s / ((double) n - p));
}

#define _USE_BLAS_

/* ||x||_2 */
double norm(double *x, int n)
{
#ifdef  _USE_BLAS_
    int one = 1;
    return(F77_CALL(dnrm2)(&n, x, &one));
#else
    int i;
    double s = 0.;
    for(i=0; i < n; i++)
	s += x[i] * x[i];
    return(sqrt(s));
#endif
}

double norm1(double *x, int n)
{
#ifdef  _USE_BLAS_
    int one = 1;
    return(F77_CALL(dasum)(&n, x, &one));
#else
    int i; double s = 0.;
    for(i=0; i < n; i++)
	s += fabs(x[i]);
    return(s);
#endif
}

/* ||x-y||_2 */
double norm_diff(double *x, double *y, int n)
{
    int i;
    double s = 0;
    for(i=0; i < n; i++)
	s += (x[i]-y[i])*(x[i]-y[i]);
    return( sqrt(s) );
}

/* ||x-y||_1 */
double norm1_diff(double *x, double *y, int n)
{
    int i;
    double s = 0;
    for(i=0; i < n; i++)
	s += fabs(x[i]-y[i]);
    return(s);
}


double MAD(double *a, int n, double center, double *b,
			double *tmp)
{
/* if center == 0 then do not center */
    int i;
/*     if( fabs(center) > 0.) { */
	for(i=0; i < n; i++)
	    b[i] = a[i] - center;
/*     } */
    return( median_abs(b,n,tmp) * 1.4826 );
}

double median(double *x, int n, double *aux)
{
    double t;
    int i;
    for(i=0; i < n; i++) aux[i]=x[i];
    if ( (n/2) == (double) n / 2 )
	t = ( kthplace(aux,n,n/2) + kthplace(aux,n,n/2+1) ) / 2.0 ;
    else t = kthplace(aux,n, n/2+1 ) ;
    return(t);
}

double median_abs(double *x, int n, double *aux)
{
    double t;
    int i;
    for(i=0; i < n; i++) aux[i]=fabs(x[i]);
    if ( (n/2) == (double) n / 2 )
	t = ( kthplace(aux,n,n/2) + kthplace(aux,n,n/2+1) ) / 2.0 ;
    else	t = kthplace(aux,n, n/2+1 ) ;
    return(t);
}



void mat_prime_mat_w(double **a, double *w, double **c,
		     int n, int m)
{
/* compute  C = A' W A (where W is a diagonal matrix given by vector w[]) */
    int i,j,k;
    for(i=0; i < m; i++) {
	for(j=0; j < m; j++) {
	    c[i][j] = 0;
	    for(k=0; k < n; k++)
		c[i][j] += a[k][i] * w[k] * a[k][j];
	}
    }
}

void mat_prime_vec(double **a, double *b, double *c, int n, int m)
{
    int i,j;
    for(i=0; i < m; i++)
	for(c[i]=0,j=0; j < n; j++) c[i] += a[j][i] * b[j];
}


void disp_vec(double *a, int n)
{
    int i;
    for(i=0; i < n; i++) Rprintf("%lf ",a[i]);
    Rprintf("\n");
}


void disp_mat(double **a, int n, int m)
{
    int i,j;
    for(i=0; i < n; i++) {
	Rprintf("\n");
	for(j=0; j < m; j++) Rprintf("%10.8f ",a[i][j]);
    }
    Rprintf("\n");
}

void R_find_D_scale(double *rr, double *kkappa, double *ttau, int *llength, 
		    double *sscale, double *cc, int *iipsi, int *ttype, double *rel_tol,
		    int *max_k, int *converged)
{
  /* compute D_scale using iterative algorithm
     type: 1: d1
           2: d2
	   3: dt1
	   4: dt2
  */
  double scale, tsum1, tsum2, w, a;

  *converged = 0;
  for (int k=0; k < *max_k; k++) {
    scale = *sscale;
    tsum1 = 0;
    tsum2 = 0;
    // calculate weights
    for(int i=0; i < *llength; i++) {
      w = wgt(rr[i] / ttau[i] / scale, *cc, *iipsi);
      switch(*ttype) {
      case 1: // d1
	     a = rr[i]/ttau[i];
	     tsum1 += a*a*w;
	     tsum2 += w; break;
      case 2: // d2
	     a = rr[i]/ttau[i]*w;
	     tsum1 += a*a;
	     tsum2 += w*w; break;
      default:
      case 3: // dt1
	     tsum1 += rr[i]*rr[i]*w;
	     tsum2 += w*ttau[i]*ttau[i]; break;
      case 4: // d2
	     a = rr[i]*w;
	     tsum1 += a*a;
	     a = ttau[i]*w;
	     tsum2 += a*a; break;
      };
    }

    *sscale = sqrt(tsum1 / tsum2 / *kkappa);

    // Rprintf("\n type = %d, scale = %10.8f \n", *ttype, *sscale);

    if (fabs(scale - *sscale) < *rel_tol * fmax2(*rel_tol, scale)) {
      *converged = 1;
      break;
    }
  }
}
