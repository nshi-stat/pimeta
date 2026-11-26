#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include "pimeta.h"
#include "cdtable.h"
// [[Rcpp::depends("RcppEigen")]]

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins("openmp")]]
#endif

using namespace Rcpp;

static inline double tabval(const std::vector<double>& px, const std::vector<double>& ty,
                            const std::vector<double>& dy, const double cpar, const double u) {
  const double x = -std::log1p(-u);
  const double v = std::exp(pimetacd::pchip_eval(px, ty, dy, x)) - cpar;
  return v < 0.0 ? 0.0 : v;
}

// Root of H(tau2) = u on a bracket grown by x100 from a typical size.
// H is increasing in tau2, so the root is the same as that of zeroin on
// [lower, mupper]; growing the bracket only saves iterations.
static void findRootGrow(const double u, const double qa, const Eigen::VectorXd& sigma,
                         const Eigen::MatrixXd& A, const int k, const double mode,
                         const int maxit1, const double eps, const double lower,
                         const double mupper, const int maxit2, const double tol,
                         double* tau2, int* status) {
  double lo = lower, hi = std::max(2.0*lower, 10.0*(1.0 + qa/(double)(k - 1)));
  while (hi < mupper) {
    const double f = fx(hi, u, 0, qa, sigma, A, k, mode, maxit1, eps);
    if (!(f < 0.0)) break; // root below hi (or NaN: reported by findRootTau2)
    lo = hi;
    hi *= 100.0;
  }
  findRootTau2(u, 0, qa, sigma, A, k, mode, maxit1, eps,
               lo, std::min(hi, mupper), maxit2, tol, tau2, status);
}

// Random numbers from the confidence distribution of tau^2.
//
// table = 0: exact inversion for each draw (the method of pimeta <= 1.1.3).
// table = 1: quantile table of size `nodes` with monotone cubic Hermite
//            interpolation; draws outside [p_1, p_{M-1}] are inverted exactly
//            (with the same tolerance as the nodes).
//
// Standardised scale: the computation uses sigma / s and (y - weighted mean) / s,
// s^2 = mean(sigma^2), and the draws are multiplied by s^2. The confidence
// distribution is exactly scale-equivariant, so nothing changes mathematically,
// but `tol` and the search bound become scale-free. (In pimeta <= 1.1.3 with
// upper = 1000 on the original scale, data with tau^2 > 1000 gave tau^2 = Inf
// for every draw, and on a tiny scale tol = 1.2e-4 exceeded tau^2.) The centring
// removes the cancellation in Q for data with a large offset.
// Upper end of the search: max(upper / s^2, 1e8, 1e6 * Q) on the standardised scale.
//
// When all y are identical (Q = 0), H(0) = 1 and every draw is 0. (Before, Q could
// be -1e-16, the distribution function failed, and the re-draws never ended.)
//
// Unlike rtau2CppWrap(), the R random number generator is never called inside
// the OpenMP parallel region. All uniform random numbers are generated before
// the parallel region, and the (rare) re-draws needed when the distribution
// function of Q could not be evaluated are done serially afterwards. Hence the
// results do not depend on the number of threads.
//
// [[Rcpp::export]]
NumericVector rtau2CDCppWrap(const int n, const Eigen::VectorXd& y, const Eigen::VectorXd& sigma,
                             const double mode, const int maxit1, const double eps,
                             const double lower, const double upper, const int maxit2,
                             const double tol, const int nthread, const int table,
                             const int nodes) {
  
  NumericVector rnd(n);
  NumericVector pvec = runif(n);
  double* prnd = rnd.begin();
  double* pp = pvec.begin();
  std::vector<int> redraw(n, 0);
  
  // standardisation
  const int k = sigma.size();
  double s2 = 0.0, sw = 0.0, swy = 0.0;
  for (int i = 0; i < k; i++) s2 += sigma[i]*sigma[i];
  s2 /= (double)k;
  const double ss = std::sqrt(s2);
  Eigen::VectorXd sn(k), yn(k);
  for (int i = 0; i < k; i++) {
    sn[i] = sigma[i]/ss;
    const double w = 1.0/(sn[i]*sn[i]);
    sw += w;
    swy += w*y[i];
  }
  for (int i = 0; i < k; i++) yn[i] = (y[i] - swy/sw)/ss;
  
  Eigen::MatrixXd A = getA(sn);
  const double qa = std::max(0.0, getqa(yn, A));
  if (!(qa > 0.0)) {
    // all y identical: Pr(Q <= 0 | tau^2) = 0, i.e. H(0) = 1
    for (int i = 0; i < n; i++) prnd[i] = 0.0;
    return rnd;
  }
  const double lowern = lower/s2;
  const double mupper = std::max(std::max(upper/s2, 1e8), 1e6*qa);
  const double tolfine = std::min(tol, 1e-10);
  const double zeroval = fx(0.0, 0.0, 0, qa, sn, A, k, mode, maxit1, eps);
  
#ifdef _OPENMP
  omp_set_num_threads(std::max(1, nthread));
#endif
  
  // quantile table (standardised scale: c = mean(sn^2) = 1)
  std::vector<double> px, ty, dy;
  const double cpar = 1.0;
  bool usetable = (table != 0) && (lower <= 0.0) && (nodes >= 4) &&
    std::isfinite(zeroval) && (zeroval < 1.0);
  
  if (usetable) {
    std::vector<double> pn;
    pimetacd::cd_nodes(zeroval, nodes, pn);
    const int nn = pn.size();
    std::vector<double> tn(nn, NA_REAL);
    std::vector<int> sn2(nn, 0);
#pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < nn; j++) {
      double t;
      int st;
      findRootGrow(pn[j], qa, sn, A, k, mode, maxit1, eps,
                   0.0, mupper, maxit2, tolfine, &t, &st);
      tn[j] = t;
      sn2[j] = st;
    }
    // keep the longest valid prefix (nodes beyond H(mupper) are dropped;
    // draws there are handled by exact inversion and become Inf)
    // interpolation is done on the transformed scale
    //   x = -log(1 - p),  y = log(tau^2 + c),
    // on which the quantile function is nearly linear in the upper tail
    // (tau^2 + c ~ q_obs c / chisq_{K-1, 1-p} for equal within-study variances).
    for (int j = 0; j < nn; j++) {
      if (sn2[j] != 0 || !std::isfinite(tn[j])) break;
      const double xj = -std::log1p(-pn[j]);
      if (!std::isfinite(xj) || (!px.empty() && !(xj > px.back()))) break;
      double t = tn[j];
      if (!ty.empty() && std::log(t + cpar) < ty.back()) t = std::exp(ty.back()) - cpar;
      px.push_back(xj);
      ty.push_back(std::log(t + cpar));
    }
    if (px.size() < 4) {
      usetable = false;
    } else {
      pimetacd::pchip_slopes(px, ty, dy);
    }
  }
  
  const double plo = usetable ? -std::expm1(-px.front()) : 0.0;
  const double phi = usetable ? -std::expm1(-px.back()) : -1.0;
  // exact inversions: tol for sampling = "exact" (as before), the node tolerance
  // for the tails of sampling = "table"
  const double tolx = usetable ? tolfine : tol;
  
#pragma omp parallel for schedule(dynamic, 256)
  for (int i = 0; i < n; i++) {
    const double u = pp[i];
    if (zeroval >= u) {
      prnd[i] = 0.0;
    } else if (usetable && u >= plo && u <= phi) {
      prnd[i] = s2*tabval(px, ty, dy, cpar, u);
    } else {
      double t;
      int st;
      findRootGrow(u, qa, sn, A, k, mode, maxit1, eps,
                   lowern, mupper, maxit2, tolx, &t, &st);
      if (st == 1) {
        prnd[i] = R_PosInf;
      } else if (st == 2) {
        prnd[i] = NA_REAL;
        redraw[i] = 1;
      } else {
        prnd[i] = s2*t;
      }
    }
  }
  
  // serial re-draws (R RNG is not thread safe)
  int nfail = 0;
  for (int i = 0; i < n; i++) {
    int iter = 0;
    while (redraw[i] == 1) {
      if (++iter > 100) {
        nfail++;
        break;
      }
      const double u = R::runif(0.0, 1.0);
      if (zeroval >= u) {
        prnd[i] = 0.0;
        redraw[i] = 0;
      } else if (usetable && u >= plo && u <= phi) {
        prnd[i] = s2*tabval(px, ty, dy, cpar, u);
        redraw[i] = 0;
      } else {
        double t;
        int st;
        findRootGrow(u, qa, sn, A, k, mode, maxit1, eps,
                     lowern, mupper, maxit2, tolx, &t, &st);
        if (st == 1) {
          prnd[i] = R_PosInf;
          redraw[i] = 0;
        } else if (st != 2) {
          prnd[i] = s2*t;
          redraw[i] = 0;
        }
      }
    }
  }
  
  int nna = 0;
  for (int i = 0; i < n; i++) {
    if (std::isnan(prnd[i])) nna++;
  }
  if (nfail > 0 || nna > 0) {
    Rcpp::warning("%d random numbers of tau^2 could not be generated (NA); consider increasing 'maxit1' or 'maxit2'.", nna);
  }
  
  return rnd;
  
}
