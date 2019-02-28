#include <Rcpp.h>
#include <RcppEigen.h>
#include "pimeta.h"
// [[Rcpp::depends("RcppEigen")]]

using namespace Rcpp;

// [[Rcpp::export]]
List dwchisqCpp(const double q, const Eigen::VectorXd& lambda, const Eigen::VectorXi& mult,
                const Eigen::VectorXd& delta, const int n, const double mode,
                const int maxit, const double eps) {
  
  double prob;
  int ifault;

  pQCpp2(q, lambda, mult, delta, n, mode, maxit, eps, &prob, &ifault);
  
  return List::create(
    Named("prob") = prob,
    Named("ifault") = ifault
  );
  
}

// [[Rcpp::export]]
List bootPICppWrap(const Eigen::VectorXd& rnd, const Eigen::VectorXd& y,
                   const Eigen::VectorXd& sigma, const double alpha) {
  
  double lbpi, ubpi, lbci, ubci;
  
  bootPICpp(rnd, y, sigma, alpha, &lbpi, &ubpi, &lbci, &ubci);
  
  return List::create(
    Named("lpi") = lbpi,
    Named("upi") = ubpi,
    Named("lci") = lbci,
    Named("uci") = ubci
  );
  
}

void bootPICpp(const Eigen::VectorXd& rnd, const Eigen::VectorXd& y,
               const Eigen::VectorXd& sigma, const double alpha,
               double* lbpi, double* ubpi, double* lbci, double* ubci) {
  
  double index, h;
  int lo, hi, n = rnd.size(), k = sigma.size();
  
  /* temporary  */
  double s1, s2, tauh2;
  Eigen::VectorXd w;
  
  w = sigma.array().pow(-2.0);
  s1 = w.sum();
  s2 = sigma.array().pow(-4.0).sum();
  tauh2 = std::max(0.0, ((w.array()*(y.array() - (w.array()*y.array()).sum()/s1).pow(2)).sum() - (k - 1))/(s1 - s2/s1));
  w = (sigma.array().pow(2) + tauh2).inverse();
  
  Eigen::VectorXd wBrs, dpt, dpn, on = Eigen::VectorXd::Ones(n), ok = Eigen::VectorXd::Ones(k);
  Eigen::MatrixXd wB, mB, vB, yB = on*y.transpose();
  
  wBrs = sigma.array().pow(2);  // The variable wBrs is temporarily used. We need VectorXd sigma^2.
  wB = (on*wBrs.transpose() + rnd*ok.transpose()).array().inverse();
  wBrs = wB.rowwise().sum();
  mB = ((yB.array()*wB.array()).rowwise().sum()).array()/wBrs.array();
  vB = ((wB.array()*((yB - mB*ok.transpose()).array().pow(2))).rowwise().sum())/(k - 1)/wBrs.array();
  
  Eigen::Map<Eigen::VectorXd> rndt = as<Eigen::Map<Eigen::VectorXd> >(Rcpp::rt(n, k - 1));
  Eigen::Map<Eigen::VectorXd> rndn = as<Eigen::Map<Eigen::VectorXd> >(Rcpp::rnorm(n, 0.0, 1.0));
  dpt = mB.array() + rndn.array()*rnd.array().sqrt() - rndt.array()*vB.array().sqrt();
  dpn = mB.array() - rndt.array()*vB.array().sqrt();
  
  // lower quantile
  index = (n - 1)*alpha*0.5;
  lo = floor(index);
  hi = ceil(index);
  std::nth_element(dpt.data(), dpt.data() + lo, dpt.data() + n);
  std::nth_element(dpn.data(), dpn.data() + lo, dpn.data() + n);
  *lbpi = dpt[lo];
  *lbci = dpn[lo];
  if (index > lo) {
    h = index - lo;
    std::nth_element(dpt.data(), dpt.data() + hi, dpt.data() + n);
    std::nth_element(dpn.data(), dpn.data() + hi, dpn.data() + n);
    *lbpi = (1.0 - h)*(*lbpi) + h*dpt[hi];
    *lbci = (1.0 - h)*(*lbci) + h*dpn[hi];
  }
  
  // upper quantile
  index = (n - 1)*(1.0 - alpha*0.5);
  lo = floor(index);
  hi = ceil(index);
  std::nth_element(dpt.data(), dpt.data() + lo, dpt.data() + n);
  std::nth_element(dpn.data(), dpn.data() + lo, dpn.data() + n);
  *ubpi = dpt[lo];
  *ubci = dpn[lo];
  if (index > lo) {
    h = index - lo;
    std::nth_element(dpt.data(), dpt.data() + hi, dpt.data() + n);
    std::nth_element(dpn.data(), dpn.data() + hi, dpn.data() + n);
    *ubpi = (1.0 - h)*(*ubpi) + h*dpt[hi];
    *ubci = (1.0 - h)*(*ubci) + h*dpn[hi];
  }
  
}

//void exactCItau2(const Eigen::VectorXd& y, const Eigen::VectorXd& sigma,
//                  const double mode, const int maxit1, const double eps,
//                  const double lower, const double upper, const int maxit2,
//                  const double tol, const double alpha) {
//   
//   int i = 0, status;
//   double tau2;
//   
//   Eigen::MatrixXd A = getA(sigma);
//   int k = sigma.size();
//   double qa = getqa(y, A);
//   double mupper = std::max(upper, qa);
//   double zeroval = fx(0.0, 0.0, 0, qa, sigma, A, k, mode, maxit1, eps);
//   
//   if (zeroval >= pvec(i)) {
//     rnd(i) = 0.0;
//   } else {
//     findRootTau2(pvec(i), 0, qa, sigma, A, k, mode, maxit1, eps, lower, mupper, maxit2, tol, &tau2, &status);
//     if (status != 2) {
//       if (status == 1) {
//         rnd(i) = R_PosInf;
//       } else {
//         rnd(i) = tau2;
//       }
//     } else {
//       rnd = R_NaN;
//     }
//   }
//   
// }

// [[Rcpp::export]]
NumericVector rtau2CppWrap(const int n, const Eigen::VectorXd& y, const Eigen::VectorXd& sigma,
                           const double mode, const int maxit1, const double eps,
                           const double lower, const double upper, const int maxit2, const double tol) {
  
  int i = 0, status;
  double tau2;
  NumericVector rnd(n), pvec(n);
  pvec = runif(n);
  
  Eigen::MatrixXd A = getA(sigma);
  int k = sigma.size();
  double qa = getqa(y, A);
  double mupper = std::max(upper, qa);
  double zeroval = fx(0.0, 0.0, 0, qa, sigma, A, k, mode, maxit1, eps);
  
  while (i < n) {
    if (zeroval >= pvec(i)) {
      rnd(i) = 0.0;
      i++;
    } else {
      findRootTau2(pvec(i), 0, qa, sigma, A, k, mode, maxit1, eps, lower, mupper, maxit2, tol, &tau2, &status);
      if (status != 2) {
        if (status == 1) {
          rnd(i) = R_PosInf;
        } else {
          rnd(i) = tau2;
        }
        i++;
      } else {
        pvec(i) = R::runif(0.0, 1.0);
      }
    }
  }
  
  return rnd;
  
}

double getqa(const Eigen::VectorXd& y, const Eigen::MatrixXd& A) {
  
  return y.transpose()*A*y.cast<double>();
  
}

Eigen::MatrixXd getA(const Eigen::VectorXd& sigma) {
  
  Eigen::VectorXd w;
  Eigen::MatrixXd A;
  
  w = sigma.array().pow(-2.0);
  A = w.asDiagonal();
  A = A - (w * w.transpose())/w.sum();
  
  return A;
  
}

void findRootTau2(const double alpha, const int side, const double qa, const Eigen::VectorXd& sigma,
                  const Eigen::MatrixXd& A, const int k, const double mode, const int maxit1, const double eps,
                  const double lower, const double upper, const int maxit2, const double tol,
                  double* tau2, int* status) {
  
  // A tlanslation of the "zeroin" by Forsythe GE, Malcolm MA, Moler CB.
  // Computer Methods for Mathematical Computations. Prentice-Hall, 1977.
  
  double a, b, fa, fb, tol1 = NA_REAL;
  a = lower;
  b = upper;
  fa = fx(a, alpha, side, qa, sigma, A, k, mode, maxit1, eps);
  fb = fx(b, alpha, side, qa, sigma, A, k, mode, maxit1, eps);
  *tau2 = NA_REAL;
  *status = 0;
  
  if (fa*fb >= 0.0) {
    // equal sign
    *status = 1;
  } else {
    
    int iter = 0;
    double c, fc, d, e, xm, p, q, r, s;
    double meps = std::numeric_limits<double>::epsilon();
    
    c = a;
    fc = fa;
    d = b - a;
    e = d;
    
    while (1) {
      
      iter += 1;
      
      if (std::isnan(fa) || std::isnan(fb)) {
        *status = 2;
        break;
      }
      
      if (fabs(fc) < fabs(fb)) {
        a = b;
        b = c;
        c = a;
        fa = fb;
        fb = fc;
        fc = fa;
      }
      
      tol1 = 2.0*meps*fabs(b) + 0.5*tol;
      xm = 0.5*(c - b);
      
      if (fabs(xm) <= tol1 || fb == (double) 0) {
        *tau2 = b;
        break;
      }
      
      if (fabs(e) < tol1 && fabs(fa) <= fabs(fb)) {
        // bisection
        d = xm;
        e = d;
      } else {
        if (a != c) {
          // quadratic interpolation
          q = fa/fc;
          r = fb/fc;
          s = fb/fa;
          p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
          q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        } else {
          // linear interpolation
          s = fb/fa;
          p = 2.0*xm*s;
          q = 1.0 - s;
        }
        if (p > 0.0) {
          q = -q;
        }
        p = fabs(p);
        if ((2.0*p) < (3.0*xm*q - fabs(tol1*q)) && p < fabs(0.5*e*q)) {
          e = d;
          d = p/q;
        } else {
          d = xm;
          e = d;
        }
      }
      
      a = b;
      fa = fb;
      if (fabs(d) > tol1) {
        b = b + d;
      } else {
        if (xm < 0.0) {
          b = b - tol1;
        } else {
          b = b + tol1;
        }
      }
      fb = fx(b, alpha, side, qa, sigma, A, k, mode, maxit1, eps);
      if (fb*(fc/fabs(fc)) > 0.0) {
        c = a;
        fc = fa;
        d = b - a;
        e = d;
      }
      
      if (iter >= maxit2) {
        *status = 3;
        break;
      }
    }
  }
  
}

double fx(const double tau2, const double alpha, const int side, const double qa,
          const Eigen::VectorXd& sigma, const Eigen::MatrixXd& A,
          const int k, const double mode, const int maxit, const double eps) {
  
  int ifault;
  double prob;
  Eigen::VectorXd tmp;
  Eigen::MatrixXd Z;
  
  tmp = (sigma.array().pow(2) + tau2).sqrt();
  Z = tmp.asDiagonal();
  Z = Z*A*Z;
  
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ZAZ(Z, false);
  Eigen::VectorXd lambda = ZAZ.eigenvalues().tail(k - 1);
  
  pQCpp(qa, lambda, k - 1, mode, maxit, eps, &prob, &ifault);
  
  if (ifault < 0 || ifault == 2 || ifault == 3) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  
  if (side == 0) {
    // lower
    return 1.0 - prob - alpha;
  } else {
    // upper
    return prob - alpha;
  }
  
}


void pQCpp(const double qa, const Eigen::VectorXd& lambda, const int n, const double mode,
           const int maxit, const double eps, double *prob, int *ifault) {
  
  // Ruben (1962), Farebrother (1984, AS 204)
  // Optimized for cochran's Q in meta-analysis
  
  if (n < 1 || qa < 0.0 || maxit < 1 || eps < 0.0) {
    *prob = -2.0;
    *ifault = 2;
  } else {
    int i, k, m, flg = 0;
    double ao, aoinv, z, beta, eps2, sum, sum1, hold, dans, lans, pans, tol = -200.0;
    double *gamma, *theta, *a, *b;
    gamma = new double[n];
    theta = new double[n];
    a = new double[maxit];
    b = new double[maxit];
    
    sum = lambda[0];
    beta = sum;
    
    for (i = 0; i < n; i++) {
      hold = lambda[i];
      if (hold <= 0.0) {
        *prob = -7.0;
        *ifault = -i;
        flg = 1;
        break;
      }
      if (beta > hold) {
        beta = hold;
      }
      if (sum < hold) {
        sum = hold;
      }
    }
    
    if (flg == 0) {
      
      if (mode > 0.0) {
        beta = mode*beta;
      } else {
        beta = 2.0/(1.0/beta + 1.0/sum);
      }
      
      k = 0;
      sum = 1.0;
      sum1 = 0.0;
      for (i = 0; i < n; i++) {
        hold = beta/lambda[i];
        gamma[i] = 1.0 - hold;
        sum = sum*hold;
        k = k + 1;
        theta[i] = 1.0;
      }
      
      ao = exp(0.5*(log(sum) - sum1));
      if (ao <= 0.0) {
        *prob = 0.0;
        *ifault = 1;
      } else {
        z = qa/beta;
        if ((k%2) == 0) {
          i = 2;
          lans = -0.5*z;
          dans = exp(lans);
          pans = 1.0 - dans;
        } else {
          i = 1;
          lans = -0.5*(z + log(z)) - 0.22579135264473;
          dans = exp(lans);
          pans = R::pnorm(sqrt(z), 0.0, 1.0, 1, 0) - R::pnorm(-sqrt(z), 0.0, 1.0, 1, 0);
        }
        k = k - 2;
        for (int j = i; j <= k; j = j + 2) {
          if (lans < tol) {
            lans = lans + log(z/(double)j);
            dans = exp(lans);
          } else {
            dans = dans*z/(double)j;
          }
          pans = pans - dans;
        }
        *prob = pans;
        eps2 = eps / ao;
        aoinv = 1.0 / ao;
        sum = aoinv - 1.0;
        
        for (m = 1; m < maxit + 1; m++) {
          sum1 = 0.0;
          
          for (i = 0; i < n; i++) {
            hold = theta[i]*gamma[i];
            theta[i] = hold;
            sum1 = sum1 + hold;
          }
          
          sum1 = 0.5*sum1;
          b[m - 1] = sum1;
          for (i = m - 1; i > 0; i--) {
            sum1 = sum1 + b[i - 1]*a[m - i - 1];
          }
          sum1 = sum1/(double)m;
          a[m - 1] = sum1;
          
          k = k + 2;
          if (lans < tol) {
            lans = lans + log(z/(double)k);
            dans = exp(lans);
          } else {
            dans = dans*z/(double)k;
          }
          
          pans = pans - dans;
          sum = sum - sum1;
          sum1 = pans*sum1;
          *prob = *prob + sum1;
          
          if (*prob < (-aoinv)) {
            *prob = -3.0;
            *ifault = 3;
            break;
          }
          if (fabs(pans*sum) < eps2) {
            if (fabs(sum1) < eps2) {
              *ifault = 0;
              flg = 1;
              break;
            }
          }
        }
        if (flg == 0) {
          *ifault = 4;
        }
        *prob = ao*(*prob);
        if (*prob < 0.0 || *prob > 1.0) {
          *ifault = *ifault + 5;
        } else if (dans < 0.0) {
          *ifault = *ifault + 6;
        }
      }
    }
    
    delete[] gamma;
    delete[] theta;
    delete[] a;
    delete[] b;
  }
  
}

void pQCpp2(const double qa, const Eigen::VectorXd& lambda, const Eigen::VectorXi& mult,
            const Eigen::VectorXd& delta, const int n, const double mode, const int maxit,
            const double eps, double *prob, int *ifault) {
  
  // Ruben (1962), Farebrother (1984, AS 204)

  if (n < 1 || qa < 0.0 || maxit < 1 || eps < 0.0) {
    *prob = -2.0;
    *ifault = 2;
  } else {
    int i, k, m, flg = 0;
    double ao, aoinv, z, beta, eps2, sum, sum1, hold, hold2, dans,
      lans, pans, tol = -200.0;
    double *gamma, *theta, *a, *b;
    gamma = new double[n];
    theta = new double[n];
    a = new double[maxit];
    b = new double[maxit];
    
    sum = lambda[0];
    beta = sum;
    
    for (i = 0; i < n; i++) {
      hold = lambda[i];
      if (hold <= 0.0 || mult[i] < 1 || delta[i] < 0.0) {
        *prob = -7.0;
        *ifault = -i;
        flg = 1;
        break;
      }
      if (beta > hold) {
        beta = hold;
      }
      if (sum < hold) {
        sum = hold;
      }
    }
    
    if (flg == 0) {
      
      if (mode > 0.0) {
        beta = mode*beta;
      } else {
        beta = 2.0/(1.0/beta + 1.0/sum);
      }
      
      k = 0;
      sum = 1.0;
      sum1 = 0.0;
      for (i = 0; i < n; i++) {
        hold = beta/lambda[i];
        gamma[i] = 1.0 - hold;
        sum = sum*pow(hold, mult[i]);
        sum1 = sum1 + delta[i];
        k = k + mult[i];
        theta[i] = 1.0;
      }
      
      ao = exp(0.5*(log(sum) - sum1));
      if (ao <= 0.0) {
        *prob = 0.0;
        *ifault = 1;
      } else {
        z = qa/beta;
        if ((k%2) == 0) {
          i = 2;
          lans = -0.5*z;
          dans = exp(lans);
          pans = 1.0 - dans;
        } else {
          i = 1;
          lans = -0.5*(z + log(z)) - 0.22579135264473;
          dans = exp(lans);
          pans = R::pnorm(sqrt(z), 0.0, 1.0, 1, 0) - R::pnorm(-sqrt(z), 0.0, 1.0, 1, 0);
        }
        k = k - 2;
        for (int j = i; j <= k; j = j + 2) {
          if (lans < tol) {
            lans = lans + log(z/(double)j);
            dans = exp(lans);
          } else {
            dans = dans*z/(double)j;
          }
          pans = pans - dans;
        }
        *prob = pans;
        eps2 = eps / ao;
        aoinv = 1.0 / ao;
        sum = aoinv - 1.0;
        
        for (m = 1; m < maxit + 1; m++) {
          sum1 = 0.0;
          
          for (i = 0; i < n; i++) {
            hold = theta[i];
            hold2 = hold*gamma[i];
            theta[i] = hold2;
            sum1 = sum1 + hold2*((double)mult[i]) + ((double)m)*delta[i]*(hold - hold2);
          }
          
          sum1 = 0.5*sum1;
          b[m - 1] = sum1;
          for (i = m - 1; i > 0; i--) {
            sum1 = sum1 + b[i - 1]*a[m - i - 1];
          }
          sum1 = sum1/(double)m;
          a[m - 1] = sum1;
          
          k = k + 2;
          if (lans < tol) {
            lans = lans + log(z/(double)k);
            dans = exp(lans);
          } else {
            dans = dans*z/(double)k;
          }
          
          pans = pans - dans;
          sum = sum - sum1;
          sum1 = pans*sum1;
          *prob = *prob + sum1;
          
          if (*prob < (-aoinv)) {
            *prob = -3.0;
            *ifault = 3;
            break;
          }
          if (fabs(pans*sum) < eps2) {
            if (fabs(sum1) < eps2) {
              *ifault = 0;
              flg = 1;
              break;
            }
          }
        }
        if (flg == 0) {
          *ifault = 4;
        }
        *prob = ao*(*prob);
        if (*prob < 0.0 || *prob > 1.0) {
          *ifault = *ifault + 5;
        } else if (dans < 0.0) {
          *ifault = *ifault + 6;
        }
      }
    }
    
    delete[] gamma;
    delete[] theta;
    delete[] a;
    delete[] b;
  }
  
}
