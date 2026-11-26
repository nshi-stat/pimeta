// Quantile table for sampling tau^2 from its confidence distribution
// (Nagashima et al. 2019).
//
// The confidence distribution of tau^2 is H(tau^2) = 1 - Pr(Q <= q_obs | tau^2).
// Instead of solving H(tau^2) = U for each of the B uniform draws U, the
// inverse function H^{-1} is evaluated exactly at M - 1 nodes and the draws are
// obtained by monotone piecewise cubic Hermite (Fritsch--Carlson / PCHIP)
// interpolation of tau^2 against p.
//
// Nodes (Chebyshev--Lobatto type, dense near both ends):
//   p_j = H(0) + (1 - H(0)) * (1 - cos(pi * j / M)) / 2,  j = 1, ..., M - 1.

#ifndef PIMETA_CDTABLE_H
#define PIMETA_CDTABLE_H

#include <vector>
#include <cmath>
#include <algorithm>

namespace pimetacd {

// nodes p_1, ..., p_{M-1}
inline void cd_nodes(const double h0, const int M, std::vector<double>& p) {
  const double pi = 3.141592653589793238462643383280;
  p.resize(M - 1);
  for (int j = 1; j < M; j++) {
    p[j - 1] = h0 + (1.0 - h0)*(1.0 - std::cos(pi*(double)j/(double)M))*0.5;
  }
}

// one-sided, shape-preserving three-point end slope (Moler, 2004)
inline double pchip_end(const double h0, const double h1,
                        const double del0, const double del1) {
  double d = ((2.0*h0 + h1)*del0 - h0*del1)/(h0 + h1);
  if (d*del0 <= 0.0) {
    d = 0.0;
  } else if (del0*del1 <= 0.0 && std::fabs(d) > std::fabs(3.0*del0)) {
    d = 3.0*del0;
  }
  return d;
}

// Fritsch--Carlson slopes; x must be strictly increasing, y non-decreasing
inline void pchip_slopes(const std::vector<double>& x, const std::vector<double>& y,
                         std::vector<double>& d) {
  const int n = x.size();
  d.assign(n, 0.0);
  if (n < 2) return;
  std::vector<double> h(n - 1), del(n - 1);
  for (int k = 0; k < n - 1; k++) {
    h[k] = x[k + 1] - x[k];
    del[k] = (y[k + 1] - y[k])/h[k];
  }
  if (n == 2) {
    d[0] = del[0];
    d[1] = del[0];
    return;
  }
  for (int k = 1; k < n - 1; k++) {
    if (del[k - 1]*del[k] > 0.0) {
      const double w1 = 2.0*h[k] + h[k - 1];
      const double w2 = h[k] + 2.0*h[k - 1];
      d[k] = (w1 + w2)/(w1/del[k - 1] + w2/del[k]);
    } else {
      d[k] = 0.0;
    }
  }
  d[0] = pchip_end(h[0], h[1], del[0], del[1]);
  d[n - 1] = pchip_end(h[n - 2], h[n - 3], del[n - 2], del[n - 3]);
}

// evaluate the interpolant at u (x[0] <= u <= x[n-1])
inline double pchip_eval(const std::vector<double>& x, const std::vector<double>& y,
                         const std::vector<double>& d, const double u) {
  const int n = x.size();
  int j = std::upper_bound(x.begin(), x.end(), u) - x.begin();
  if (j < 1) j = 1;
  if (j > n - 1) j = n - 1;
  const int i = j - 1;
  const double hh = x[j] - x[i];
  const double t = (u - x[i])/hh;
  const double t1 = 1.0 - t;
  const double h00 = (1.0 + 2.0*t)*t1*t1;
  const double h10 = t*t1*t1;
  const double h01 = t*t*(3.0 - 2.0*t);
  const double h11 = t*t*(t - 1.0);
  return h00*y[i] + h10*hh*d[i] + h01*y[j] + h11*hh*d[j];
}

} // namespace pimetacd

#endif
