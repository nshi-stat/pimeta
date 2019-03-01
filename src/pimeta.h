#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

extern List pwchisqCpp(const double, const Eigen::VectorXd&, const Eigen::VectorXi&,
                       const Eigen::VectorXd&, const int, const double, const int, const double);

extern List bootPICppWrap(const Eigen::VectorXd&, const Eigen::VectorXd&,
                          const Eigen::VectorXd&, const double);

extern void bootPICpp(const Eigen::VectorXd&, const Eigen::VectorXd&,
                      const Eigen::VectorXd&, const double,
                      double*, double*, double*, double*);

extern NumericVector rtau2CppWrap(const int, const Eigen::VectorXd&, const Eigen::VectorXd&,
                                  const double, const int, const double,
                                  const double, const double, const int, const double);

extern double getqa(const Eigen::VectorXd&, const Eigen::MatrixXd&);

extern Eigen::MatrixXd getA(const Eigen::VectorXd&);

extern void findRootTau2(const double, const int, const double,
                         const Eigen::VectorXd&, const Eigen::MatrixXd&,
                         const int, const double, const int, const double,
                         const double, const double, const int, const double,
                         double*, int*);

extern double fx(const double, const double, const int, const double,
                 const Eigen::VectorXd&, const Eigen::MatrixXd&,
                 const int, const double, const int, const double);

extern void pQCpp(const double, const Eigen::VectorXd&, const int,
                  const double, const int, const double, double*, int*);

extern void pQCpp2(const double, const Eigen::VectorXd&, const Eigen::VectorXi&,
                   const Eigen::VectorXd&, const int, const double, const int,
                   const double, double*, int*);
