// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;
#include <Rcpp.h>
#include <boost/math/distributions/non_central_t.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

// [[Rcpp::export]]
double rcpp_qt(double p, double nu, double delta){
  return boost::math::quantile(boost::math::non_central_t(nu, delta), p);
}

// [[Rcpp::export]]
double rcpp_qchisq(double p, double nu, double delta){
  return boost::math::quantile(boost::math::non_central_chi_squared(nu, delta), p);
}


double integrand(double u, double l, double p, double k, double d){
  double x = u/(1-u);
  double q = boost::math::quantile(boost::math::non_central_chi_squared(1, x*x), p);
  double pchisq = boost::math::cdf(boost::math::complement(boost::math::chi_squared(l), l*q/k/k));
  double w = 2*boost::math::pdf(boost::math::normal(0, d), x);
  return w * pchisq / ((1-u)*(1-u));
}

class Integrand: public Func
{
private:
  double l;
  double p;
  double k;
  double d;
public:
  Integrand(double l_, double p_, double k_, double d_) : l(l_), p(p_), k(k_), d(d_) {}
  
  double operator()(const double& u) const
  {
    return integrand(u, l, p, k, d);
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector integral(double l, double p, double k, double d){
  Integrand f(l, p, k, d);
  double err_est;
  int err_code;
  const double res = integrate(f, 0, 1, err_est, err_code);
  Rcpp::NumericVector out = Rcpp::NumericVector::create(res);
  out.attr("err_est") = err_est;
  out.attr("err_code") = err_code;
  return out;
}