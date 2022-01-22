# include <RcppArmadillo.h>
# include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
int sgn(const double x) {
  return (x > 0) - (x < 0);
}

// [[Rcpp::export]]
double mad(const arma::vec& x) {
  return 1.482602 * arma::median(arma::abs(x - arma::median(x)));
}

// [[Rcpp::export]]
void updateL2(const arma::mat& Z, const arma::vec& res, arma::vec& grad, const double n1) {
  grad = -n1 * Z.t() * res;
}

// [[Rcpp::export]]
void updateExpectile(const arma::mat& Z, const arma::vec& res, const double tau, arma::vec& der, arma::vec& grad, const int n, const double rob, const double n1) {
  for (int i = 0; i < n; i++) {
    double cur = res(i);
    if (cur > rob) {
      der(i) = -2 * tau * rob;
    } else if (cur > 0) {
      der(i) = -2 * tau * cur;
    } else if (cur > -rob) {
      der(i) = 2 * (tau - 1) * cur;
    } else {
      der(i) = 2 * (1 - tau) * rob;
    }
  }
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
void updateHuber(const arma::mat& Z, const arma::vec& res, arma::vec& der, arma::vec& grad, const int n, const double rob, const double n1) {
  for (int i = 0; i < n; i++) {
    double cur = res(i);
    if (cur > rob) {
      der(i) = -rob;
    } else if (cur > -rob) {
      der(i) = -cur;
    } else {
      der(i) = rob;
    }
  }
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
arma::vec l2Reg(const arma::mat& Z, const arma::vec& Y, arma::vec& gradOld, arma::vec& gradNew, const double n1, const double tol = 0.0001, 
                const int iteMax = 5000) {
  updateL2(Z, Y, gradOld, n1);
  arma::vec beta = -gradOld, betaDiff = -gradOld;
  arma::vec res = Y - Z * beta;
  updateL2(Z, res, gradNew, n1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    updateL2(Z, res, gradNew, n1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  return beta;
}

// [[Rcpp::export]]
arma::vec expectile(const arma::mat& Z, const arma::vec& Y, const double tau, arma::vec& der, arma::vec& gradOld, arma::vec& gradNew, const int n, 
                    const double n1, const double tol = 0.0001, const double constTau = 1.345, const int iteMax = 5000) {
  double rob = constTau * mad(Y);
  updateExpectile(Z, Y, tau, der, gradOld, n, rob, n1);
  arma::vec beta = -gradOld, betaDiff = -gradOld;
  arma::vec res = Y - Z * beta;
  rob = constTau * mad(res);
  updateExpectile(Z, res, tau, der, gradNew, n, rob, n1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    rob = constTau * mad(res);
    updateExpectile(Z, res, tau, der, gradNew, n, rob, n1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  return beta;
}

// [[Rcpp::export]]
arma::vec huberReg(const arma::mat& Z, const arma::vec& Y, arma::vec& der, arma::vec& gradOld, arma::vec& gradNew, const int n, 
                   const double n1, const double tol = 0.0001, const double constTau = 1.345, const int iteMax = 5000) {
  double rob = constTau * mad(Y);
  updateHuber(Z, Y, der, gradOld, n, rob, n1);
  arma::vec beta = -gradOld, betaDiff = -gradOld;
  arma::vec res = Y - Z * beta;
  rob = constTau * mad(res);
  updateHuber(Z, res, der, gradNew, n, rob, n1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    rob = constTau * mad(res);
    updateHuber(Z, res, der, gradNew, n, rob, n1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  return beta;
}

// [[Rcpp::export]]
arma::mat standardize(arma::mat X, const arma::rowvec& mx, const arma::vec& sx1, const int p) {
  for (int i = 0; i < p; i++) {
    X.col(i) = (X.col(i) - mx(i)) * sx1(i);
  }
  return X;
}

// [[Rcpp::export]]
void updateGauss(const arma::mat& Z, const arma::vec& res, arma::vec& der, arma::vec& grad, const double tau, const double n1, const double h1) {
  der = arma::normcdf(-res * h1) - tau;
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
Rcpp::List twoStep(const arma::mat& X, arma::vec Y, const double alpha = 0.1, double h = 0.0, const double constTau = 1.345, 
                   const double tol = 0.0001, const int iteMax = 5000) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  if (h <= 0.05) {
    h = std::max(std::pow((std::log(n) + p) / n, 0.4), 0.05);
  }
  const double n1 = 1.0 / n;
  const double h1 = 1.0 / h;
  // standardiza design for more efficient computation
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec der(n);
  arma::vec gradOld(p + 1), gradNew(p + 1);
  // initialization via expectile
  arma::vec beta = expectile(Z, Y, alpha, der, gradOld, gradNew, n, n1, tol, constTau, iteMax);
  arma::vec quant = {alpha};
  beta(0) = arma::as_scalar(arma::quantile(Y - Z.cols(1, p) * beta.rows(1, p), quant));
  // first step: a smoothed quantile estimator
  arma::vec res = Y - Z * beta;
  updateGauss(Z, res, der, gradOld, alpha, n1, h1);
  beta -= gradOld;
  arma::vec betaDiff = -gradOld;
  res -= Z * betaDiff;
  updateGauss(Z, res, der, gradNew, alpha, n1, h1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double step = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      step = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -step * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    updateGauss(Z, res, der, gradNew, alpha, n1, h1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  // second step: an estimator for expected shortfall
  arma::vec w = Y - Z * beta;
  w = arma::min(w, arma::zeros(n));
  double mw = arma::mean(w);
  w -= mw;
  arma::vec theta = l2Reg(Z, w, gradOld, gradNew, n1, tol, iteMax);
  // transform back to the original scale
  beta.rows(1, p) %= sx1;
  beta(0) += my - arma::as_scalar(mx * beta.rows(1, p));
  theta.rows(1, p) %= sx1;
  theta(0) += mw - arma::as_scalar(mx * theta.rows(1, p));
  theta = theta / alpha + beta;
  return Rcpp::List::create(Rcpp::Named("beta") = beta, Rcpp::Named("theta") = theta);
}

// [[Rcpp::export]]
Rcpp::List twoStepNonstd(const arma::mat& X, arma::vec Y, const double alpha = 0.2, double h = 0.0, const double constTau = 1.345, 
                         const double tol = 0.0001, const int iteMax = 5000) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  if (h <= 0.05) {
    h = std::max(std::pow((std::log(n) + p) / n, 0.4), 0.05);
  }
  const double n1 = 1.0 / n;
  const double h1 = 1.0 / h;
  arma::mat Z = arma::join_rows(arma::ones(n), X);
  arma::vec der(n);
  arma::vec gradOld(p + 1), gradNew(p + 1);
  // initialization via expectile
  arma::vec beta = expectile(Z, Y, alpha, der, gradOld, gradNew, n, n1, tol, constTau, iteMax);
  arma::vec quant = {alpha};
  beta(0) = arma::as_scalar(arma::quantile(Y - Z.cols(1, p) * beta.rows(1, p), quant));
  // first step: a smoothed quantile estimator
  arma::vec res = Y - Z * beta;
  updateGauss(Z, res, der, gradOld, alpha, n1, h1);
  beta -= gradOld;
  arma::vec betaDiff = -gradOld;
  res -= Z * betaDiff;
  updateGauss(Z, res, der, gradNew, alpha, n1, h1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double step = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      step = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -step * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    updateGauss(Z, res, der, gradNew, alpha, n1, h1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  // second step: an estimator for expected shortfall
  arma::vec w = Y - Z * beta;
  w = arma::min(w, arma::zeros(n));
  arma::vec theta = arma::solve(Z, w);
  theta = theta / alpha + beta;
  return Rcpp::List::create(Rcpp::Named("beta") = beta, Rcpp::Named("theta") = theta);
}

// [[Rcpp::export]]
Rcpp::List oracle(const arma::mat& X, arma::vec Y, const arma::vec& beta, const double alpha = 0.1, const double constTau = 1.345, 
                  const double tol = 0.0001, const int iteMax = 5000) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const double n1 = 1.0 / n;
  arma::mat Z = arma::join_rows(arma::ones(n), X);
  arma::vec w = Y - Z * beta;
  w = arma::min(w, arma::zeros(n));
  // standardiza design for more efficient computation
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double mw = arma::mean(w);
  w -= mw;
  arma::vec der(n);
  arma::vec gradOld(p + 1), gradNew(p + 1);
  arma::vec theta = l2Reg(Z, w, gradOld, gradNew, n1, tol, iteMax);
  // transform back to the original scale
  theta.rows(1, p) %= sx1;
  theta(0) += mw - arma::as_scalar(mx * theta.rows(1, p));
  theta = theta / alpha + beta;
  return Rcpp::List::create(Rcpp::Named("beta") = beta, Rcpp::Named("theta") = theta);
}

// [[Rcpp::export]]
Rcpp::List twoStepRob(const arma::mat& X, arma::vec Y, const double alpha = 0.2, double h = 0.0, const double constTau = 1.345, 
                      const double tol = 0.0001, const int iteMax = 5000) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  if (h <= 0.05) {
    h = std::max(std::pow((std::log(n) + p) / n, 0.4), 0.05);
  }
  const double n1 = 1.0 / n;
  const double h1 = 1.0 / h;
  // standardiza design for more efficient computation
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  double my = arma::mean(Y);
  Y -= my;
  arma::vec der(n);
  arma::vec gradOld(p + 1), gradNew(p + 1);
  // initialization via expectile
  arma::vec beta = expectile(Z, Y, alpha, der, gradOld, gradNew, n, n1, tol, constTau, iteMax);
  arma::vec quant = {alpha};
  beta(0) = arma::as_scalar(arma::quantile(Y - Z.cols(1, p) * beta.rows(1, p), quant));
  // first step: a smoothed quantile estimator
  arma::vec res = Y - Z * beta;
  updateGauss(Z, res, der, gradOld, alpha, n1, h1);
  beta -= gradOld;
  arma::vec betaDiff = -gradOld;
  res -= Z * betaDiff;
  updateGauss(Z, res, der, gradNew, alpha, n1, h1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax) {
    double step = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      step = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -step * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    updateGauss(Z, res, der, gradNew, alpha, n1, h1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  // second step: an estimator for expected shortfall
  arma::vec w = Y - Z * beta;
  w = arma::min(w, arma::zeros(n));
  double mw = arma::mean(w);
  w -= mw;
  arma::vec theta = huberReg(Z, w, der, gradOld, gradNew, n, n1, tol, constTau, iteMax);
  // transform back to the original scale
  beta.rows(1, p) %= sx1;
  beta(0) += my - arma::as_scalar(mx * beta.rows(1, p));
  theta.rows(1, p) %= sx1;
  theta(0) += mw - arma::as_scalar(mx * theta.rows(1, p));
  theta = theta / alpha + beta;
  return Rcpp::List::create(Rcpp::Named("beta") = beta, Rcpp::Named("theta") = theta);
}

