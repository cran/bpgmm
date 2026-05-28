// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include "calculate_post_lambda_psi.h"
using namespace arma;

namespace {

bool constraint_matches(const arma::vec &constraint,
                        int lambda_constraint,
                        int psi_constraint,
                        int isotropic_constraint) {
  return constraint(0) == lambda_constraint &&
         constraint(1) == psi_constraint &&
         constraint(2) == isotropic_constraint;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List calculate_post_lambda_psi_native(int m,
                                   int p,
                                   Rcpp::S4 hparam,
                                   Rcpp::List cxy_list,
                                   Rcpp::S4 theta_y_list,
                                   arma::vec q_vec,
                                   arma::vec constraint) {

  validate_positive_int(m, "m");
  validate_positive_int(p, "p");
  validate_q_vec(q_vec, m);
  validate_constraint_vec(constraint);

  double alpha2 = get_positive_finite_slot(hparam, "alpha2");
  double bbeta  = get_positive_finite_slot(hparam, "bbeta");
  double delta  = get_positive_finite_slot(hparam, "delta");

  Rcpp::List Cxxk      = cxy_list["Cxxk"];
  Rcpp::List Cxyk      = cxy_list["Cxyk"];
  Rcpp::List Cyyk      = cxy_list["Cyyk"];
  Rcpp::List Cytytk    = cxy_list["Cytytk"];
  Rcpp::List Cxtytk    = cxy_list["Cxtytk"];
  Rcpp::List CxL1k     = cxy_list["CxL1k"];
  Rcpp::List Cxmyk     = cxy_list["Cxmyk"];

  arma::mat sumCxmyk  = cxy_list["sumCxmyk"];
  arma::mat sumCyyk   = cxy_list["sumCyyk"];

  Rcpp::List A   = cxy_list["A"];
  arma::vec nVec = cxy_list["nVec"];

  Rcpp::List M   = theta_y_list.slot("M");
  Rcpp::List psy = theta_y_list.slot("psy");

  if (Cxxk.size() < m || Cxyk.size() < m || Cyyk.size() < m ||
      Cytytk.size() < m || Cxtytk.size() < m || CxL1k.size() < m ||
      Cxmyk.size() < m || A.size() < m || M.size() < m || psy.size() < m) {
    Rcpp::stop("cxy_list and theta_y_list components must have length at least m");
  }
  if (nVec.n_elem < static_cast<arma::uword>(m)) {
    Rcpp::stop("nVec length must be at least m");
  }
  validate_finite_matrix(sumCxmyk, "sumCxmyk");
  validate_finite_matrix(sumCyyk, "sumCyyk");
  if (!nVec.is_finite()) {
    Rcpp::stop("nVec must contain only finite values");
  }

  Rcpp::List lambda(m);

  Rcpp::Environment mvtnorm = Rcpp::Environment::namespace_env("mvtnorm");
  Rcpp::Function rmvnorm = mvtnorm["rmvnorm"];

  Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
  Rcpp::Function kronecker = base["kronecker"];
  Rcpp::Function c = base["c"];
  Rcpp::Function matrix = base["matrix"];

  if (constraint_matches(constraint, 1, 1, 1)) {
    // Model CCC
    sumCxmyk = sumCxmyk.zeros();
    sumCyyk = sumCyyk.zeros();

    for (int k=0; k<m; ++k) {
      arma::mat Cxmyk_k = Cxmyk[k];
      arma::mat Cyyk_k = Cyyk[k];
      sumCxmyk = sumCxmyk + Cxmyk_k;
      arma::mat alpha2_eye(q_vec[k], q_vec[k], arma::fill::eye);
      alpha2_eye = alpha2 / m * alpha2_eye;  
      sumCyyk = sumCyyk + Cyyk_k + alpha2_eye;
    }


    
    for (int k=0; k<m; ++k) {
      if (k == 0) {

        Rcpp::NumericVector mean_vec = c(sumCxmyk * sumCyyk.i());

        Rcpp::NumericMatrix sigma_mat = kronecker(sumCyyk.i(), psy[k]);


        lambda[k] = rmvnorm(Rcpp::Named("n", 1),
                            Rcpp::Named("mean", mean_vec),
                            Rcpp::Named("sigma", sigma_mat));
        Rcpp::NumericMatrix lambdak = lambda[k];
        lambda[k] = matrix(lambdak, p, q_vec[k]);

      } else{
        lambda[k] = lambda[0];
      }
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    Rcpp::List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];

      arma::vec m_ka = m_k;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      lambda_ka.insert_cols(0, m_ka);
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    Rcpp::List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);

    for (int k=0; k<m; ++k) {


      shapePara += 0.5 * p * (nVec[k] + q_vec[k] / m + (2 * delta - 2) / (m * p) + 1);

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];


      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);


      bbeta_eye = 2 * bbeta/(m * p) * bbeta_eye;


      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka / m) * trans(tildaLambda_ka) + bbeta_eye;
      ratePara_vec += arma::diagvec(ratePara_k) * 0.5;

    }
    shapePara += 1;
    double ratePara = arma::sum(ratePara_vec);
    double scalePara = 1 / ratePara;
    double invpsy = sum(arma::randg( 1, distr_param(shapePara, scalePara) ));


    for (int k=0; k<m; ++k) {
      arma::mat post_psy_eye(p, p, arma::fill::eye);

      post_psy(k) = 1/invpsy * post_psy_eye;
    }
    Rcpp::List res = Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                                  Rcpp::Named("psy")    = post_psy);

    return(res);

  } else if (constraint_matches(constraint, 1, 1, 0)) {

    // Model CCU

    sumCxmyk = sumCxmyk.zeros();
    sumCyyk = sumCyyk.zeros();

    for (int k=0; k<m; ++k) {
      arma::mat Cxmyk_k = Cxmyk[k];
      arma::mat Cyyk_k = Cyyk[k];
      sumCxmyk = sumCxmyk + Cxmyk_k;
      arma::mat alpha2_eye(q_vec[k], q_vec[k], arma::fill::eye);
      alpha2_eye = alpha2 / m * alpha2_eye;  
      sumCyyk = sumCyyk + Cyyk_k + alpha2_eye;
    }

    for (int k=0; k<m; ++k) {
      if (k == 0) {

        Rcpp::NumericVector mean_vec = c(sumCxmyk * sumCyyk.i());

        Rcpp::NumericMatrix sigma_mat = kronecker(sumCyyk.i(), psy[k]);


        lambda[k] = rmvnorm(Rcpp::Named("n", 1),
                            Rcpp::Named("mean", mean_vec),
                            Rcpp::Named("sigma", sigma_mat));
        Rcpp::NumericMatrix lambdak = lambda[k];
        lambda[k] = matrix(lambdak, p, q_vec[k]);

      } else{
        lambda[k] = lambda[0];
      }
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    Rcpp::List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];

      arma::vec m_ka = m_k;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      lambda_ka.insert_cols(0, m_ka);
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    Rcpp::List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);

    for (int k=0; k<m; ++k) {


      shapePara += 0.5 * (nVec[k] + q_vec[k] / m + (2 * delta - 2)/m + 1);

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];


      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = (2 * bbeta / m) * bbeta_eye;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka/m) * trans(tildaLambda_ka) + bbeta_eye;
      ratePara_vec += arma::diagvec(ratePara_k) * 0.5;

    }
    shapePara += 1;
    arma::vec ratePara = ratePara_vec;
    arma::vec scalePara = 1 / ratePara;

    arma::vec invpsy(p);
    for (int j=0; j<p; ++j) {
      invpsy[j] = sum(arma::randg( 1, distr_param(shapePara, scalePara[j]) ));

    }


    for (int k=0; k<m; ++k) {
      // arma::mat post_psy_eye(p, p, arma::fill::eye);

      // post_psy(k) = trans(invpsy) * post_psy_eye;
      post_psy(k) = diagmat(1/invpsy);

    }
    Rcpp::List res = Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                                  Rcpp::Named("psy")    = post_psy);

    return(res);
  } else if (constraint_matches(constraint, 1, 0, 1)) {
    arma::mat Cxmyk_0 = Cxmyk[0];
    arma::mat Cyyk_0 = Cyyk[0];

    arma::mat sumPhiCxy(Cxmyk_0.n_rows, Cxmyk_0.n_cols, fill::zeros);
    arma::mat sumPhiCyy(Cyyk_0.n_rows, Cyyk_0.n_cols, fill::zeros);


    for (int k=0; k<m; ++k) {
      arma::mat psy_k = psy[k];
      Rcpp::NumericMatrix Cxmyk_k = Cxmyk[k];
      Rcpp::NumericMatrix Cyyk_k = Cyyk[k];

      arma::mat Cxmyk_ka = Rcpp::as<arma::mat>(Cxmyk_k);
      arma::mat Cyyk_ka = Rcpp::as<arma::mat>(Cyyk_k);


      arma::mat alpha2_eye(q_vec[k], q_vec[k], arma::fill::eye);
      alpha2_eye = (alpha2 / m) * alpha2_eye;


      sumPhiCxy += Cxmyk_ka / psy_k(1, 1);
      sumPhiCyy += ( Cyyk_ka + alpha2_eye ) / psy_k(1, 1);

    }

    for (int k=0; k<m; ++k) {
      if (k == 0) {

        Rcpp::NumericVector mean_vec = c(sumPhiCxy * sumPhiCyy.i());

        arma::mat p_eye(p, p, arma::fill::eye);
        Rcpp::NumericMatrix sigma_mat = kronecker(sumPhiCyy.i(), p_eye);


        lambda[k] = rmvnorm(Rcpp::Named("n", 1),
                            Rcpp::Named("mean", mean_vec),
                            Rcpp::Named("sigma", sigma_mat));
        Rcpp::NumericMatrix lambdak = lambda[k];
        lambda[k] = matrix(lambdak, p, q_vec[k]);

      } else{
        lambda[k] = lambda[0];
      }
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    Rcpp::List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];

      arma::vec m_ka = m_k;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      lambda_ka.insert_cols(0, m_ka);
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    Rcpp::List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);

    for (int k=0; k<m; ++k) {


      shapePara = 0.5 * p * (nVec[k] + q_vec[k] / m + (2 * delta - 2) / p + 1) + 1;

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];


      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * (bbeta / p) * bbeta_eye;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka) * trans(tildaLambda_ka) + bbeta_eye;
      ratePara_vec = arma::diagvec(ratePara_k) * 0.5;


      double ratePara = arma::sum(ratePara_vec);
      double scalePara = 1 / ratePara;
      double invpsy = sum(arma::randg( 1, distr_param(shapePara, scalePara) ));
      arma::mat post_psy_eye(p, p, arma::fill::eye);
      post_psy(k) = 1/invpsy * post_psy_eye;

    }



    Rcpp::List res = Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                                  Rcpp::Named("psy")    = post_psy);

    return(res);

  } else if (constraint_matches(constraint, 1, 0, 0)) {
    arma::mat sumVar;
    arma::mat B;

    for (int k=0; k<m; ++k) {
      arma::mat psy_k = psy[k];

      Rcpp::NumericMatrix Cxmyk_k = Cxmyk[k];
      Rcpp::NumericMatrix Cyyk_k = Cyyk[k];

      arma::mat Cxmyk_ka = Rcpp::as<arma::mat>(Cxmyk_k);
      arma::mat Cyyk_ka = Rcpp::as<arma::mat>(Cyyk_k);

      arma::mat alpha2_eye(q_vec[k], q_vec[k], arma::fill::eye);
      alpha2_eye = (alpha2 / m) * alpha2_eye;


      Rcpp::NumericMatrix sumVar_plus = kronecker(Cyyk_ka + alpha2_eye, psy_k.i());
      arma::mat sumVar_plusa = Rcpp::as<arma::mat>(sumVar_plus);
      if (k == 0) {
        sumVar = sumVar_plusa;
        B = psy_k.i() * Cxmyk_ka;
      } else {
        sumVar += sumVar_plusa;
        B += psy_k.i() * Cxmyk_ka;

      }

    }

    arma::mat lambdaVar = sumVar.i();
    Rcpp::NumericVector c_B = c(B);
    arma::vec c_Ba = Rcpp::as<arma::vec>(c_B);


    arma::mat lambdaMean = trans(c_Ba) * lambdaVar;

    for (int k=0; k<m; ++k) {
      if (k == 0) {
        lambda[k] = rmvnorm(Rcpp::Named("n", 1),
                            Rcpp::Named("mean", lambdaMean),
                            Rcpp::Named("sigma", lambdaVar));
        Rcpp::NumericMatrix lambdak = lambda[k];
        lambda[k] = matrix(lambdak, p, q_vec[k]);

      } else{
        lambda[k] = lambda[0];
      }
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    Rcpp::List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];

      arma::vec m_ka = m_k;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      lambda_ka.insert_cols(0, m_ka);
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    Rcpp::List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);

    for (int k=0; k<m; ++k) {


      shapePara = 0.5 * (nVec[k] + q_vec[k] / m + 2 * delta - 1) + 1;

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];


      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * bbeta * bbeta_eye;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka / m) * trans(tildaLambda_ka) + bbeta_eye;

      ratePara_k = arma::diagvec(ratePara_k) * 0.5;
      arma::vec scalePara = 1 / ratePara_k;

      arma::vec invpsy(p);
      for (int j=0; j<p; ++j) {


        invpsy[j] = sum(arma::randg( 1, distr_param(shapePara, scalePara[j]) ));;

      }

      //

      post_psy(k) = diagmat(1/invpsy);

    }

    // for (int k=0; k<m; ++k) {
    //   arma::mat post_psyk = post_psy(k);
    //
    //
    // }

    Rcpp::List res = Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                                  Rcpp::Named("psy")    = post_psy);

    return(res);

  } else if (constraint_matches(constraint, 0, 1, 1)) {
    for (int k=0; k<m; ++k) {


      Rcpp::NumericMatrix Cxmyk_k = Cxmyk[k];
      Rcpp::NumericMatrix Cyyk_k = Cyyk[k];

      arma::mat Cxmyk_ka = Rcpp::as<arma::mat>(Cxmyk_k);
      arma::mat Cyyk_ka = Rcpp::as<arma::mat>(Cyyk_k);
      arma::mat alpha2_eye(q_vec[k], q_vec[k], arma::fill::eye);

      alpha2_eye *= alpha2;

      Cyyk_ka += alpha2_eye;

      Rcpp::NumericVector mean_vec = c(Cxmyk_ka * Cyyk_ka.i());


      Rcpp::NumericMatrix sigma_mat = kronecker(Cyyk_ka.i(), psy[k]);


      lambda[k] = rmvnorm(Rcpp::Named("n", 1),
                          Rcpp::Named("mean", mean_vec),
                          Rcpp::Named("sigma", sigma_mat));
      Rcpp::NumericMatrix lambdak = lambda[k];
      lambda[k] = matrix(lambdak, p, q_vec[k]);
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    Rcpp::List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];

      arma::vec m_ka = m_k;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      lambda_ka.insert_cols(0, m_ka);
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    Rcpp::List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);

    for (int k=0; k<m; ++k) {


      shapePara += 0.5 * p * (nVec[k] + q_vec[k] + (2 * delta - 2)/(m * p) + 1);

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];


      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * (bbeta / (m * p)) * bbeta_eye;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka) * trans(tildaLambda_ka) + bbeta_eye;
      ratePara_vec += arma::diagvec(ratePara_k) * 0.5;

    }
    shapePara += 1;
    double ratePara = arma::sum(ratePara_vec);
    double scalePara = 1 / ratePara;
    double invpsy = sum(arma::randg( 1, distr_param(shapePara, scalePara) ));


    for (int k=0; k<m; ++k) {
      arma::mat post_psy_eye(p, p, arma::fill::eye);

      post_psy(k) = 1/invpsy * post_psy_eye;
    }
    Rcpp::List res = Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                                  Rcpp::Named("psy")    = post_psy);

    return(res);


  } else if (constraint_matches(constraint, 0, 1, 0)) {

    for (int k=0; k<m; ++k) {


      Rcpp::NumericMatrix Cxmyk_k = Cxmyk[k];
      Rcpp::NumericMatrix Cyyk_k = Cyyk[k];

      arma::mat Cxmyk_ka = Rcpp::as<arma::mat>(Cxmyk_k);
      arma::mat Cyyk_ka = Rcpp::as<arma::mat>(Cyyk_k);
      arma::mat alpha2_eye(q_vec[k], q_vec[k], arma::fill::eye);

      alpha2_eye *= alpha2;

      Cyyk_ka += alpha2_eye;

      Rcpp::NumericVector mean_vec = c(Cxmyk_ka * Cyyk_ka.i());


      Rcpp::NumericMatrix sigma_mat = kronecker(Cyyk_ka.i(), psy[k]);


      lambda[k] = rmvnorm(Rcpp::Named("n", 1),
                          Rcpp::Named("mean", mean_vec),
                          Rcpp::Named("sigma", sigma_mat));
      Rcpp::NumericMatrix lambdak = lambda[k];
      lambda[k] = matrix(lambdak, p, q_vec[k]);
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    Rcpp::List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];

      arma::vec m_ka = m_k;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      lambda_ka.insert_cols(0, m_ka);
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    Rcpp::List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);

    for (int k=0; k<m; ++k) {


      shapePara += 0.5 * (nVec[k] + q_vec[k] + (2 * delta - 2)/m + 1);

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];


      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * (bbeta / m) * bbeta_eye;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka) * trans(tildaLambda_ka) + bbeta_eye;
      ratePara_vec += arma::diagvec(ratePara_k) * 0.5;

    }
    shapePara += 1;

    arma::vec scalePara = 1 / ratePara_vec;
    arma::vec invpsy(p);
    for (int j=0; j<p; ++j) {

      invpsy[j] = sum(arma::randg( 1, distr_param(shapePara, scalePara[j]) ));;

    }

    for (int k=0; k<m; ++k) {

      post_psy(k) = diagmat(1/invpsy);
    }

    Rcpp::List res = Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                                  Rcpp::Named("psy")    = post_psy);

    return(res);

  } else if (constraint_matches(constraint, 0, 0, 1)) {


    for (int k=0; k<m; ++k) {


      Rcpp::NumericMatrix Cxmyk_k = Cxmyk[k];
      Rcpp::NumericMatrix Cyyk_k = Cyyk[k];

      arma::mat Cxmyk_ka = Rcpp::as<arma::mat>(Cxmyk_k);
      arma::mat Cyyk_ka = Rcpp::as<arma::mat>(Cyyk_k);
      arma::mat alpha2_eye(q_vec[k], q_vec[k], arma::fill::eye);

      alpha2_eye *= alpha2;

      Cyyk_ka += alpha2_eye;

      Rcpp::NumericVector mean_vec = c(Cxmyk_ka * Cyyk_ka.i());


      Rcpp::NumericMatrix sigma_mat = kronecker(Cyyk_ka.i(), psy[k]);


      lambda[k] = rmvnorm(Rcpp::Named("n", 1),
                          Rcpp::Named("mean", mean_vec),
                          Rcpp::Named("sigma", sigma_mat));
      Rcpp::NumericMatrix lambdak = lambda[k];
      lambda[k] = matrix(lambdak, p, q_vec[k]);
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    Rcpp::List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];

      arma::vec m_ka = m_k;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      lambda_ka.insert_cols(0, m_ka);
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    Rcpp::List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);

    for (int k=0; k<m; ++k) {

      shapePara = 0.5 * p * (nVec[k] + q_vec[k] +  (2 * delta - 2)/p + 1) + 1;

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];

      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * (bbeta / p) * bbeta_eye;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka) * trans(tildaLambda_ka) + bbeta_eye;
      ratePara_vec = arma::diagvec(ratePara_k) * 0.5;
      double ratePara = arma::sum(ratePara_vec);

      double scalePara = 1 / ratePara;
      double invpsy = sum(arma::randg( 1, distr_param(shapePara, scalePara) ));

      arma::mat post_psy_eye(p, p, arma::fill::eye);


      post_psy(k) = 1/invpsy * post_psy_eye;

    }

    Rcpp::List res = Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                                  Rcpp::Named("psy")    = post_psy);

    return(res);

  } else if (constraint_matches(constraint, 0, 0, 0)) {

    for (int k=0; k<m; ++k) {


      Rcpp::NumericMatrix Cxmyk_k = Cxmyk[k];
      Rcpp::NumericMatrix Cyyk_k = Cyyk[k];

      arma::mat Cxmyk_ka = Rcpp::as<arma::mat>(Cxmyk_k);
      arma::mat Cyyk_ka = Rcpp::as<arma::mat>(Cyyk_k);
      arma::mat alpha2_eye(q_vec[k], q_vec[k], arma::fill::eye);

      alpha2_eye *= alpha2;

      Cyyk_ka += alpha2_eye;

      Rcpp::NumericVector mean_vec = c(Cxmyk_ka * Cyyk_ka.i());


      Rcpp::NumericMatrix sigma_mat = kronecker(Cyyk_ka.i(), psy[k]);


      lambda[k] = rmvnorm(Rcpp::Named("n", 1),
                          Rcpp::Named("mean", mean_vec),
                          Rcpp::Named("sigma", sigma_mat));
      Rcpp::NumericMatrix lambdak = lambda[k];
      lambda[k] = matrix(lambdak, p, q_vec[k]);
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    Rcpp::List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];

      arma::vec m_ka = m_k;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      lambda_ka.insert_cols(0, m_ka);
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    Rcpp::List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);

    for (int k=0; k<m; ++k) {

      shapePara = 0.5 * (nVec[k] + q_vec[k] + 2 * delta - 1) + 1;

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];

      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * bbeta * bbeta_eye;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka) * trans(tildaLambda_ka) + bbeta_eye;
      Rcpp::NumericVector ratePara = c(arma::diagvec(ratePara_k) * 0.5);

      // arma:vec scalePara_vec = 1.0 / ratePara_vec;

      // arma::vec invpsy(p);
      // for (int j=0; j<p; ++j) {


        // invpsy[j] = sum(arma::randg( 1, distr_param(shapePara, scalePara_vec[j]) ));;

      // }

      // Rcpp::NumericVector invpsy = rgamma(Name("n", p),
      //                                     Name("shape", shapePara),
      //                                     Name("rate", ratePara));

      arma::vec invpsy(p);
      for (int j=0; j<p; ++j) {
        double ratePara_j = ratePara[j];
        double scalePara_j = 1.0/ratePara_j;
        invpsy[j] = sum(Rcpp::rgamma(1, shapePara, scalePara_j));

      }

      post_psy(k) = diagmat(1.0/invpsy);

    }

    Rcpp::List res = Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                                  Rcpp::Named("psy")    = post_psy);
    return(res);
  }

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                                Rcpp::Named("psy")    = psy);
  return(res);
}
