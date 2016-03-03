
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double SMOOTH = 0.00000001;
const double ALPHA = 0.8;
const double BETA = 0.6;

// function to solve basic Least square with network cohesion, when covariates are provided.
arma::mat RNC_LS_X(arma::mat X, arma::sp_mat L, arma::mat Y, double lambda, arma::sp_mat W ){

    int n = X.n_rows;
    int p = X.n_cols;
    arma::sp_mat I_sp = arma::speye<arma::sp_mat>(n,n);
    arma::mat I = arma::eye<arma::mat>(n,n);

    //arma::superlu_opts settings;
    //settings.symmetric = true;
    //arma::mat inverse = spsolve(I_sp + lambda * L, I, "superlu", settings);
    // uncomment above if superlu is available. We use the less efficient but by default available lapack.
    arma::mat inverse = spsolve(W + lambda * L, I,"lapack");
    arma::mat b = X.t()*(I-inverse)*Y;
    arma::mat Amat = X.t()*(I-W*inverse)*W*X;
    arma::mat beta = solve(Amat,b);
    arma::mat alpha = inverse*W*(Y-X*beta);
    arma::mat result = arma::join_cols(alpha,beta);

    return result;

}

arma::mat RNC_LS_noX(arma::sp_mat L, arma::mat Y, double lambda, arma::sp_mat W ){

    //arma::superlu_opts settings;
    //settings.symmetric = true;
    //arma::mat inverse = spsolve(I_sp + lambda * L, I, "superlu", settings);
    // uncomment above if superlu is available. We use the less efficient but by default available lapack.
    arma::mat alpha = spsolve(W + lambda * L, W*Y,"lapack");
    return alpha;

}



// The currently available 'as' functions from R sparse matrix to sp_mat all have some numerical issues.
// So we currently only support using the normal matrix as the format for Laplacian. Notice that this
// sacrifices some speed.
// [[Rcpp::export]]
extern "C" SEXP rnc_solver_X(SEXP X, SEXP L, SEXP Y, SEXP lambda, SEXP W){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Ldsmat=as<arma::mat>(L);
    arma::sp_mat Lspmat=arma::sp_mat(Ldsmat);
    arma::mat Ymat=as<arma::mat>(Y);
    arma::mat Wdsmat=as<arma::mat>(W);
    arma::sp_mat Wspmat=arma::sp_mat(Wdsmat);
    double lambda_num= as<double>(lambda);

    arma::mat result = RNC_LS_X(Xmat,Lspmat,Ymat,lambda_num,Wspmat);
    return(wrap(result));
}
// [[Rcpp::export]]
extern "C" SEXP rnc_solver_noX(SEXP L, SEXP Y, SEXP lambda, SEXP W){

    arma::mat Ldsmat=as<arma::mat>(L);
    arma::sp_mat Lspmat=arma::sp_mat(Ldsmat);
    arma::mat Ymat=as<arma::mat>(Y);
    arma::mat Wdsmat=as<arma::mat>(W);
    arma::sp_mat Wspmat=arma::sp_mat(Wdsmat);
    double lambda_num= as<double>(lambda);

    arma::mat result = RNC_LS_noX(Lspmat,Ymat,lambda_num,Wspmat);
    return(wrap(result));
}



// function to solve basic Least square with network cohesion, when covariates are provided.
arma::mat RNC_LS_Naive(arma::mat X, arma::sp_mat L, arma::mat Y, double lambda, arma::sp_mat W ){

    int n = X.n_rows;
    int p = X.n_cols;
    arma::sp_mat I_sp = arma::speye<arma::sp_mat>(n,n);
    arma::mat I = arma::eye<arma::mat>(n,n);
    arma::mat X_tilde = arma::join_rows(I,X);
    arma::sp_mat X_tilde_sp = arma::sp_mat(X_tilde);
    arma::sp_mat core = X_tilde_sp.t()*W*X_tilde_sp+lambda*L;

    arma::mat result = spsolve(core,X_tilde.t()*W*Y,"lapack");

    return result;

}

// [[Rcpp::export]]
extern "C" SEXP rnc_solver_naive(SEXP X, SEXP L, SEXP Y, SEXP lambda, SEXP W){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Ldsmat=as<arma::mat>(L);
    int n = Xmat.n_rows;
    int p = Xmat.n_cols;
    arma::mat Omega = arma::zeros<arma::mat>(n+p,n+p);
    Omega(arma::span(0,n-1),arma::span(0,n-1)) = Ldsmat;
    arma::sp_mat Lspmat=arma::sp_mat(Omega);
    arma::mat Ymat=as<arma::mat>(Y);
    arma::mat Wdsmat=as<arma::mat>(W);
    arma::sp_mat Wspmat=arma::sp_mat(Wdsmat);
    double lambda_num= as<double>(lambda);

    arma::mat result = RNC_LS_Naive(Xmat,Lspmat,Ymat,lambda_num,Wspmat);
    return(wrap(result));
}


double loglike_bernoulli(arma::mat Y, arma::mat P){
    double result = 0;
    int n = Y.n_rows;
    for (int i = 0; i<n; i++) {
        if (P(i,0) < SMOOTH) {
            P(i,0) = SMOOTH;
        }
        if (P(i,0)>1-SMOOTH) {
            P(i,0) = 1-SMOOTH;
        }
        result += Y(i,0)*log(P(i,0)) + (1-Y(i,0))*log(1-P(i,0));
    }
    return result;
}

arma::mat logit_p(arma::mat eta){
    int n = eta.n_rows;
    arma::mat expeta = exp(eta);
    arma::mat ones = arma::ones(n,1);
    arma::mat result = expeta/(expeta+ones);
    return result;
}


// function to solve logistic regression with network cohesion, when covariates are provided.

arma::mat RNC_Logit_X(arma::mat X, arma::sp_mat L, arma::mat Y, double lambda, arma::mat theta_init, double tol, int max_iter, bool verbose){

    // begin newton method
    int n = X.n_rows;
    int p = X.n_cols;
    int iter = 0;
    double err = 0;
    arma::mat I = arma::eye<arma::mat>(n,n);

    arma::mat X_tilde = arma::join_rows(I,X);
    arma::mat eta, theta_old, theta_new;
    double ell;
    theta_old = theta_init;
    eta = X_tilde * theta_old;
    arma::mat one_n = arma::ones(n,1);
    arma::mat P = logit_p(eta);
    ell = loglike_bernoulli(Y,P);
    arma::mat residual = Y - P;
    arma::mat w_vec = P%(one_n - P);
    arma::mat z = eta + residual/w_vec;
    arma::sp_mat W = arma::sp_mat(n,n);
    W.diag() = w_vec;
    bool converge = false;
    while (!converge) {
        iter += 1;
        theta_new = RNC_LS_Naive(X, L, z, lambda, W );
        err = arma::norm(theta_new - theta_old,2)/(arma::norm(theta_old,2)+SMOOTH);
        if (err < tol) {
            converge = true;
        }
        theta_old = theta_new;
        eta = X_tilde * theta_old;
        P = logit_p(eta);
        ell = loglike_bernoulli(Y,P);
        residual = Y - P;
        w_vec = P%(one_n - P);
        z = eta + residual/w_vec;
        W.diag() = w_vec;
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteraction reached before converge!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finish iteration " << iter << " with loglikelihood " << ell << std::endl;
        }
        if ((P.max()>1-SMOOTH) || (P.min() < SMOOTH)) {
            converge = true;
        }
    }

    return theta_old;

}


// [[Rcpp::export]]
extern "C" SEXP rnc_logistic_fit(SEXP X, SEXP L, SEXP Y, SEXP lambda, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Xmat=as<arma::mat>(X);
    arma::mat Ldsmat=as<arma::mat>(L);
    int n = Xmat.n_rows;
    int p = Xmat.n_cols;
    arma::mat Omega = arma::zeros<arma::mat>(n+p,n+p);
    Omega(arma::span(0,n-1),arma::span(0,n-1)) = Ldsmat;
    arma::sp_mat Lspmat=arma::sp_mat(Omega);
    //arma::sp_mat Lspmat=arma::sp_mat(Ldsmat);
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    arma::mat result = RNC_Logit_X(Xmat,Lspmat,Ymat,lambda_num,thetamat,tol_num,iter_max_num,verbose_ind);
    return(wrap(result));
}

// function to solve logistic regression with network cohesion, when covariates are not provided.

arma::mat RNC_Logit_noX(arma::sp_mat L, arma::mat Y, double lambda, arma::mat theta_init, double tol, int max_iter, bool verbose){

    // begin newton method
    int n = Y.n_rows;
    int iter = 0;
    double err = 0;
    arma::mat I = arma::eye<arma::mat>(n,n);

    arma::mat X_tilde = I;
    arma::mat eta, theta_old, theta_new;
    double ell;
    theta_old = theta_init;
    eta = X_tilde * theta_old;
    arma::mat one_n = arma::ones(n,1);
    arma::mat P = logit_p(eta);
    ell = loglike_bernoulli(Y,P);
    arma::mat residual = Y - P;
    arma::mat w_vec = P%(one_n - P);
    arma::mat z = eta + residual/w_vec;
    arma::sp_mat W = arma::sp_mat(n,n);
    W.diag() = w_vec;
    bool converge = false;
    while (!converge) {
        iter += 1;
        theta_new = RNC_LS_noX(L, z, lambda, W );
        err = arma::norm(theta_new - theta_old,2)/(arma::norm(theta_old,2)+SMOOTH);
        if (err < tol) {
            converge = true;
        }
        theta_old = theta_new;
        eta = X_tilde * theta_old;
        P = logit_p(eta);
        ell = loglike_bernoulli(Y,P);
        residual = Y - P;
        w_vec = P%(one_n - P);
        z = eta + residual/w_vec;
        W.diag() = w_vec;
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteraction reached before converge!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Finish iteration " << iter << " with loglikelihood " << ell << std::endl;

        }
        if ((P.max()>1-SMOOTH) || (P.min() < SMOOTH)) {
            converge = true;
        }
    }

    return theta_old;

}


// [[Rcpp::export]]
extern "C" SEXP rnc_logistic_fit_noX(SEXP L, SEXP Y, SEXP lambda, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Ldsmat=as<arma::mat>(L);
    int n = Ldsmat.n_rows;
    arma::mat Omega = Ldsmat;
    arma::sp_mat Lspmat=arma::sp_mat(Omega);
    arma::mat Ymat=as<arma::mat>(Y);
    double lambda_num= as<double>(lambda);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    arma::mat result = RNC_Logit_noX(Lspmat,Ymat,lambda_num,thetamat,tol_num,iter_max_num,verbose_ind);
    return(wrap(result));
}


double cox_partialloglike(arma::mat delta, arma::mat eta, arma::mat PairMat){
    int n = eta.n_rows;
    arma::mat expeta = exp(eta);
    double result = 0;
    for (int i=0; i<n; i++) {
        double tmp = 0;
        if (delta(i,0)==1) {
            tmp = eta(i,0);
            double tmpsum = 0;
            arma::mat summat = PairMat.row(i)*expeta;
            tmpsum = summat(0,0);
            tmp = tmp - log(tmpsum);
            result = result + tmp;
        }
    }
    return result;
}

void cox_newton_terms_calculation(arma::mat &grad, arma::mat &hess, arma::mat &response_g, arma::mat delta, arma::mat eta, arma::mat PairMat){

    int n = eta.n_rows;
    arma::mat expeta = exp(eta);
    
    for (int i = 0; i<n; i++) {
        double tmp_grad=delta(i,0);
        double tmp_hess = 0;
        for (int j=0; j<n; j++) {
            if ((int)PairMat(j,i)==1 && (int)delta(j,0)==1) {
                arma::mat denom_mat = PairMat.row(j)*expeta;
                double denom = denom_mat(0,0);
                tmp_grad = tmp_grad - expeta(i,0)/(denom+SMOOTH);
                tmp_hess = tmp_hess + expeta(i,0)*(denom-expeta(i,0))/(denom*denom+SMOOTH);
            }
        }
        grad(i,0) = tmp_grad;
        hess(i,0) = tmp_hess;
        response_g(i,0) = eta(i,0) + grad(i,0)/(hess(i,0)+SMOOTH);
    }
    return;
}



// function to fit Cox's proportional hazard model with network cohesion, when covariates are provided. The Breslow approximation is used for ties. The detailed algorithm formulation can be found in Generalized Additive Models, by Hastie and Tibshirani
// There are many operations to calculate sets, which can be done more efficiently in R. So we will assume such sets as the input.

arma::mat RNC_Cox_X(arma::mat X, arma::sp_mat L, arma::mat Y, arma::mat delta, double lambda, arma::mat theta_init, double tol, int max_iter, bool verbose){

    // begin newton method
    int n = X.n_rows;
    int p = X.n_cols;
    
    // First find the hazard set relationship
    // For each row i, the jth position is 1 if yj >= yi;
    // For each column j, the ith position is 1 if yi <= yj;
    arma::mat PairMat = arma::zeros<arma::mat>(n,n);
    for (int i=0; i<n; i++) {
        double dr = 0;
        for (int j=0; j<n; j++) {
            if (Y(i,0)<=Y(j,0)) {
                PairMat(i,j) = 1;
            }
        }
        
    }
    //Rcout << "Finish pairwise comparison " << std::endl;
    
    arma::mat grad(n,1);
    arma::mat hess=arma::zeros<arma::mat>(n,1);
    arma::mat response_g(n,1);
    
    
    int iter = 0;
    double err = 0;
    arma::mat I = arma::eye<arma::mat>(n,n);

    arma::mat X_tilde = arma::join_rows(I,X);
    arma::mat eta, theta_old, theta_new;
    double ell;
    theta_old = theta_init;
    eta = X_tilde * theta_old;
    arma::mat one_n = arma::ones(n,1);
    // calculate terms needed for newton
    //Rcout << "Begin first newton calculation " << std::endl;
    cox_newton_terms_calculation(grad, hess, response_g, delta, eta, PairMat);
    
    arma::sp_mat H = arma::sp_mat(n,n);
    H.diag() = hess;
    bool converge = false;
    while (!converge) {
        iter += 1;
        theta_new = RNC_LS_Naive(X, L, response_g, lambda, H );
        err = arma::norm(theta_new - theta_old,2)/(arma::norm(theta_old,2)+SMOOTH);
        if (err < tol) {
            converge = true;
        }
        theta_old = theta_new;
        eta = X_tilde * theta_old;
        ell = cox_partialloglike(delta, eta, PairMat);
        cox_newton_terms_calculation(grad, hess, response_g, delta, eta, PairMat);
        H.diag() = hess;
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteraction reached before converge!" << std::endl;
            }

            break;
        }
        if (verbose) {
            Rcout << "Begin iteration " << iter << " with partial loglikelihood " << ell << std::endl;
        }

    }

    return theta_old;

}


// [[Rcpp::export]]
extern "C" SEXP rnc_cox_fit(SEXP X, SEXP L, SEXP Y, SEXP delta, SEXP lambda, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    //Rcout << "Begin converting objects " << std::endl;
    arma::mat Xmat=as<arma::mat>(X);
    //Rcout << "Begin converting L " << std::endl;
    arma::mat Ldsmat=as<arma::mat>(L);
    int n = Xmat.n_rows;
    int p = Xmat.n_cols;
    arma::mat Omega = arma::zeros<arma::mat>(n+p,n+p);
    //Rcout << "Begin incorporating L " << std::endl;
    Omega(arma::span(0,n-1),arma::span(0,n-1)) = Ldsmat;
    arma::sp_mat Lspmat=arma::sp_mat(Omega);
    //arma::sp_mat Lspmat=arma::sp_mat(Ldsmat);
    arma::mat Ymat=as<arma::mat>(Y);
    arma::mat deltamat=as<arma::mat>(delta);
    double lambda_num= as<double>(lambda);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    //Rcout << "Begin C calculation " << std::endl;
    arma::mat result = RNC_Cox_X(Xmat,Lspmat,Ymat,deltamat,lambda_num,thetamat,tol_num,iter_max_num,verbose_ind);
    return(wrap(result));
}




arma::mat RNC_Cox_noX(arma::sp_mat L, arma::mat Y, arma::mat delta, double lambda, arma::mat theta_init, double tol, int max_iter, bool verbose){
    
    // begin newton method
    int n = Y.n_rows;
    
    // First find the hazard set relationship
    // For each row i, the jth position is 1 if yj >= yi;
    // For each column j, the ith position is 1 if yi <= yj;
    arma::mat PairMat = arma::zeros<arma::mat>(n,n);
    for (int i=0; i<n; i++) {
        double dr = 0;
        for (int j=0; j<n; j++) {
            if (Y(i,0)<=Y(j,0)) {
                PairMat(i,j) = 1;
            }
        }
        
    }
    //Rcout << "Finish pairwise comparison " << std::endl;
    
    arma::mat grad(n,1);
    arma::mat hess=arma::zeros<arma::mat>(n,1);
    arma::mat response_g(n,1);
    
    
    int iter = 0;
    double err = 0;
    arma::mat eta, theta_old, theta_new;
    double ell;
    theta_old = theta_init;
    eta = theta_old;
    // calculate terms needed for newton
    //Rcout << "Begin first newton calculation " << std::endl;
    cox_newton_terms_calculation(grad, hess, response_g, delta, eta, PairMat);
    
    arma::sp_mat H = arma::sp_mat(n,n);
    H.diag() = hess;
    bool converge = false;
    while (!converge) {
        iter += 1;
        theta_new = RNC_LS_noX(L, response_g, lambda, H );
        err = arma::norm(theta_new - theta_old,2)/(arma::norm(theta_old,2)+SMOOTH);
        if (err < tol) {
            converge = true;
        }
        theta_old = theta_new;
        eta = theta_old;
        ell = cox_partialloglike(delta, eta, PairMat);
        cox_newton_terms_calculation(grad, hess, response_g, delta, eta, PairMat);
        H.diag() = hess;
        if(iter == max_iter){
            if (verbose) {
                Rcout << "Maximum iteraction reached before converge!" << std::endl;
            }
            
            break;
        }
        if (verbose) {
            Rcout << "Begin iteration " << iter << " with partial loglikelihood " << ell << std::endl;
        }
        
    }
    
    return theta_old;
    
}



// [[Rcpp::export]]
extern "C" SEXP rnc_cox_fit_noX(SEXP L, SEXP Y, SEXP delta, SEXP lambda, SEXP theta_init, SEXP tol, SEXP max_iter, SEXP verbose){
    arma::mat Ldsmat=as<arma::mat>(L);
    int n = Ldsmat.n_rows;
    int p = 0;
    arma::mat Omega = arma::zeros<arma::mat>(n+p,n+p);
    //Rcout << "Begin incorporating L " << std::endl;
    Omega(arma::span(0,n-1),arma::span(0,n-1)) = Ldsmat;
    arma::sp_mat Lspmat=arma::sp_mat(Omega);
    //arma::sp_mat Lspmat=arma::sp_mat(Ldsmat);
    arma::mat Ymat=as<arma::mat>(Y);
    arma::mat deltamat=as<arma::mat>(delta);
    double lambda_num= as<double>(lambda);
    double tol_num= as<double>(tol);
    bool verbose_ind= as<bool>(verbose);
    int iter_max_num= as<int>(max_iter);
    arma::mat thetamat=as<arma::mat>(theta_init);
    //Rcout << "Begin C calculation " << std::endl;
    arma::mat result = RNC_Cox_noX(Lspmat,Ymat,deltamat,lambda_num,thetamat,tol_num,iter_max_num,verbose_ind);
    return(wrap(result));
}




double Cox_Partial(arma::mat eta, arma::mat Y, arma::mat delta){
    
    // begin newton method
    int n = Y.n_rows;
    
    // First find the hazard set relationship
    // For each row i, the jth position is 1 if yj >= yi;
    // For each column j, the ith position is 1 if yi <= yj;
    arma::mat PairMat = arma::zeros<arma::mat>(n,n);
    for (int i=0; i<n; i++) {
        double dr = 0;
        for (int j=0; j<n; j++) {
            if (Y(i,0)<=Y(j,0)) {
                PairMat(i,j) = 1;
            }
        }
        
    }
    double ell;
    ell = cox_partialloglike(delta, eta, PairMat);

    return ell;
    
}




// [[Rcpp::export]]
extern "C" SEXP cox_pll(SEXP eta, SEXP Y, SEXP delta){
    //Rcout << "Begin converting objects " << std::endl;
    arma::mat etamat=as<arma::mat>(eta);
    //arma::sp_mat Lspmat=arma::sp_mat(Ldsmat);
    arma::mat Ymat=as<arma::mat>(Y);
    arma::mat deltamat=as<arma::mat>(delta);
   //arma::mat result = Ldsmat(ind,col_ind);
    double result = Cox_Partial(etamat, Ymat, deltamat);
   return(wrap(result));
}

