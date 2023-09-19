//-----------------------------------------#
// This code has been done by Adam N. Smith and downloaded on 05/24/2023 from 
// https://github.com/adam-n-smith/demand-models-random-partitions
//-----------------------------------------#

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// --------------------------------------------------------------- //
// SECTION 1: core functions
// --------------------------------------------------------------- //

// compute similarity matrix
// [[Rcpp::export]]
mat compsm(mat const& draws) {
  int n = draws.n_cols;
  mat sm = zeros<mat>(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      uvec wch = find(draws.col(i)==draws.col(j));
      sm(i,j) = wch.size();
    }
  }
  return sm/draws.n_rows;
}

// function for random sampling from a categorical distribution
int rcat(vec pvec) {
  if (sum(pvec) != 1) pvec = pvec / sum(pvec);
  vec bins = cumsum(pvec);
  double u = as_scalar(randu(1));
  int draw = min(find(u<bins));
  return draw;
}

// density of non-standardized student-t
double dst(double const& x, double const& nu, double const& mu, double const& var){
  double constant = 0.5*nu*log(nu)+lgamma((nu+1)/2)-0.5*log(M_PI)-lgamma(nu/2);
  return exp(constant-log(sqrt(var))-((nu+1)/2.0)*log(nu+(x-mu)*(x-mu)/var));
}

// function to sample from Wishart distribution
List rwishart(int const& nu, mat const& V){
  int m = V.n_rows;
  mat T = zeros(m,m);
  for(int i=0; i<m; i++) {
    T(i,i) = sqrt(rchisq(1,nu-i)[0]);
  }
  for(int j=0; j<m; j++) {  
    for(int i=j+1; i<m; i++) {    
      T(i,j) = rnorm(1)[0]; 
    }
  }
  mat C = trans(T)*chol(V);
  mat CI = solve(trimatu(C),eye(m,m)); 
  
  return List::create(
    Named("W") = trans(C) * C,
    Named("IW") = CI * trans(CI),
    Named("C") = C,
    Named("CI") = CI);
}

// function to evaluate log multivariate normal density (where root'root=Sigma^-1)
double lndmvn(vec const& y, vec const& mu, mat const& root){
  vec z = root*(y-mu);
  return -0.5*mu.size()*log(2*M_PI)+sum(log(diagvec(root)))-.5*as_scalar(trans(z)*z);
}

// --------------------------------------------------------------- //
// SECTION 2: functions for LSP, LSPx, and DP distributions
// --------------------------------------------------------------- //

// random draw from DP partition distribution
// [[Rcpp::export]]
mat rdp(int const& n, double const& alpha){
  
  uvec wch;
  
  // vector to store draws
  vec pi = zeros<vec>(n);
  pi(0) = 1;
  
  for (int i=1;i<n;i++){
    
    // number of groups
    vec uniq = unique(pi.subvec(0,i-1));
    int J = uniq.n_rows;
    
    // counts matrix
    mat C = zeros<mat>(J+1,2);
    C(span(0,J-1),0) = uniq;
    
    // probability vector
    vec pvec = zeros<vec>(J+1);
    
    // fill in counts matrix and probability vector
    for (int j=0;j<J;j++){
      wch = find(pi==C(j,0));
      C(j,1) = wch.n_elem;
      pvec(j) = C(j,1)/(alpha+i-1);
    }
    pvec(J) = alpha/(alpha+i-1);
    
    // categorical variable with probability pvec
    int draw = rcat(pvec/sum(pvec));
    
    // draw of pi
    if(draw==J){
      pi(i) = J+1;
    }
    else {pi(i) = C(draw,0);
    }
  }
  
  return(pi);
}

// density of DP partition distribution
// [[Rcpp::export]]
double ddp(vec const& pi,double const& alpha){
  
  uvec wch;
  
  // vector to store draws 
  int n = pi.n_elem;
  
  // matrix of probabilities
  mat out_pvec = zeros<mat>(n,n);
  out_pvec(0,0) = 1;
  
  for (int i=1;i<n;i++){
    
    // number of groups
    vec uniq = unique(pi.subvec(0,i-1));
    int J = uniq.n_rows;
    
    // counts matrix
    mat C = zeros<mat>(J+1,2);
    C(span(0,J-1),0) = uniq;
    
    // probability vector
    vec pvec = zeros<vec>(J+1);
    
    // fill in counts matrix and probability vector
    for (int j=0;j<J;j++){
      wch = find(pi==C(j,0));
      C(j,1) = wch.n_elem;
      pvec(j) = C(j,1)/(alpha+i-1);
    }
    pvec(J) = alpha/(alpha+i-1);
    int end = pvec.n_rows-1;
    vec pnorm = pvec/sum(pvec);
    out_pvec(i,span(0,end)) = pnorm.t();
  }
  
  vec out = zeros<vec>(n);
  for(int i=0;i<n;i++){
    out(i) = out_pvec(i,pi(i)-1);
  }
  return(prod(out));
}

// similarity function
double similarity(int const& si,vec const& Sk,int const& C,double const& tau){
  
  // create tau vector
  vec tauvec = zeros<vec>(C);
  tauvec.fill(tau);
  
  // count number of elements in Sk equal to si
  vec nSkc = zeros<vec>(C);
  for(int c=0;c<C;c++){
    uvec count = find(Sk == c+1);
    nSkc(c) = count.n_elem;
  }
  
  // create taustar
  vec tau_star = tauvec + nSkc;
  double out = tau_star(si-1)/sum(tau_star);
  
  return out;
}

// similarity function w/ covariates
// [[Rcpp::export]]
double similarityx(int const& si,vec const& Sk,int const& C,double const& tau,vec const& xi,mat const& Xk,vec const& kappa0,double const& wk){
  
  // location-scale
  vec tauvec = zeros<vec>(C);
  tauvec.fill(tau);
  vec nSkc = zeros<vec>(C);
  for(int c=0;c<C;c++){
    uvec count = find(Sk == c+1);
    nSkc(c) = count.n_elem;
  }
  vec tau_star = tauvec + nSkc;
  double sim = tau_star(si-1)/sum(tau_star);
  
  // covariates
  int p = xi.n_elem;
  vec simx = ones<vec>(p);
  double nki,kappatilde,mutilde,betatilde,var;
  double mu0 = 0;
  double alpha0 = 0.5;
  double beta0 = 0.5;
  if (wk==1)nki = Xk.n_cols;
  else nki = 0;
  vec Xkbar = mean(Xk,1);
  double alphatilde = alpha0 + 0.5*nki;
  for(int j=0;j<p;j++){
    kappatilde = kappa0(j) + nki;
    mutilde = (kappa0(j)*mu0 + nki*Xkbar(j))/(kappa0(j)+nki);
    betatilde = beta0 + wk*0.5*(xi(j)-Xkbar(j))*(xi(j)-Xkbar(j)) + 0.5*nki*kappa0(j)/(kappa0(j)+nki)*(Xkbar(j)-mu0)*(Xkbar(j)-mu0);
    var = betatilde*(kappatilde+1)/(alphatilde*kappatilde);
    simx(j) = dst(xi(j),2*alphatilde,mutilde,var);
  }
  
  return sim*prod(simx);
}

// random draw from LSP distribution
// [[Rcpp::export]]
mat rlsp(int const& R,vec const& rho,double const& tau){
  
  mat C;
  vec uniq,pvec,pvecnorm;
  uvec wch;
  int r, i,Ki,Ci,draw;
  int N = rho.n_elem;
  mat pidraw = ones<mat>(R,N);
  vec pi = ones<vec>(N);
  
  for(r=0;r<R;r++){
    vec vCi = vec(1);
    
    for (i=1;i<N;i++){
      
      // number of groups
      uniq = unique(pi(span(0,i-1)));
      Ki = uniq.n_elem;
      
      // number of groups in rho
      Ci = max(rho(span(0,i-1)));
      vCi.fill(Ci);
      
      // probability vector
      pvec = zeros<vec>(Ki+1);
      
      // fill in counts matrix and probability vector
      for (int j=0;j<Ki;j++){
        wch = find(pi(span(0,i-1))==uniq(j));
        pvec(j) = similarity(rho(i),rho(wch),Ci+1,tau);
      }
      pvec(Ki) = similarity(rho(i),vCi+1,Ci+1,tau);
      pvecnorm = pvec/sum(pvec);
      draw = rcat(pvecnorm);
      
      // draw of gi
      if(draw==Ki){
        pi(i) = Ki+1;
      }
      else {pi(i) = uniq(draw);}
    }
    pidraw.row(r) = trans(pi);
  }
  
  return pidraw;
}

// random draw from LSPx distribution
// [[Rcpp::export]]
mat rlspx(int const& R, vec const& rho,double const& tau, mat const& X, vec const& kappa0){
  
  mat C;
  vec uniq,pvec;
  uvec wch;
  int i,j,r,Ki,Ci,draw;
  int N = rho.n_elem;
  vec pi = ones<vec>(N);
  vec vCi = vec(1);
  mat pidraw = ones<mat>(R,N);
  
  for (r=0;r<R;r++){
    
    for (i=1;i<N;i++){
      
      // number of groups
      uniq = unique(pi(span(0,i-1)));
      Ki = uniq.n_elem;
      
      // number of groups in rho
      Ci = max(rho(span(0,i-1)));
      vCi.fill(Ci);
      
      // probability vector
      pvec = zeros<vec>(Ki+1);
      
      // fill in counts matrix and probability vector
      for (j=0;j<Ki;j++){
        wch = find(pi(span(0,i-1))==uniq(j));
        pvec(j) = similarityx(rho(i),rho(wch),Ci+1,tau,trans(X.row(i)),trans(X.rows(wch)),kappa0,1);
      }
      pvec(Ki) = similarityx(rho(i),vCi+1,Ci+1,tau,trans(X.row(i)),zeros<mat>(X.n_cols,X.n_cols),kappa0,0);
      draw = rcat(pvec/sum(pvec));
      
      // draw of gi
      if(draw==Ki){
        pi(i) = Ki+1;
      }
      else {pi(i) = uniq(draw);}
    }
    pidraw.row(r) = trans(pi);
  }
  
  return pidraw;
}

// random draw from LSP distribution starting at item m
// [[Rcpp::export]]
vec rblocklsp(vec const& rho,double const& tau,int start,int end){
  
  vec pi = rho;
  int N = rho.n_elem;
  vec vCi = vec(1);
  if(start>1){
    start = start - 1;
  }
  
  for (int i=start;i<end;i++){
    
    // number of groups
    vec uniq = unique(pi(span(0,i-1)));
    int Ki = uniq.n_elem;
    
    // number of groups in rho
    int Ci = max(rho(span(0,i-1)));
    vCi.fill(Ci);
    
    // probability vector
    vec pvec = zeros<vec>(Ki+1);
    
    // fill in counts matrix and probability vector
    for (int j=0;j<Ki;j++){
      uvec wch = find(pi(span(0,i-1))==uniq(j));
      pvec(j) = similarity(rho(i),rho(wch),Ci+1,tau);
    }
    pvec(Ki) = similarity(rho(i),vCi+1,Ci+1,tau);
    vec pvecnorm = pvec/sum(pvec);
    int draw = rcat(pvecnorm);
    
    // draw of gi
    if(draw==Ki){
      pi(i) = Ki+1;
    }
    else {pi(i) = uniq(draw);
    }
  }
  
  // reorder
  if(end<N){
    // case (ii)
    vec notblock = rho(span(end,N-1));
    int maxpiblock = max(pi(span(0,end-1)));
    int maxrhoblock = max(rho(span(0,end-1)));
    if(maxpiblock>maxrhoblock){
      if(max(notblock)>maxpiblock){
        uvec wch = find(notblock>maxrhoblock);
        notblock(wch) = notblock(wch)+1;
        pi(span(end,N-1)) = notblock;
      }
    }
    
    // case (i)
    if(maxrhoblock>maxpiblock){
      for(int i=end;i<N;i++){
        if(pi(i)>=maxpiblock){
          if(i==end){
            pi(i) = max(pi(span(0,i-1)))+1;
          }
          else{
            vec pii = pi(span(end,i-1));
            vec rhoi = rho(span(end,i-1));
            if(any(rhoi==rho(i))){
              pi(i) = min(pii(find(rhoi==rho(i))));
            }
            else{
              pi(i) = max(pi(span(0,i-1)))+1;
            } 
          }
        }
      }
    }
    
  }
  return pi;
}

// evaluate PMF of LSP distribution
// [[Rcpp::export]]
double dlsp(vec const& pi,vec const& rho,double const& tau){
  
  mat C;
  vec uniq,vCi,pvec,pvecnorm;
  uvec wch;
  int Ki,Ci;
  
  int N = rho.n_elem;
  vec lprob = zeros<vec>(N);
  vCi = vec(1);
  
  for (int i=1;i<N;i++){
    
    // number of groups
    uniq = unique(pi(span(0,i-1)));
    Ki = uniq.n_elem;
    
    // number of groups in rho
    Ci = max(rho(span(0,i-1)));
    vCi.fill(Ci);
    
    // probability vector
    pvec = zeros<vec>(Ki+1);
    
    // fill in counts matrix and probability vector
    for (int j=0;j<Ki;j++){
      wch = find(pi(span(0,i-1))==uniq(j));
      pvec(j) = similarity(rho(i),rho(wch),Ci+1,tau);
    }
    pvec(Ki) = similarity(rho(i),vCi+1,Ci+1,tau);
    pvecnorm = pvec/sum(pvec);
    lprob(i) = log(pvecnorm(pi(i)-1));
  }
  return exp(sum(lprob));
}

// evaluate PMF of LSPx distribution
// [[Rcpp::export]]
double dlspx(vec const& pi,vec const& rho,double const& tau,mat const& X, vec const& kappa0){
  
  mat C;
  vec uniq,vCi,pvec,pvecnorm;
  uvec wch;
  int i,j,Ki,Ci;
  
  int N = rho.n_elem;
  vec lprob = zeros<vec>(N);
  vCi = vec(1);
  
  for (i=1;i<N;i++){
    
    // number of groups
    uniq = unique(pi(span(0,i-1)));
    Ki = uniq.n_elem;
    
    // number of groups in rho
    Ci = max(rho(span(0,i-1)));
    vCi.fill(Ci);
    
    // probability vector
    pvec = zeros<vec>(Ki+1);
    
    // fill in counts matrix and probability vector
    for (j=0;j<Ki;j++){
      wch = find(pi(span(0,i-1))==uniq(j));
      pvec(j) = similarityx(rho(i),rho(wch),Ci+1,tau,trans(X.row(i)),trans(X.rows(wch)),kappa0,1);
    }
    pvec(Ki) = similarityx(rho(i),vCi+1,Ci+1,tau,trans(X.row(i)),zeros<mat>(X.n_cols,X.n_cols),kappa0,0);
    pvecnorm = pvec/sum(pvec);
    lprob(i) = log(pvecnorm(pi(i)-1));
  }
  return exp(sum(lprob));
}

// evaluate PMF of LSP distribution starting at item m
// [[Rcpp::export]]
double dblocklsp(vec const& pi,vec const& rho,double const& tau, int start,int end){
  
  // initialize starting items
  int N = rho.n_elem;
  vec vCi = vec(1);
  vec lprob = zeros<vec>(N);
  if(start>1){
    start = start - 1;
  }
  
  for (int i=start;i<end;i++){
    
    // number of groups
    vec uniq = unique(pi(span(0,i-1)));
    int Ki = uniq.n_elem;
    
    // number of groups in rho
    int Ci = max(rho(span(0,i-1)));
    vCi.fill(Ci);
      
    // probability vector
    vec pvec = zeros<vec>(Ki+1);
    
    // fill in counts matrix and probability vector
    for (int j=0;j<Ki;j++){
      uvec wch = find(pi(span(0,i-1))==uniq(j));
      pvec(j) = similarity(rho(i),rho(wch),Ci+1,tau);
    }
    pvec(Ki) = similarity(rho(i),vCi+1,Ci+1,tau);
    vec pvecnorm = pvec/sum(pvec);
    lprob(i) = log(pvecnorm(pi(i)-1));
  }
  return exp(sum(lprob));
}


//-----------------------------------------#
// This MCMC sampling has been written by Sally Paganin
// https://github.com/adam-n-smith/demand-models-random-partitions
//-----------------------------------------#


// mcmc sampler for isolated demand model w/ random partition and LSP prior
// [[Rcpp::export]]
List isomvreglspmcmc(mat const& Q, mat const& P, List const& Zlist,
                     double const& abeta, double const& betab, mat const& Apsi, vec const& psibar, int const& nu, mat const& V, vec const& rho, double const& tau,
                     double const& pi_step,
                     int const& Rep, int const& keep, int const& nprint){
  
  // number of products, observations, parameters
  double n = Q.n_cols;
  double T = P.n_rows;
  
  // data pre-processing
  int nvar = 0;
  vec indreg = ones(n+1);
  vec nk = zeros(n);
  mat Zstar = Zlist[0];
  indreg(0) = nvar + 1;
  nvar = nvar + Zstar.n_cols;
  nk(0) = Zstar.n_cols;
  for(int i=1;i<n;i++){
    indreg(i) = nvar + 1;
    mat Z = Zlist[i];
    nvar = nvar + Z.n_cols;
    nk(i) = Z.n_cols;
    Zstar = join_rows(Zstar, Z);
  }
  indreg(n) = nvar+1;
  vec cumnk = cumsum(nk);
  mat ZspZs = trans(Zstar)*Zstar;
  
  // storage matrices
  mat betadraw(Rep/keep,n*n);
  mat pidraw(Rep/keep,n);
  mat psidraw(Rep/keep,nvar);
  mat Sigmadraw(Rep/keep,n*n);
  vec loglikedraw(Rep/keep);
  
  // initialize miscellaneous objects
  double logold, lognew, logold_prior, lognew_prior, logold_prop, lognew_prop, ldiff, u;
  int i, j, mkeep, rep;
  vec psitilde, yti, btilde, beta_c, pi_c;
  mat RR, Y, Zpsi, IR, ucholinv, EEVinv, Ztipyti, Ydti, Z;
  mat B_c, PB_c, Y_c, RR_c;
  List rwout;
  
  // initial values
  vec pi = ones(n);
  double K = max(pi);
  mat Pi = zeros(n,K);
  for(i=0;i<n;i++){
    Pi(i,pi(i)-1)=1;
  }
  mat Bmult = Pi*eye<mat>(K,K)*trans(Pi);
  double nbeta = as_scalar(sum(Pi,0)*trans(sum(Pi,0)));
  vec psi = zeros(nvar);
  mat B = zeros(n,n);
  vec beta = B(find(Bmult==1));
  mat PB = P*B;
  mat E = zeros(T,n);
  mat Sigma = eye(n,n);
  mat Sigmainv = inv(Sigma);
  
  // prior
  vec betabar = betab*ones(nbeta);
  mat Abeta = abeta*eye(nbeta,nbeta);
  mat rootAbeta = chol(Abeta);
  
  // accept/reject objects
  vec acpt = ones<vec>(2);
  double naccept = 0;
  
  // print time
  time_t itime;
  itime = time(NULL);
  if(nprint>0){
    Rcout << " MCMC Iteration (estimated time to end in minutes) \n";
  }
  
  // MCMC loop
  for (rep=0; rep<Rep; rep++){
    
    // -------------------------------------------------------------- //
    // 1. draw Sigma, psi
    // -------------------------------------------------------------- //
    
    Y = Q - PB;
    
    // draw psi | Sigma --------------------------------------------- //
    
    // compute Xtilde'Xtilde
    mat ZtipZti = zeros<mat>(sum(nk), sum(nk));
    for (i=0; i<n; i++){
      for (j=0; j<n; j++){
        ZtipZti(span(cumnk[i]-nk[i],cumnk[i]-1), span(cumnk[j]-nk[j],cumnk[j]-1)) =
          Sigmainv(i,j) * ZspZs(span(cumnk[i]-nk[i],cumnk[i]-1), span(cumnk[j]-nk[j],cumnk[j]-1));
      }
    }
    
    // now compute Xtilde'ytilde
    Ydti = Y*Sigmainv;
    Ztipyti = trans(as<mat>(Zlist[0]))*Ydti(span::all,0);
    for (i=1; i<n; i++){
      mat Z = Zlist[i];
      Ztipyti = join_cols(Ztipyti, trans(Z)*Ydti(span::all,i));
    }
    
    // draw psi
    IR = solve(trimatu(chol(ZtipZti + Apsi)), eye(nvar,nvar));
    psitilde = (IR*trans(IR))*(Ztipyti + Apsi*psibar);
    psi = psitilde + IR*as<arma::vec>(rnorm(nvar));
    
    // draw Sigma | psi --------------------------------------------- //
    
    Z = as<mat>(Zlist[0]);
    Zpsi = Z*psi(span(indreg[0]-1,indreg[1]-2));
    E(span::all,0) = Y.col(0) - Z*psi(span(indreg[0]-1,indreg[1]-2));
    for (i=1; i<n; i++){
      Z = as<mat>(Zlist[i]);
      Zpsi = join_cols(Zpsi,Z*psi(span(indreg[i]-1,indreg[i+1]-2)));
      E.col(i) = Y.col(i) - Z*psi(span(indreg[i]-1,indreg[i+1]-2));
    }
    Zpsi.set_size(T,n);
    
    // compute the inverse of E'E+V
    ucholinv = solve(trimatu(chol(trans(E)*E+V)), eye(n,n));
    EEVinv = ucholinv*trans(ucholinv);
    
    // draw of Sigma, Sigmainv
    rwout = rwishart(nu+T, EEVinv);
    Sigma = as<mat>(rwout["IW"]);
    Sigmainv = as<mat>(rwout["W"]);
    
    // -------------------------------------------------------------- //
    // 2. draw pi, beta
    // -------------------------------------------------------------- //
    
    // propose partition
    vec pi_c = vectorise(rlsp(1,pi,pi_step));
    double K_c = max(pi_c);
    mat Pi_c = zeros(n,K_c);
    for(int i=0;i<n;i++){
      Pi_c(i,pi_c(i)-1) = 1;
    }
    double nbeta_c = as_scalar(sum(Pi_c,0)*trans(sum(Pi_c,0)));
    mat Bmult_c = Pi_c*eye(K_c,K_c)*trans(Pi_c);
    
    // prior
    vec betabar_c = betab*ones(nbeta_c);
    mat Abeta_c = abeta*eye(nbeta_c,nbeta_c);
    mat rootAbeta_c = sqrt(abeta)*eye(nbeta_c,nbeta_c);
    
    // update old proposal
    List bout = drawbeta(Q,P,Zpsi,Sigmainv,n,pi,Pi,nbeta,Bmult,Abeta,betabar);
    mat root = bout["root"];
    mat IR = bout["IR"];
    vec btilde = bout["btilde"];
    
    // form new proposal
    List bout_c = drawbeta(Q,P,Zpsi,Sigmainv,n,pi_c,Pi_c,nbeta_c,Bmult_c,Abeta_c,betabar_c);
    mat root_c = bout_c["root"];
    mat IR_c = bout_c["IR"];
    vec btilde_c = bout_c["btilde"];
    
    // draw beta
    beta_c = btilde_c + IR_c*as<arma::vec>(rnorm(nbeta_c));
    
    // B matrix
    B_c = zeros(n,n);
    B_c(find(Bmult_c==1)) = beta_c;
    PB_c = P*B_c;
    
    // likelihood, prior, proposal
    logold = mvregll(Q,PB+Zpsi,Sigmainv);
    lognew = mvregll(Q,PB_c+Zpsi,Sigmainv);
    logold_prior = lndmvn(beta,betabar,rootAbeta) + log(dlsp(pi,rho,tau));
    lognew_prior = lndmvn(beta_c,betabar_c,rootAbeta_c) + log(dlsp(pi_c,rho,tau));
    logold_prop = lndmvn(beta,btilde,root) + log(dlsp(pi,pi_c,pi_step));
    lognew_prop = lndmvn(beta_c,btilde_c,root_c) + log(dlsp(pi_c,pi,pi_step));
    
    // accept/reject
    ldiff = lognew + lognew_prior + logold_prop - logold - logold_prior - lognew_prop;
    acpt(1) = exp(ldiff);
    u = as_scalar(randu(1));
    if (u < acpt.min()){
      pi = pi_c;
      Pi = Pi_c;
      K = K_c;
      nbeta = nbeta_c;
      Abeta = Abeta_c;
      rootAbeta = rootAbeta_c;
      betabar = betabar_c;
      Bmult = Bmult_c;
      beta = beta_c;
      B = B_c;
      PB = PB_c;
      btilde = btilde_c;
      IR = IR_c;
      logold = lognew;
      naccept = naccept + 1;
    }
    
    // -------------------------------------------------------------- //
    // print time and store draws
    // -------------------------------------------------------------- //
    
    if (nprint>0) if ((rep+1)%nprint==0){
      time_t ctime = time(NULL);
      char buf[32];
      double timetoend = difftime(ctime, itime) / 60.0 * (Rep - rep - 1) / (rep+1);
      sprintf(buf, " %d (%.1f)\n", rep+1, timetoend);
      Rcout <<  buf;
    }
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraw(mkeep-1, span::all) = trans(vectorise(B));
      pidraw(mkeep-1, span::all) = trans(pi);
      psidraw(mkeep-1, span::all) = trans(psi);
      Sigmadraw(mkeep-1, span::all) = vectorise(Sigma,1);
      loglikedraw(mkeep-1) = logold;
    }
    
  }
  
  // print total time elapsed
  if (nprint>0){
    time_t ctime = time(NULL);
    char buf[32];
    sprintf(buf, " Total Time Elapsed: %.2f \n", difftime(ctime, itime) / 60.0);
    Rcout << buf;
    itime = 0;
  }
  
  // output
  return List::create(
    Named("betadraw") = betadraw,
    Named("pidraw") = pidraw,
    Named("psidraw") = psidraw,
    Named("Sigmadraw") = Sigmadraw,
    Named("loglikedraw") = loglikedraw,
    Named("acptrate") = naccept/Rep
  );
  
}
