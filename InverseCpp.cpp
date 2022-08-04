#include <RcppArmadillo.h>
#include <string>
#include <unordered_map>
#include <vector>

// [[Rcpp::plugins(cpp17)]] 

using namespace arma;
using namespace Rcpp;
using namespace std;

class HashCSRMatrix{
  unordered_map<u64, u32> M;
  unordered_map<u64, u32>::iterator it;
  vector<double> values;
  u32 pos = 0;
  int nnz = 0;
  
public: 
  void insert(const u32 & i, const u32 & j, const double & val){
    u64 ij;
    if(i >= j){
      ij = ((u64)i) << 32 | j;
    } else if (i < j){
      ij = ((u64)j) << 32 | i;
    }
    it     = M.find(ij);
    if(it == M.end()){
      if(i >= j){
        pos   = nnz;
        M[ij] = pos;
        values.push_back(val);
        nnz++;
      }
    } else if (it != M.end()){
      if(i >= j){
        pos = M[ij];
        values[pos] += val;
      }
    }
  }
  double getval(const u32 & i, const u32 & j) {
    u64    ij;
    double val = 0.0;
    if(i >= j){
      ij = ((u64)i) << 32 | j;
    } else if (i < j){
      ij = ((u64)j) << 32 | i;
    }
    it     = M.find(ij);
    if(it != M.end()){
      val = values[it->second];
    }
    return(val);
  }
  int nelem() const { 
    return M.size(); 
  }
  void printvals(){
    for (auto & it: M) {
      u32    i   = (u32)(it.first >> 32); 
      u32    j   = (u32)it.first;
      double val = values[it.second];
      Rcpp::Rcout << i << " " << j << " " << val << "\n";
    }
  }
};

// [[Rcpp::depends(RcppArmadillo)]]

/*vec Inbreeding(const ivec & Sire, const ivec & Dam){
  int nanim = Sire.n_rows;
  int S     = 0;
  int D     = 0;
  int t     = 0;
  int i     = 0;
  int j     = 0;
  
  ivec Anc = zeros<ivec>(nanim + 1);  
  ivec Lap = zeros<ivec>(nanim + 1); 
  vec Diag = zeros<vec>(nanim + 1);
  vec F    = zeros<vec>(nanim + 1);
  vec B    = zeros<vec>(nanim + 1);
  vec L    = zeros<vec>(nanim + 1);
  
  F(0)   =-1; 
  Lap(0) =-1;
  for(i = 1, t = -1; i <= nanim; i++) { 
    S = Sire(i - 1); 
    D = Dam(i - 1);
    Lap(i) = ((Lap(S) < Lap(D)) ? Lap(D) : Lap(S)) + 1;
    if (Lap(i) > t) t = Lap(i);
  }
  ivec St = zeros<ivec>(t + 1);
  ivec Mi = zeros<ivec>(t + 1);
  
  for(i = 1; i <= nanim; i++) {
    S    = Sire(i - 1); 
    D    = Dam(i - 1);
    B(i) = 0.5 - 0.25*(F(S) + F(D)); 
    for (j = 0; j < Lap(i); j++) {
      ++St(j); 
      ++Mi(j);
    } 
    if (S == 0 || D == 0) {
      F(i) = L(i) = 0; 
      continue;
    }
    if(S == Sire(i - 2) && D == Dam(i - 2)) {
      F(i) = F(i - 1); 
      L(i) = L(i - 1); 
      continue;
    }
    F(i)         = -1; 
    L(i)         = 1; 
    t            = Lap(i); 
    Anc(Mi(t)++) = i; 
    while(t > -1) {
      j = Anc(--Mi(t)); 
      S = Sire(j - 1); 
      D = Dam(j - 1); 
      if (S) {
        if (!L(S)){
          Anc(Mi(Lap(S))++) = S;
        }  
        L(S) += 0.5*L(j); 
      }
      if (D) {
        if (!L(D)){
          Anc(Mi(Lap(D))++) = D;
        } 
        L(D) += 0.5*L(j);
      }
      F(i) += L(j)*L(j)*B(j);
      L(j) = 0;
      if (Mi(t) == St(t)) --t;
    } 
  }
  return(F.tail(nanim));
}*/

vec Dvec(const ivec & Sire, const ivec & Dam){
  int nanim = Sire.n_rows;
  int S     = 0;
  int D     = 0;
  int t     = 0;
  int i     = 0;
  int j     = 0;
  
  ivec Anc = zeros<ivec>(nanim + 1);  
  ivec Lap = zeros<ivec>(nanim + 1); 
  vec Diag = zeros<vec>(nanim + 1);
  vec F    = zeros<vec>(nanim + 1);
  vec B    = zeros<vec>(nanim + 1);
  vec L    = zeros<vec>(nanim + 1);
  
  F(0)   =-1; 
  Lap(0) =-1;
  for(i = 1, t = -1; i <= nanim; i++) { 
    S = Sire(i - 1); 
    D = Dam(i - 1);
    Lap(i) = ((Lap(S) < Lap(D)) ? Lap(D) : Lap(S)) + 1;
    if (Lap(i) > t) t = Lap(i);
  }
  ivec St = zeros<ivec>(t + 1);
  ivec Mi = zeros<ivec>(t + 1);
  
  for(i = 1; i <= nanim; i++) {
    S    = Sire(i - 1); 
    D    = Dam(i - 1);
    B(i) = 0.5 - 0.25*(F(S) + F(D)); 
    for (j = 0; j < Lap(i); j++) {
      ++St(j); 
      ++Mi(j);
    } 
    if (S == 0 || D == 0) {
      F(i) = L(i) = 0; 
      continue;
    }
    if(S == Sire(i - 2) && D == Dam(i - 2)) {
      F(i) = F(i - 1); 
      L(i) = L(i - 1); 
      continue;
    }
    F(i)         = -1; 
    L(i)         = 1; 
    t            = Lap(i); 
    Anc(Mi(t)++) = i; 
    while(t > -1) {
      j = Anc(--Mi(t)); 
      S = Sire(j - 1); 
      D = Dam(j - 1); 
      if (S) {
        if (!L(S)){
          Anc(Mi(Lap(S))++) = S;
        }  
        L(S) += 0.5*L(j); 
      }
      if (D) {
        if (!L(D)){
          Anc(Mi(Lap(D))++) = D;
        } 
        L(D) += 0.5*L(j);
      }
      F(i) += L(j)*L(j)*B(j);
      L(j) = 0;
      if (Mi(t) == St(t)) --t;
    } 
  }
  for(int i = 1; i <= nanim; i++){
    Diag(i) = B(i);
  }
  return(Diag.tail(nanim));
}

void A_times_v(vec & w, const ivec & v, const ivec & Sire, const ivec & Dam, const vec & F){
  int      n    = w.n_rows;
  int      s    = 0;
  int      d    = 0;
  double   dii  = 0.0;
  double   tmp  = 0.0;
  double   Fs   = 0;
  double   Fd   = 0;
  vec      q    = zeros<vec>(n);
  
  for(int i = n - 1; i >= 0; i--){
    q(i) = q(i) + v(i);
    s    = Sire(i);
    d    = Dam(i);
    if(s != 0){
      q(s - 1) += q(i)*0.5;
    }
    if(d != 0){
      q(d - 1) += q(i)*0.5;
    }
  }
  for(int i = 0; i < n; i++){
    s    = Sire(i);
    d    = Dam(i);
    if(s == 0){
      Fs = 0;
    } else if(s > 0){
      Fs = F(s - 1);
    }
    if(d == 0){
      Fd = 0;
    } else if(d > 0){
      Fd = F(d - 1);
    }
    dii  = ((s == 0)*1.0 + (d == 0)*1.0 + 2.0)/4.0 - 0.25*(Fs + Fd);
    tmp  = 0.0;
    if(s != 0){
      tmp += w(s - 1);
    }
    if(d != 0){
      tmp += w(d - 1);
    }
    w(i) = 0.5*tmp;
    w(i) += dii*q(i);
  }
}

mat Get_A22(const ivec & Sire, const ivec & Dam, const ivec & NodeID, const vec & F){
  int  n    = Sire.n_rows;
  int  nInd = NodeID.n_rows; 
  ivec v    = zeros<ivec>(n);
  vec  w    = zeros<vec>(n);
  mat  A22  = zeros<mat>(nInd, nInd);
  
  for(int i = 0; i < nInd; i++){
    v = zeros<ivec>(n);
    v(NodeID(i) - 1) = 1;
    A_times_v(w, v, Sire, Dam, F);
    for(int j = 0; j < nInd; j++){
      A22(j, i) = w(NodeID(j) - 1);
    }
  }
  return(A22);
}

/*
int TestingMat(const ivec & Row, const ivec & Col, const vec & Val){
  HashCSRMatrix X;
  int           n   = Row.n_rows;
  double        val = 0.0;
  u32           tt  = 0;
  
  for(int i = 0; i < n; i++){
    u32 j = Row(i);
    u32 k = Col(i);
    val   = Val(i);
    X.insert(j, k, val);
  }
  Rcpp::Rcout << "Element 0,0: " << X.getval(tt, tt) << "\n";
  return(0);
}*/

vec GetL(const mat & U){
  vec    L   = zeros<vec>(2);
  double det = U(1, 1)*U(2, 2) - U(2, 1)*U(1, 2);
  L(0)       = (U(2, 2)*U(1, 0) - U(2, 0)*U(1, 2))/det;
  L(1)       = (U(1, 1)*U(2, 0) - U(1, 0)*U(1, 2))/det;    
  return(L);
}

double GetCVar(const mat & U){
  double detu = 1.0/(U(1, 1)*U(2, 2) - U(1, 2)*U(1, 2));
  double var  = U(0, 0) - U(0, 1)*(U(2, 2)*U(1, 0) - U(2, 0)*U(1, 2))*detu - 
    U(0, 2)*(U(1, 1)*U(2, 0) - U(1, 0)*U(1, 2))*detu;
  return(var);
}

// [[Rcpp::export]]
sp_mat MakeKInv(const ivec & ID, const ivec & Sire, const ivec & Dam, const ivec & GenID,
                const mat Z, const int & scale){
  /*
   * Integer variables
   */
  int    nGeno   = GenID.n_rows;
  int    nAnim   = ID.n_rows;
  int    m       = Z.n_cols;
  int    m2      = Z.n_rows;
  int    iid     = 0;
  int    sid     = 0;
  int    did     = 0;
  int    gii     = 0;
  int    gis     = 0;
  int    gid     = 0;
  /*
   * Double precision variables
   */
  double val     = 0.0;
  double kk      = scale;
  double tt      = 0;
  /*
   * Others
   */
  wall_clock timer;
  /*
   * Vectors (integer, doubles) and matrices (doubles, and sparse)
   */
  ivec   pednode = zeros<ivec>(3);
  vec    L       = {1.0, -0.5, -0.5};
  vec    D       = Dvec(Sire, Dam);
  vec    Dg      = zeros<vec>(nGeno);
  mat    U       = zeros<mat>(3, 3);
  mat    Zi      = zeros<mat>(3, m);
  sp_mat Kinv    = zeros<sp_mat>(nAnim, nAnim);
  if(nGeno != m2){
    Rcpp::Rcout << "Number of genotype individuals not equal to dimensions of Z.\n";
  }
  
  /*
   * Set up A Inverse
   */
  timer.tic();
  for(int i = 0; i < nAnim; i++){
    pednode = {ID(i), Sire(i), Dam(i)};
    for(int j = 0; j < 3; j++){
      for(int k = 0; k < 3; k++){
        if(pednode(j) != 0 && pednode(k) != 0){
          int ii = pednode(j) - 1;
          int jj = pednode(k) - 1;
          val    = L(j)*L(k)/D(i);
          Kinv(ii, jj) += val;
        }
      }
    }
  }
  tt = timer.toc();
  Rcpp::Rcout << "A inverse completed in:      " << tt << " seconds\n";
  Rcpp::Rcout << "Number of animals:           " << nAnim << "\n";
  /*
   * Recoding Genotyped IDs to find rows of Z
   */
  std::unordered_map<int, int> NewID;
  std::unordered_map<int, int>::iterator it;
  for(int i = 0; i < nGeno; i++){
    it = NewID.find(GenID(i));
    if(it == NewID.end()){
      NewID.insert({GenID(i), i});
    }
  }
  /*
   * Adding G inverse 
   */
  timer.tic();
  for(int i = 0; i < nGeno; i++){
    Zi  = zeros<mat>(3, m);
    iid = GenID(i); 
    it  = NewID.find(iid);
    if(it != NewID.end()){
      gii       = it->second;
      Zi.row(0) = Z.row(gii);
    }
    sid = Sire(GenID(i) - 1);
    it  = NewID.find(sid);
    if(it != NewID.end()){
      gis       = it->second;
      Zi.row(1) = Z.row(gis);
    }
    did = Dam(GenID(i) - 1);
    it  = NewID.find(did);
    if(it != NewID.end()){
      gid       = it->second;
      Zi.row(2) = Z.row(gid);
    }
    mat XpX = Zi*Zi.t()/kk;
    it      = NewID.find(iid);
    if(it == NewID.end()){
      XpX(0, 0) = 1.0;
      gii       = 0;
    }
    it      = NewID.find(sid);
    if(it == NewID.end()){
      XpX(1, 1) = 1.0;
      gis       = 0;
    }
    it      = NewID.find(did);
    if(it == NewID.end()){
      XpX(2, 2) = 1.0;
      gid       = 0;
    }
    L(span(1, 2)) = -GetL(XpX);
    Dg(i)         = GetCVar(XpX);
    pednode       = {GenID(i), Sire(GenID(i) - 1), Dam(GenID(i) - 1)};
    for(int j = 0; j < 3; j++){
      for(int k = 0; k < 3; k++){
        if((pednode(j) != 0) && (pednode(k) != 0)){
          int ii = pednode(j) - 1;
          int jj = pednode(k) - 1;
          val    = L(j)*L(k)/Dg(i);
          Kinv(ii, jj) += val; 
        }
      }
    }
  }
  tt = timer.toc();
  Rcpp::Rcout << "G inverse completed in:      " << tt << " seconds\n";
  Rcpp::Rcout << "Number of genotyped animals: " << nGeno << "\n";
  Rcpp::Rcout << "Number of markers used     : " << m << "\n";
  return(Kinv);
}