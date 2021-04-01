// This function has been simplified and adapted from the sfreemap R package (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1554-7)

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(.dwelling_times)]]
NumericMatrix dwelling_times(List tree,
                             List Q_eigen,
                             NumericMatrix Q,
                             List fl
) {
  
  NumericVector qvec_vals = Q_eigen["values"];
  NumericMatrix qvec = Q_eigen["vectors"];
  NumericMatrix qvec_inv = Q_eigen["vectors_inv"];
  
  NumericVector edge_weights = tree["edge.length"];
  NumericVector edges = tree["edge"];
  int n_edges = edge_weights.size();
  
  NumericMatrix flF = fl["F"];
  NumericMatrix flG = fl["G"];
  NumericMatrix flS = fl["S"];
  double likelihood = fl["L"];
  
  int i, j, k, r, l, b, s, t;
  double diff, Iij;
  double TOL = 1e-8;
  int n_states = qvec.ncol();
  int n_comb = n_states*n_states;
  
  NumericMatrix multiplier(n_states, n_states);
  NumericMatrix results(n_edges, n_states);
  
  NumericMatrix S(n_states, n_comb);
  NumericMatrix Si(n_states, n_comb);
  NumericMatrix H(n_edges, n_states*n_states);
  NumericVector prm(n_edges);
  
  for (s=0; s<n_states; s++){
    
    if (s == 0) {
      t = 1;
    } else {
      t = 0;
    }
    
    multiplier(s,t) = Q(s,t);
    
    // calculate S
    
    for (i=0; i<n_states; i++) {
      for (j=0; j<n_states; j++){
        for (k=0; k<n_states; k++){
          S(i, (n_states*j)+k) = qvec(j,i)*qvec_inv(i,k);
        }
      }
    }
    
    // calculate Si
    for (i=0; i<n_states; i++) {
      for (j=0; j<n_states; j++) {
        for (k=0; k<n_states; k++) {
          Si(i, (n_states*j)+k) = S(i, (n_states*j)+0) * multiplier(0,k) + 
            S(i, (n_states*j)+1) * multiplier(1,k);
        }
      }
    }
    
    // calculate H
    
    std::fill( H.begin(), H.end(), 0 ) ;
    
    for (k=0; k<n_edges; k++){
      for (i=0; i<n_states; i++) {
        for (j=0; j<n_states; j++) {
          diff = qvec_vals(i) - qvec_vals(j);
          
          if (std::abs(diff) < TOL) {
            Iij = edge_weights(k) * (std::exp(qvec_vals(i)*edge_weights(k)));
          } else {
            Iij = (std::exp(qvec_vals(i)*edge_weights(k)) - std::exp(qvec_vals(j)*edge_weights(k))) / (diff);
          }
          H(k, 0) += Iij * (Si(i,0)*S(j,0) + Si(i,1)*S(j,2));
          H(k, 1) += Iij * (Si(i,0)*S(j,1) + Si(i,1)*S(j,3));
          H(k, 2) += Iij * (Si(i,2)*S(j,0) + Si(i,3)*S(j,2));
          H(k, 3) += Iij * (Si(i,2)*S(j,1) + Si(i,3)*S(j,3));
        }
      }
    }
    
    // calculate prm
    std::fill( prm.begin(), prm.end(), 0 ) ;
    
    for (k=0; k<n_edges; k++){
      l = edges(k,0) - 1;
      r = edges(k,1) - 1;
      
      if (k%2==0){
        b = edges(k+1,1)-1;
      } else {
        b = edges(k-1,1)-1;
      }
      
      for (i=0; i<n_states; i++) {
        for (j=0; j<n_states; j++) {
          prm(k) += flG(l,i) * flS(b,i) * flF(r,j) * H(k, j + i*2);
        }
      }
    }
    
    // calculate ev
    for (k=0; k<n_edges; k++){
      results(k, s) = prm(k) / likelihood;
    }
    multiplier(s,t) = 0;
  }
  
  return results;
}