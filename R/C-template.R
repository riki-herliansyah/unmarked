# The following functions are created via Rcpp (C++) package for more efficiency


library(Rcpp)

# It works in a similar way to rowSums available in basic R
src <- "
NumericVector rowsumC(NumericMatrix x) {
  int J = x.cols();
  int N = x.rows();
  NumericVector out(J);
  for(int i = 0; i < N; ++i) {
    double sum = 0;
    for (int j=0; j<J; j++){
      sum += x(i,j);
    }
    out[i] = sum;
  }
  return(out);
}"
cppFunction(src)

# It works in a similar way to colSums available in basic R
src <- "
NumericVector colsumC(NumericMatrix x) {
  int J = x.cols();
  int N = x.rows();
  NumericVector out(J);
  for(int j = 0; j < J; ++j) {
    double sum = 0;
    for (int i=0; i<N; i++){
      sum += x(i,j);
    }
    out[j] = sum;
  }
  return(out);
}"
cppFunction(src)

# The function is created to compute the Euclidean distance matrix given input S (activity centres)
# and X (the trap locations)
src <- "
NumericMatrix eudistC(NumericMatrix S, NumericMatrix X) {
  int N = S.rows();
  int J = X.rows();
  NumericMatrix d(N,J);
  for (int j=0; j<J; j++){
          d(_, j) = sqrt(pow(S(_,0)-X(j,0),2) + pow(S(_,1)-X(j,1),2));
        }
  return(d);
}"
cppFunction(src)

# The lamC function is used to compute the detection matrix, lambda_ij for given euclidean distance d,
# and the scale parameter sigma
src <- "
NumericMatrix lamC(NumericMatrix d, double sigma) {
  int N = d.rows();
  int J = d.cols();
  NumericMatrix lam(N,J);
  for (int j=0; j<J; j++){
      lam(_, j) = exp(-d(_,j)*d(_,j)/(2*sigma*sigma));
  }
  return(lam);
}"
cppFunction(src)

# It works in a similar way to dpois available in basic R
# I am not sure if it is needed but just to make it more efficient.
src <- "
NumericVector dpoisC(NumericVector y, NumericVector lambda, int give_log) {
  int J = y.size();
  NumericVector logans(J);
  for (int j=0; j<J; j++){
      logans[j] = y[j]*log(lambda[j]) - lambda[j] - lgamma(y[j]+1);
  }
  if(give_log) return logans; 
  else return exp(logans);
}"
cppFunction(src)

# This function is created to obtain the mode 
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
