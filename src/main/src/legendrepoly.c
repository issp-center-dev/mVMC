/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Legendre Polynomial
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

double LegendrePoly(const double x, const int n){
  double P01, P02, P03;
  int i;
  if(n<=0) return 1.0;
  else if(n==1) return x;
  else{
    P01 = 1.0;
    P02 = x;
    for(i=2;i<=n;i++){ /* recurrence relation */
      P03 = 1.0/(1.0*i) * ( (2.0*i-1.0)*x*P02 - (1.0*i-1.0)*P01 );
      P01 = P02;
      P02 = P03;
    }
    return P03;
  };
}
