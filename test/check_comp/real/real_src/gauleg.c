/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Gauss-Legendre Quadrature
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#define GaussLeg_EPS 5.0e-14

/* calculate n points x[] and weight w[] for integration from x1 to x2. */
void GaussLeg(const double x1, const double x2, double *x, double *w, const int n){
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;

  m =(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  #pragma omp parallel for default(shared)        \
    private(i,j,z,p1,p2,p3,pp,z1)
  for(i=0;i<m;i++){
    z = cos( 3.14159265358979323846 * (i+0.75)/(n+0.5) );
    do{
      p1=1.0;
      p2=0.0;
      for(j=1;j<=n;j++){
        p3=p2;
        p2=p1;
        p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while(fabs(z-z1) > GaussLeg_EPS);
    x[i]=xm-xl*z;
    x[n-i-1]=xm+xl*z;
    w[i]=2.0*xl/((1.0-z*z)*pp*pp);
    w[n-i-1]=w[i];
  }

  return;
}
