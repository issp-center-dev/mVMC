/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Quantum Projection
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

void InitQPWeight();
double complex CalculateLogIP_fcmp(double complex * const pfM, const int qpStart, const int qpEnd, MPI_Comm comm);
double complex CalculateIP_fcmp(double complex * const pfM, const int qpStart, const int qpEnd, MPI_Comm comm);
void UpdateQPWeight();

/* calculate SPGLCos, SPGLSin and QPFullWeight */
void InitQPWeight() {
  int i,j,idx;
  double *beta;
  double *weight;
  double w;

  RequestWorkSpaceDouble(2*NSPGaussLeg);
  beta   = GetWorkSpaceDouble(NSPGaussLeg);
  weight = GetWorkSpaceDouble(NSPGaussLeg);

  if(NSPGaussLeg==1) {
    beta[0] = 0.0;
    weight[0] = 1.0;
    SPGLCos[0] = 1.0;
    SPGLSin[0] = 0.0;
    SPGLCosSin[0] = 0.0;
    SPGLCosCos[0] = 1.0;
    SPGLSinSin[0] = 0.0;

    for(j=0;j<NMPTrans;j++) {
      QPFixWeight[j] = ParaQPTrans[j];
    }

  } else {
    GaussLeg(0, M_PI, beta, weight, NSPGaussLeg);

    #pragma omp parallel for default(shared) private(i,j,w,idx)
    for(i=0;i<NSPGaussLeg;i++) {
      SPGLCos[i]    = cos(0.5 * beta[i])+0*I;    //TBC
      SPGLSin[i]    = sin(0.5 * beta[i])+0*I;    //TBC
      SPGLCosSin[i] = SPGLCos[i]*SPGLSin[i]+0*I; //TBC
      SPGLCosCos[i] = SPGLCos[i]*SPGLCos[i]+0*I; //TBC
      SPGLSinSin[i] = SPGLSin[i]*SPGLSin[i]+0*I; //TBC
      
      w = 0.5*sin(beta[i])*weight[i]*LegendrePoly(cos(beta[i]), NSPStot);
      
      for(j=0;j<NMPTrans;j++) {
        idx = i + j*NSPGaussLeg;
        QPFixWeight[idx] = w * ParaQPTrans[j];
      }
    }
  }

  UpdateQPWeight();

  ReleaseWorkSpaceDouble();
  return;
}


/* Calculate logarithm of inner product <phi|L|x> */
double complex CalculateLogIP_fcmp(double complex * const pfM, const int qpStart, const int qpEnd, MPI_Comm comm) {
  const int qpNum = qpEnd-qpStart;
  double complex ip=0.0+0.0*I, ip2;
  int qpidx;
  int size;
  MPI_Comm_size(comm,&size);

  #pragma loop noalias
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    ip += QPFullWeight[qpidx+qpStart] * pfM[qpidx];
  }
  if(size>1) {
    MPI_Allreduce(&ip, &ip2, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);//TBC
    ip = ip2;
  }
  return clog(ip);
}

/* Calculate inner product <phi|L|x> */
double complex CalculateIP_fcmp(double complex * const pfM, const int qpStart, const int qpEnd, MPI_Comm comm) {
  const int qpNum = qpEnd-qpStart;
  double complex ip=0.0+0.0*I, ip2;
  int qpidx;
  int size;
  MPI_Comm_size(comm,&size);

  #pragma loop noalias
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    ip += QPFullWeight[qpidx+qpStart] * pfM[qpidx];
  }
  if(size>1) {
    MPI_Allreduce(&ip, &ip2, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
    ip = ip2;
  }
  return ip;
}

void UpdateQPWeight() {
  int i,j,offset;
  double complex tmp; //TBC

  if(FlagOptTrans>0) {
    for(i=0;i<NOptTrans;i++) {
      offset = i*NQPFix;
      tmp = OptTrans[i];
      for(j=0;j<NQPFix;j++) {
        QPFullWeight[offset+j] = tmp * QPFixWeight[j];
      }
    }
  } else {
    for(j=0;j<NQPFix;j++) {
      QPFullWeight[j] = QPFixWeight[j];
    }
  }

  return;
}
