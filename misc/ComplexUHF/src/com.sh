#icc  UHFmain.c  -O3 -static  -L/opt/intel/mkl/10.0.3.020/lib/em64t -lmkl_lapack -lmkl_em64t -lguide -lpthread -o UHF
#icc UHFmain.c  -O3 -static  -L/opt/intel/mkl/8.1/lib/em64t -lmkl_lapack -lmkl_em64t -lguide -lpthread -o UHF
#icc  UHFmain.c  -openmp -openmp-report2 -O3 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread -o UHF
icc UHFmain.c  -O3 -LKLPATH -IKLINCLUDE -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o UHF
