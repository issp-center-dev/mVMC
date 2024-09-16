#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <chrono>

#include "colmaj.hh"
#include "invert.tcc"

struct info_t {
    double err;
    long long tim;
};

template <typename T>
info_t test_decomp2inv(uplo_t uplo, int n)
{
    using namespace std::chrono;

    colmaj<T> A(new T[n * n], n);
    colmaj<T> M(new T[n * n], n);
    colmaj<T> X(new T[n * n], n);
    int *ipiv = new int[n];
    T *vT = new T[n-1];

    assert(&A(0, 0));
    assert(&M(0, 0));
    assert(&X(0, 0));
    assert(ipiv);
    assert(vT);

    for (int j = 0; j < n; ++j) {
        A(j, j) = 0.0;
        X(j, j) = 0.0;

        for (int i = j+1; i < n; ++i) {
            A(i, j) = rand();
            A(i, j) /= RAND_MAX;

            A(j, i) = -A(i, j);
            switch (uplo) {
            case BLIS_LOWER:
                X(i, j) = A(i, j); break;
            case BLIS_UPPER:
                X(j, i) = A(j, i); break;
            default:
                break;
            }
        }
    }

    auto start = high_resolution_clock::now();

    sktrf<T>(uplo, n, &X(0, 0), n, ipiv, &M(0, 0), n*n);
    switch (uplo) {
    case BLIS_LOWER:
        ltl2inv<T>(n, &X(0, 0), n, ipiv, vT, &M(0, 0), n); break;
    case BLIS_UPPER:
        utu2inv<T>(n, &X(0, 0), n, ipiv, vT, &M(0, 0), n); break;
    default:
        break;
    }

    auto tim = duration_cast<nanoseconds>(high_resolution_clock::now() - start).count();

    gemm<T>(BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE,
            n, n, n,
            1.0,
            &A(0, 0), n,
            &X(0, 0), n,
            0.0,
            &M(0, 0), n);

#ifdef VERBOSE
    printf("iPiv = [ ");
    for (int i = 0; i < n; ++i)
        printf("%d ", ipiv[i]);
    printf("]\n");
    printf("vT = [ ");
    for (int i = 0; i < n-1; ++i)
        printf("%f ", vT[i]);
    printf("]\n");

    printf("A*inv(A) =\n");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            printf("%10.6f ", M(i, j));
        printf("\n");
    }
#endif

    double err = 0.0;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            err += std::abs(M(i, j) - double(!(i - j)));
    err /= n; // sqrt(n*n);

    delete[] &A(0, 0);
    delete[] &M(0, 0);
    delete[] &X(0, 0);
    delete[] ipiv;
    delete[] vT;

    return info_t{ err, tim };
}

int main(void)
{
    info_t info; int n;
    n= 10; info = test_decomp2inv<double>  (BLIS_LOWER, n); printf("double   n=%5d lower err=%e at %6lfgflops.\n", n, info.err, 1.33*std::pow(10 ,3)/info.tim);
    n= 10; info = test_decomp2inv<float>   (BLIS_LOWER, n); printf("float    n=%5d lower err=%e at %6lfgflops.\n", n, info.err, 1.33*std::pow(10 ,3)/info.tim);
    n= 10; info = test_decomp2inv<ccdcmplx>(BLIS_LOWER, n); printf("dcomplex n=%5d lower err=%e at %6lfgflops.\n", n, info.err, 5.33*std::pow(10 ,3)/info.tim);
    n=200; info = test_decomp2inv<double>  (BLIS_LOWER, n); printf("double   n=%5d lower err=%e at %6lfgflops.\n", n, info.err, 1.33*std::pow(200,3)/info.tim);
    n=200; info = test_decomp2inv<float>   (BLIS_LOWER, n); printf("float    n=%5d lower err=%e at %6lfgflops.\n", n, info.err, 1.33*std::pow(200,3)/info.tim);
    n=200; info = test_decomp2inv<ccdcmplx>(BLIS_LOWER, n); printf("dcomplex n=%5d lower err=%e at %6lfgflops.\n", n, info.err, 5.33*std::pow(200,3)/info.tim);

    n= 10; info = test_decomp2inv<double>  (BLIS_UPPER, n); printf("double   n=%5d upper err=%e at %6lfgflops.\n", n, info.err, 1.33*std::pow(10 ,3)/info.tim);
    n= 10; info = test_decomp2inv<float>   (BLIS_UPPER, n); printf("float    n=%5d upper err=%e at %6lfgflops.\n", n, info.err, 1.33*std::pow(10 ,3)/info.tim);
    n= 10; info = test_decomp2inv<ccdcmplx>(BLIS_UPPER, n); printf("dcomplex n=%5d upper err=%e at %6lfgflops.\n", n, info.err, 5.33*std::pow(10 ,3)/info.tim);
    n=200; info = test_decomp2inv<double>  (BLIS_UPPER, n); printf("double   n=%5d upper err=%e at %6lfgflops.\n", n, info.err, 1.33*std::pow(200,3)/info.tim);
    n=200; info = test_decomp2inv<float>   (BLIS_UPPER, n); printf("float    n=%5d upper err=%e at %6lfgflops.\n", n, info.err, 1.33*std::pow(200,3)/info.tim);
    n=200; info = test_decomp2inv<ccdcmplx>(BLIS_UPPER, n); printf("dcomplex n=%5d upper err=%e at %6lfgflops.\n", n, info.err, 5.33*std::pow(200,3)/info.tim);

    return 0;
}

