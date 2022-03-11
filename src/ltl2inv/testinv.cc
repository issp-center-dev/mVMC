#include <cstdio>
#include <cstdlib>

#include "colmaj.hh"
#include "invert.tcc"

template <typename T>
double test_ltl2inv(int n)
{
    colmaj<T> A(new T[n * n], n);
    colmaj<T> M(new T[n * n], n);
    colmaj<T> X(new T[n * n], n);
    int *ipiv = new int[n];
    T *vT = new T[n-1];

    for (int j = 0; j < n; ++j) {
        A(j, j) = 0.0;
        X(j, j) = 0.0;

        for (int i = j+1; i < n; ++i) {
            A(i, j) = rand();
            A(i, j) /= RAND_MAX;

            A(j, i) = -A(i, j);
            X(i, j) = A(i, j);
        }
    }

    sktrf<T>(BLIS_LOWER, n, &X(0, 0), n, ipiv, &M(0, 0), n*n);
    ltl2inv<T>(n, &X(0, 0), n, ipiv, vT, &M(0, 0), n);

    gemm<T>(BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE,
            n, n, n,
            1.0,
            &A(0, 0), n,
            &X(0, 0), n,
            0.0,
            &M(0, 0), n);

#ifdef VERBOSE
    for (int i = 0; i < n; ++i)
        printf("%d ", ipiv[i]);
    printf("\n");

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

    return err;
}

int main(void)
{
    printf("double 10x10 err=%e.\n", test_ltl2inv<double>(10));
    printf("float 10x10 err=%e.\n", test_ltl2inv<float>(10));
    printf("dcomplex 10x10 err=%e.\n", test_ltl2inv<ccdcmplx>(10));
    printf("double 200x200 err=%e.\n", test_ltl2inv<double>(200));
    printf("float 200x200 err=%e.\n", test_ltl2inv<float>(200));
    printf("dcomplex 200x200 err=%e.\n", test_ltl2inv<ccdcmplx>(200));

    return 0;
}

