#ifndef tridiag_H
#define tridiag_H

using namespace std;

void init_toeplitz(size_t n, double a, double b, double *alpha, double *beta)
{
    for(size_t i = 0; i<n; i++)
    {
        alpha[i] = a;
        if (i<n-1)
            beta[i] = b;
    }
}

void tridiag_sym(size_t n, double *alpha, double *beta, double *b)
{
    double temp = beta[0];
    
    for (size_t k = 1; k<n; k++)
    {
        temp = beta[k-1];
        beta[k-1] = temp/alpha[k-1];
        alpha[k] -= temp*beta[k-1];
    }
    
    for (size_t k = 1; k<n; k++)
    {
        b[k] -= beta[k-1]*b[k-1];
    }
    
    b[n-1] /= alpha[n-1];
    
    for (size_t k = n-2; k>0; k--)
    {
        b[k] = b[k]/alpha[k] - beta[k]*b[k+1];
    }
}

#endif //tridiag_H