
#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>

/* Parameters */
#define N 5
#define NRHS 3
#define LDA N
#define LDB NRHS

/* Main program */
int main() {
        /* Locals */
        lapack_int n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
        /* Local arrays */
        lapack_int ipiv[N];
        double a[LDA*N] = {
            6.80, -6.05, -0.45,  8.32, -9.67,
           -2.11, -3.30,  2.58,  2.71, -5.14,
            5.66, 5.36, -2.70,  4.35, -7.26,
            5.97, -4.44,  0.27, -7.17, 6.08,
            8.23, 1.08,  9.04,  2.14, -6.87
        };
        double b[LDB*N] = {
            4.02, -1.56, 9.81,
            6.19,  4.00, -4.09,
           -8.22, -8.67, -4.57,
           -7.57,  1.75, -8.61,
           -3.03,  2.86, 8.99
        };

       double aNorm;
       double rcond;
       char ONE_NORM = '1';
       lapack_int NROWS = n;
       lapack_int NCOLS = n;
       lapack_int LEADING_DIMENSION_A = n;

        /* Executable statements */
        printf( "LAPACKE_dgecon Example Program Results\n" );
        aNorm = LAPACKE_dlange(LAPACK_ROW_MAJOR, ONE_NORM, NROWS, NCOLS, a, LEADING_DIMENSION_A);
        info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, NROWS, NCOLS, a, LEADING_DIMENSION_A, ipiv);
        info = LAPACKE_dgecon(LAPACK_ROW_MAJOR, ONE_NORM, n, a, LEADING_DIMENSION_A, aNorm, &rcond); // aNorm should be 35.019999999999996
        double work[4*N];
        int iwork[N];
        //info = LAPACKE_dgecon_work(LAPACK_ROW_MAJOR, ONE_NORM, n, a, LEADING_DIMENSION_A, aNorm, &rcond, work, iwork); // aNorm should be 35.019999999999996
        //dgecon_( &ONE_NORM, &n, a, &LEADING_DIMENSION_A, &aNorm, &rcond, work, iwork, &info );
        /* Check for the exact singularity */
              if (info == 0)
              {
                     printf("LAPACKE_dgecon completed SUCCESSFULLY...\n");
              }
              else if ( info < 0 )
              {
            printf( "Element %d of A had an illegal value\n", -info );
            exit( 1 );
        }
              else
              {
            printf( "Unrecognized value of INFO = %d\n", info );
            exit( 1 );
              }

        /* Print solution */
       printf("LAPACKE_dlange / One-norm of A = %lf\n", aNorm);
        printf("LAPACKE_dgecon / RCOND of A    = %f\n", rcond);
        exit( 0 );
} /* End of LAPACKE_dgesv Example */