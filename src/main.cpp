#include "slau_gradient.hh"

const double SML = 1.0e-4;

int
main()
{
    dmat A = generateMatrix(1000, 1000);
    dvec ans = generateVector(1000, 1000);
    dvec B = matrix_vec(A, ans);
    dvec omp_solved = omp_solve(A, B);
#ifdef _OPENMP
    dvec solved = omp_solve(A, B);
#else
    dvec solved = solve(A, B);
#endif
}