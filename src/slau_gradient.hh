// Copyright 2024 Bulgakov Daniil


#ifndef SRC_SLAU_GRADIENT_HH_
#define SRC_SLAU_GRADIENT_HH_


#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

const double SMOL = 1.0e-10;

using dvec = std::vector<double>;
using dmat = std::vector<dvec>;

//----------------------------- Test Data Generation --------------------------

/**
 * @brief generate random simmetric positive defined (spd) matrix
 *
 * @param size size of matrix
 * @param seed seed
 * @return matrix (dmat, std::vector<std::vector<double>>)
 */
dmat
generateMatrix(int size, unsigned int seed);

/**
 * @brief generate random vector
 *
 * @param size vector size
 * @param seed seed
 * @return vector
 */
dvec
generateVector(int size, unsigned int seed);

//----------------------------- Sequential Computing Methods ------------------

/**
 * @brief scalar multiplication of vectors
 *
 * @param a first vector
 * @param b second vector
 * @return result value
 */
double
vec_vec(const dvec &a, const dvec &b);

/**
 * @brief matrix-vector multiplication
 *
 * @param a matrix (dmat, std::vector<std::vector<double>>)
 * @param b vector
 * @return result vector
 */
dvec
matrix_vec(const dmat &a, const dvec &b);

/**
 * @brief linear combination of vectors
 *
 * @param ma scalar value mult for first vector
 * @param a first vector
 * @param mb scalar value mult for second vector
 * @param b second vector
 * @return result vector
 */
dvec
vec_vec_comb(double ma, const dvec &a, double mb, const dvec &b);

/**
 * @brief vector normalizing
 *
 * @param a vector
 * @return magnitude
 */
double
vec_norm(const dvec &a);

//----------------------------- Parallel Computing Methods --------------------

/**
 * @brief parallel matrix-vector multiplication
 *
 * @param a matrix (dmat, std::vector<std::vector<double>>)
 * @param b vector
 * @return result vector
 */
dvec
omp_matrix_vec(const dmat &a, const dvec &b);

/**
 * @brief parallel scalar multiplication of vectors
 *
 * @param a first vector
 * @param b second vector
 * @return result value
 */
double
omp_vec_vec(const dvec &a, const dvec &b);

//----------------------------- Conjugate Method Algorithm --------------------

/**
 * @brief conjugate gradient method solver
 *
 * @param a matrix (dmat, std::vector<std::vector<double>>)
 * @param b vector
 * @return result
 */
dvec
solve(const dmat &a, const dvec &b);

/**
 * @brief parallel conjugate gradient method solver
 *
 * @param a matrix (dmat, std::vector<std::vector<double>>)
 * @param b vector
 * @return result
 */
dvec
omp_solve(const dmat &a, const dvec &b);

#endif  // SRC_SLAU_GRADIENT_HH_