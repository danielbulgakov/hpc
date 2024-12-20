#ifndef SRC_SLAU_GRADIENT_ORIG_HH_
#define SRC_SLAU_GRADIENT_ORIG_HH_

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>

using dvec = std::vector<double>;
using dmat = std::vector<dvec>;

/**
 * @brief generate random simmetric positive defined (spd) matrix
 * 
 * @param size size of matrix
 * @param seed seed
 * @return matrix (dmat, std::vector<std::vector<double>>)
 */
dmat generateMatrix(int size, unsigned int seed);

/**
 * @brief generate random vector
 * 
 * @param size vector size
 * @param seed seed
 * @return vector
 */
dvec generateVector(int size, unsigned int seed);
/**
 * @brief scalar multiplication of vectors
 * 
 * @param a first vector
 * @param b second vector
 * @return result value 
 */
double vec_vec(const dvec &a, const dvec &b);

/**
 * @brief matrix-vector multiplication
 * 
 * @param a matrix (dmat, std::vector<std::vector<double>>)
 * @param b vector
 * @return result vector
 */
dvec matrix_vec(const dmat &a, const dvec &b);

/**
 * @brief linear combination of vectors
 * 
 * @param ma scalar value mult for first vector
 * @param a first vector
 * @param mb scalar value mult for second vector
 * @param b second vector
 * @return result vector
 */
dvec vec_vec_comb(double ma, const dvec& a, double mb, const dvec& b);

/**
 * @brief vector normalizing
 * 
 * @param a vector
 * @return magnitude 
 */
double vec_norm(const dvec& a);

/**
 * @brief conjagurate gradient method solver
 * 
 * @param a matrix (dmat, std::vector<std::vector<double>>)
 * @param b vector
 * @return result 
 */
dvec solve(const dmat &a, const dvec& b);

#endif  // SRC_SLAU_GRADIENT_ORIG_HH_