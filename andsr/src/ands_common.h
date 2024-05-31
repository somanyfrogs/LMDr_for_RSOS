#ifndef ANDS_COMMON_H
#define ANDS_COMMON_H

#include <vector>
#include <type_traits>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <Eigen/Dense>
#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

// Type checker of std::vector.
template <typename T>
struct is_vec : std::false_type{};

template <typename T>
struct is_vec< std::vector<T> > : std::true_type{};

template <typename T>
bool is_vector(T x) {
    int res = is_vec<decltype(x)>::value;
    return res != 0;
}

// Implementation for index subsetting of matrix or vector in Eigen.
// these functions have already been implemented in Eigen 3.4.
// but RcppEigen has not yet adopted to Eigen 3.4.
// Thus, these will be removed with the future update of RcppEigen for Eigen 3.4.
template <typename VectorType>
VectorType vecSubset(const VectorType& V, const std::vector<size_t>& idx) {
    size_t m = idx.size();
    VectorType V_(m);

    for (size_t i = 0; i < m; ++i) V_(i) = V(idx[i]);
    return V_;
}

template <typename MatrixType>
MatrixType matSubset(const MatrixType& M, const std::vector<size_t>& idx, bool byrow = true) {
    size_t m = idx.size();
    MatrixType M_;

    if (byrow) {
        M_.resize(m, M.cols());
        for (size_t i = 0; i < m; ++i) M_.row(i) = M.row(idx[i]);
    } else {
        M_.resize(M.rows(), m);
        for (size_t i = 0; i < m; ++i) M_.col(i) = M.col(idx[i]);
    }

    return M_;
}

template <typename MatrixType>
MatrixType matSubset(const MatrixType& M, const std::vector<size_t>& idxRow, const std::vector<size_t>& idxCol) {
    size_t m = idxRow.size(), n = idxCol.size();
    MatrixType M_(m, n);

    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            M_(i, j) = M(idxRow[i], idxCol[j]);

    return M_;
}

// Define R 'seq' like function.
// from, to: the starting and end value of sequences.
// len: desired length of the sequence.
template<typename T>
std::vector<double> getSeq(T from, T to, size_t len) {
    std::vector<double> seq;
    double from_ = static_cast<double>(from), to_ = static_cast<double>(to), len_ = static_cast<double>(len);

    if (len == 0) { return seq; }

    seq.push_back(from_);
    if (len == 1) { return seq; }

    double delta = (to_ - from_) / (len_ - 1);
    for (size_t i = 1; i < len - 1; ++i) seq.push_back(from_ + delta * i);

    return seq;
}

// forward declarations of functions and classes.
std::vector<size_t> sortID(const std::vector<double>& vec, const std::vector<size_t>& idx);

double medianMPD(double beta, size_t h = 10);

Eigen::VectorXd getPred(const Eigen::MatrixXd& C, const Eigen::MatrixXd& X);

double getRMSE(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);

double getMAE(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);

double getRho(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);

Eigen::MatrixXd getNEDist(const Eigen::MatrixXd& X);

Eigen::MatrixXd getGeodesDist(const Eigen::MatrixXd& dmat, size_t k);

Eigen::MatrixXd getLMDist(const Eigen::MatrixXd& X, const Eigen::MatrixXd& theta);

// Function objects for increasing sorting.
class PrInc
{
    public:
        bool operator()(const std::pair<size_t, double> left, const std::pair<size_t, double> right){
            return left.second < right.second;
        }
};

// Class for boost RND with uniform_real.
class RngUnif
{
public:
    // Constructor
    RngUnif();
    RngUnif(int seed);
    RngUnif(double min, double max);
    RngUnif(int seed, double min, double max);

    double operator()();

private:
    typedef boost::uniform_real<> Dist;
    boost::variate_generator<boost::mt19937, Dist> rng;
};

// Class for boost RND with normal_distribution.
class RngNorm
{
public:
    // Constructor
    RngNorm();
    RngNorm(int seed);
    RngNorm(double mean, double sigma);
    RngNorm(int seed, double mean, double sigma);

    double operator()();

private:
    typedef boost::normal_distribution<> Dist;
    boost::variate_generator<boost::mt19937, Dist> rng;
};

#endif
