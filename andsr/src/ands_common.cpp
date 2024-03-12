#include "ands_common.h"

using namespace std;
using namespace boost;
using namespace Eigen;
using namespace Rcpp;
using namespace RcppEigen;

// Sort index by ascending order of vec's values.
vector<size_t> sortID(const vector<double>& vec, const vector<size_t>& idx) {
    vector<size_t> tmp = idx;
    sort(tmp.begin(), tmp.end(), [&vec](size_t i1, size_t i2) { return vec[i1] < vec[i2]; });
    return tmp;
}

// Calculate median value of Marchenko-Pasture (MP) Distribution.
// beta: shape parameter of MPD.
// h: iteration limit.
double medianMPD(double beta, size_t h) {
    double bot = pow(1.0 - sqrt(beta), 2.0), top = pow(1.0 + sqrt(beta), 2.0);

    auto mp = [bot, top, beta](double x) {
        double tmp = (top - x) * (x - bot);
        return tmp > 0.0 ? sqrt(tmp) / (2.0 * boost::math::constants::pi<double>() * beta * x) : 0.0;
    };

    double lbnd = bot, ubnd = top;
    bool change = true;

    while (change & (ubnd - lbnd > 0.001)) {
        change = false;
        vector<double> x; // = getSeq(lbnd, ubnd, h);
        vector<double> lhalf, uhalf;

        for (size_t i = 0; i < h; ++i) {
            double y = 1.0 - boost::math::quadrature::trapezoidal(mp, x[i], top);

            if (y < 0.5) lhalf.push_back(x[i]);
            if (y > 0.5) uhalf.push_back(x[i]);
        }

        if (lhalf.size() > 0) {
            lbnd = *max_element(lhalf.begin(), lhalf.end());
            change = true;
        }

        if (uhalf.size() > 0) {
            ubnd = *min_element(uhalf.begin(), uhalf.end());
            change = true;
        }
    }

    return (ubnd + lbnd) / 2.0;
}

// Return prediction of mapping.
// C: coefficient matrix.
// X: data matrix.
// [[Rcpp::export]]
VectorXd getPred(const MatrixXd& C, const MatrixXd& X) {
    return C.cwiseProduct(X).rowwise().sum();
}

// Return RMSE between v1 and v2.
// [[Rcpp::export]]
double getRMSE(const VectorXd& v1, const VectorXd& v2) {
    return sqrt(pow((v1 - v2).array(), 2.0).mean());
}

// Return MAE between v1 and v2.
// [[Rcpp::export]]
double getMAE(const VectorXd& v1, const VectorXd& v2) {
    return abs((v1 - v2).array()).mean();
}

// Return correlation coefficient between v1 and v2.
// [[Rcpp::export]]
double getRho(const VectorXd& v1, const VectorXd& v2) {
    ArrayXd a1 = v1.array() - v1.mean(), a2 = v2.array() - v2.mean();
    return (a1 * a2).sum() / sqrt(pow(a1, 2.0).sum() * pow(a2, 2.0).sum());
}

// Calculate Normal Euclidean Distance (NE-Dist) matrix.
// [[Rcpp::export]]
MatrixXd getNEDist(const MatrixXd& X) {
    size_t m = X.rows();
    MatrixXd dist = MatrixXd::Zero(m, m);

    for (size_t i = 0; i < m - 1; ++i) {
        for (size_t j = i + 1; j < m; ++j) {
            auto d = sqrt((X.row(i) - X.row(j)).squaredNorm());
            dist(i, j) = dist(j, i) = d;
        }
    }

    return dist;
}

// Calculate Geodesic Distance matrix.
// [[Rcpp::export]]
Eigen::MatrixXd getGeodesDist(const Eigen::MatrixXd& dmat, size_t k) {
    size_t m = dmat.rows(), cnt;
    Eigen::MatrixXd kmat = Eigen::MatrixXd::Zero(m, m);

    // calculate KNN by data points
    for(size_t i = 0; i < m; ++i) {
        // sort by distance
        std::vector<std::pair<size_t, double>> knn(m);

        for(size_t j = 0; j < m; ++j) {
            knn[j] = std::pair<size_t, double>(j, dmat(i, j));
        }

        std::sort(knn.begin(), knn.end(), PrInc());

        // crate a matrix
        cnt = 0;

        for(size_t j = 0; j < m && cnt < k; ++j) {
            if(i != knn[j].first) {
                kmat(i, knn[j].first) = knn[j].second;
                cnt++;
            }
        }
    }

    // symmetrizatino of matrix
    for(size_t i = 0; i < m - 1; ++i) {
        for(size_t j = i + 1; j < m; ++j) {
            if(kmat(i, j) == 0 && kmat(j, i) != 0) {
                kmat(i, j) = kmat(j, i);
            } else if(kmat(i, j) != 0 && kmat(j, i) == 0) {
                kmat(j, i) = kmat(i, j);
            }
        }
    }


    // calculate minimum distance between each point (Warshall-Floyd algorithm)
    for(size_t i = 0; i < m; ++i) {
        for(size_t j = 0; j < m; ++j) {
            for(size_t k = 0; k < m; ++k) {
                if(kmat(j, i) != 0 && kmat(i, k) != 0) {
                    if(j != k && ((kmat(j, k) == 0) || (kmat(j, i) + kmat(i, k) < kmat(j, k)))) {
                        kmat(j, k) = kmat(j, i) + kmat(i, k);
                    }
                }
            }
        }
    }

    return kmat;
}

// Calculate Local Manifold Distance (LM-Dist) matrix.
// Locality parameters are controlled by theta<MatrixXd>
// [[Rcpp::export]]
MatrixXd getLMDist(const MatrixXd& X, const MatrixXd& theta) {
    size_t m = X.rows();
    MatrixXd X_ = X * theta, dist = MatrixXd::Zero(m, m);

    for (size_t i = 0; i < m - 1; ++i) {
        for (size_t j = i + 1; j < m; ++j) {
            auto d = sqrt((X_.row(i) - X_.row(j)).squaredNorm());
            dist(i, j) = dist(j, i) = d;
        }
    }

    return dist;
}

// Implementation of RngUnif class.
RngUnif::RngUnif() : rng(boost::mt19937(time(NULL)), Dist(0.0, 1.0)) {}
RngUnif::RngUnif(int seed) : rng(boost::mt19937(seed), Dist(0.0, 1.0)) {}
RngUnif::RngUnif(double min, double max) : rng(boost::mt19937(time(NULL)), Dist(min, max)) {}
RngUnif::RngUnif(int seed, double min, double max) : rng(boost::mt19937(seed), Dist(min, max)) {}

double RngUnif::operator()() { return rng(); }

// Implementation of RngNorm class.
RngNorm::RngNorm() : rng(boost::mt19937(time(NULL)), Dist(0.0, 1.0)) {}
RngNorm::RngNorm(int seed) : rng(boost::mt19937(seed), Dist(0.0, 1.0)) {}
RngNorm::RngNorm(double mean, double sigma) : rng(boost::mt19937(time(NULL)), Dist(mean, sigma)) {}
RngNorm::RngNorm(int seed, double mean, double sigma) : rng(boost::mt19937(seed), Dist(mean, sigma)) {}

double RngNorm::operator()() { return rng(); }

