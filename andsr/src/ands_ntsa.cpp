#include "ands_ntsa.h"

using namespace std;
using namespace boost;
using namespace Eigen;
using namespace Rcpp;
using namespace RcppEigen;

// Calculate optimal hard threashold for singular values.
// Code is adopted from Gavish & Danoho (2014), IEEE Trans Inf Theor.
// [[Rcpp::export]]
size_t getOptHT(const VectorXd& sv, size_t m, size_t n) {
    double beta = m < n ? (double) m / n : (double) n / m;
    double coef = sqrt(2 * (beta + 1) + 8 * beta / (beta + 1 + sqrt(pow(beta, 2) + 14 * beta + 1)));
    double cut = coef / sqrt(medianMPD(beta)) * (n % 2 == 1 ? sv((n - 1) / 2) : (sv(n / 2 - 1) + sv(n / 2)) / 2);
    return (sv.array() > cut).count();
}

// Generate local seedds for multi-thread computing.
vector<int> genSeeds(size_t threadNo, int gSeed) {
    boost::mt19937 mt = IntegerVector::is_na(gSeed) ? boost::mt19937(time(NULL)) : boost::mt19937(gSeed);

    vector<int> seeds(threadNo);
    for (int& s : seeds) s = mt();
    return seeds;
}

// Calculate distance average of NED.
VectorXd getDBar(const MatrixXd& X, const vector<size_t>& idx) {
    MatrixXd dist = getNEDist(X);
    dist = matSubset(dist, idx, false);
    return dist.rowwise().mean();
}

// Calculate weight matrix.
MatrixXd getWMat(const MatrixXd& X, const VectorXd& dbar, const MatrixXd& theta) {
    MatrixXd wmat = getLMDist(X.leftCols(theta.cols()), theta);
    wmat = dbar.cwiseInverse().asDiagonal() * wmat;
    wmat = exp(-wmat.array());
    wmat.diagonal().setZero();

    return wmat;
}

// Perform S-map with SVD.
// Note: Add intersept term to the last column of 'X' before run.
// [[Rcpp::export]]
MatrixXd smapSVD(const MatrixXd& X, const RowVectorXd& y, const MatrixXd& wmat) {
    size_t m = wmat.rows(), n = X.cols(), mx = X.rows(), r;
    MatrixXd X_, C(m, n), Sv(n, n);
    VectorXd y_, sv(n);

    for (size_t tar = 0; tar < m; ++tar) {
        X_ = wmat.row(tar).asDiagonal() * X;
        y_ = wmat.row(tar).cwiseProduct(y);

        // Perform truncated SVD
        JacobiSVD<MatrixXd> svd(X_, ComputeThinU | ComputeThinV);
        sv = svd.singularValues();
        r = (sv.array() > sv(0) * 1e-5).count();    // r = getOptHT(sv, m, n);
        Sv = (sv.array().inverse()).matrix().asDiagonal();
        C.row(tar) = svd.matrixV().block(0, 0, n, r) * Sv.block(0, 0, r, r) * svd.matrixU().block(0, 0, mx, r).transpose() * y_;
    }

    return C;
}

// Perform S-map with L2 regularization.
// Note:: Add intersept term to the last column of 'X' before run.
// [[Rcpp::export]]
MatrixXd smapRidge(const MatrixXd& X, const RowVectorXd& y, const MatrixXd& wmat, double lambda) {
    size_t m = wmat.rows(), n = X.cols();
    MatrixXd X_(n, n), C(m, n), Lambda = lambda * MatrixXd::Identity(n, n);
    VectorXd y_(X.rows());

    for (size_t tar = 0; tar < m; ++tar) {
        X_ = X.transpose() * wmat.row(tar).asDiagonal() * X + Lambda;
        y_ = wmat.row(tar).cwiseProduct(y);

        FullPivLU<MatrixXd> lu(X_);
        C.row(tar) = lu.solve(X.transpose()) * y_;
    }

    return C;
}

// LMD-based prediction with simplex projection
// [[Rcpp::export]]
VectorXd lmdSmplx(const MatrixXd& X, const RowVectorXd& y, const MatrixXd& dmat,
                  const vector<size_t>& idxLib, const vector<size_t>& idxPrd, size_t nns) {
    // Implementation of NNP (Nearest Neighbour Points) search
    auto find_nn_pts = [&idxLib, &nns](size_t tar, const vector<double>& dist) {
        vector<size_t> idx = sortID(dist, idxLib);
        vector<size_t> nn_pts;

        for (size_t i : idx) {
            if (i != tar) nn_pts.push_back(i);
            if (nn_pts.size() == nns) break;
        }

        double maxD = nn_pts.back();

        // Check whether same distance exists after nns'th state
        for (size_t i = nns; i < idx.size(); ++i) {
            if (dist[idx[i]] != maxD) break;
            nn_pts.push_back(idx[i]);
        }

        return nn_pts;
    };

    size_t xm = X.rows();
    vector<double> pred;

    for (size_t tar : idxPrd) {
        vector<double> dist(xm);
        Map<VectorXd>(&dist[0], xm) = dmat.row(tar);

        vector<size_t> nn_pts = find_nn_pts(tar, dist);
        double min_dist = dist[nn_pts[0]], tmp = 0.0;
        size_t len = nn_pts.size();
        VectorXd weights = VectorXd::Zero(len);

        if (min_dist == 0.0) {
            for (size_t i = 0; i < len; ++i) {
                if (dist[nn_pts[i]] == min_dist) {
                    weights(i) = 1.0;
                    tmp += y(nn_pts[i]);
                } else {
                    break;
                }
            }
        } else {
            for (size_t i = 0; i < len; ++i) {
                weights(i) = exp(-dist[nn_pts[i]] / min_dist);
                tmp += weights(i) * y(nn_pts[i]);
            }
        }

        pred.push_back(tmp / weights.sum());
    }

    return Map<VectorXd>(&pred[0], pred.size());
}

// LMD-based prediction with S-map.
// [[Rcpp::export]]
MatrixXd lmdSMap(const MatrixXd& X, const RowVectorXd& y,
                 const vector<size_t>& idxLib, const vector<size_t>& idxPrd,
                 int method, const MatrixXd& theta, double lambda) {
    size_t m = idxPrd.size(), n = X.cols();
    // MatrixXd Xl = X(idxLib, all), Xp = X(idxPrd, all), C(m, n);
    // RowVectorXd yl = y(idxLib);
    MatrixXd Xl = matSubset(X, idxLib, true), Xp = matSubset(X, idxPrd, true), C(m, n);
    RowVectorXd yl = vecSubset(y, idxLib);
    VectorXd dbar = getDBar(X.leftCols(n - 1), idxLib);

    // Make weight matrix
    MatrixXd wmat = getWMat(X.leftCols(n - 1), dbar, theta);
    wmat = matSubset(wmat, idxPrd, idxLib); // wmat = wmat(idxPrd, idxLib);

    if (method == 0) {
        C = smapSVD(Xl, yl, wmat);
    } else {
        C = smapRidge(Xl, yl, wmat, lambda);
    }

    return C;
}

// Imprementation of LMD parameter search with Temperature-Parallel Simulated Annealing (TPSA).
// [[Rcpp::export]]
List funcTPSA(int gSeed, size_t threadNo, size_t iterLim, size_t tsLength, int criterion,
                  int method, int diag, const vector<double>& sigmas, const vector<double>& temps,
                  const MatrixXd& X, const RowVectorXd& y, const vector<size_t>& idxLib, const vector<size_t>& idxPrd) {
    // Return probability of state exchange
    auto prob = [](double curT, double nexT, double curE, double nexE) {
        double dT = nexT - curT, dE = nexE - curE;
        return dT * dE < 0.0 ? 1.0 : exp(-dT * dE / (nexT * curT));
    };

    // Preparation for TPSA analysis
    RngUnif rUnif(gSeed, 0.0, 1.0);
    vector<int> seeds = genSeeds(threadNo, gSeed);
    vector<SolveSA> vecSA;
    vector<double> vecTemp = getSeq(temps[1], temps[0], threadNo);

    // Preparation of library and prediction dataset
    size_t n = X.cols();
    // MatrixXd Xl = X(idxLib, all), Xp = X(idxPrd, all);
    // RowVectorXd yl = y(idxLib), yp = y(idxPrd);
    MatrixXd Xl = matSubset(X, idxLib, true), Xp = matSubset(X, idxPrd, true);
    RowVectorXd yl = vecSubset(y, idxLib), yp = vecSubset(y, idxPrd);
    VectorXd dbar = getDBar(X.leftCols(n - 1), idxLib);

    // Cost evaluation of initial parameters
    double iniLambda = 0.0;
    MatrixXd iniTheta = MatrixXd::Identity(n - 1, n - 1);
    MatrixXd wmat = matSubset(getWMat(X.leftCols(n - 1), dbar, iniTheta), idxPrd, idxLib), C;

    VectorXd pred;

    if (method == 0) {
        C = smapSVD(Xl, yl, wmat);
    } else {
        C = smapRidge(Xl, yl, wmat, iniLambda);
    }

    double eval = SolveSA::evalCost(criterion, C, Xp, yp);

    // Make SA Solver
    for (size_t i = 0; i < threadNo; ++i) {
        SolveSA sa(seeds[i], sigmas, iniTheta, iniLambda, vecTemp[i], eval);
        vecSA.push_back(sa);
    }

    // Perform TPSA
    for (size_t iter = 0; iter < iterLim; ++iter) {
        // Run threads
        vector<thread> threads;

        for (size_t i = 0; i < threadNo; ++i) {
            threads.push_back(thread(ref(vecSA[i]), tsLength, criterion, method, diag, ref(X),
                                     ref(Xl), ref(Xp), ref(yl), ref(yp),
                                     ref(dbar), ref(idxLib), ref(idxPrd)));
        }

        for (thread& th : threads) th.join();

        // State exchange
        for (size_t i = 0; i < threadNo - 1; ++i) {
            double curT = vecSA[i].getTemp(), nexT = vecSA[i + 1].getTemp();
            double curE = vecSA[i].getEval(), nexE = vecSA[i + 1].getEval();

            if (rUnif() <= prob(curT, nexT, curE, nexE)) vecSA[i].paramSwap(vecSA[i + 1]);
        }
    }

    // Finalization of TPSA
    vector<double> threadEval;
    for (size_t i = 0; i < threadNo; ++i) threadEval.push_back(vecSA[i].getEval());

    vector<double>::iterator minIter = min_element(threadEval.begin(), threadEval.end());
    size_t minID = std::distance(threadEval.begin(), minIter);
    List L = List::create(Named("Eval") = vecSA[minID].getEval(), Named("Theta") = vecSA[minID].getTheta(), Named("Lambda") = vecSA[minID].getLambda());
    return L;
}

// Implementation of SolveSA class.
SolveSA::SolveSA(int seed, const vector<double>& sigmas, const MatrixXd& iniTheta, double iniLambda, double iniTemp, double iniEval) :
    curTheta(iniTheta), curLambda(iniLambda), curTemp(iniTemp), curEval(iniEval)
{
    // RNG initialization
    rUnif = RngUnif(seed);
    rNorm0 = RngNorm(seed, 0.0, sigmas[0]);
    rNorm1 = RngNorm(seed, 0.0, sigmas[1]);
}

// Mutation event of theta.
MatrixXd SolveSA::mutTheta(int diag) {
    MatrixXd theta = curTheta;
    size_t n = theta.cols();

    // diag == 0: mutation events occur only at diagonal elements
    if (diag == 0) {
        for (size_t i = 0; i < n; ++i) theta(i, i) += rNorm0();
    } else {
        for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) theta(i, j) += rNorm0();
    }

    return theta.array().abs();
}

// Mutation event of lambda.
double SolveSA::mutLambda() {
    double lambda = curLambda;
    lambda += rNorm1();
    return lambda < 0.0 ? 0.0 : lambda;
}

// Cost evaluation function.
double SolveSA::evalCost(int criterion, const MatrixXd& C, const MatrixXd& Xp, const RowVectorXd& yp) {
    VectorXd pred = getPred(C, Xp);
    double eval;

    switch(criterion) {
        case 0:
            eval = getRMSE(pred, yp);
            break;
        case 1:
            eval = getMAE(pred, yp);
            break;
        default:
            eval = -getRho(pred, yp);
            break;
    }

    return eval;
}

// Exchange of parameters with other SolveSA object.
void SolveSA::paramSwap(SolveSA& sa) {
    MatrixXd tmpTheta(curTheta);
    double tmpLambda(curLambda);
    double tmpEval(curEval);

    this->setTheta(sa.getTheta());
    this->setLambda(sa.getLambda());
    this->setEval(sa.getEval());

    sa.setTheta(tmpTheta);
    sa.setLambda(tmpLambda);
    sa.setEval(tmpEval);
}

// Implementation of SA process.
void SolveSA::operator()(size_t len, int criterion, int method, int diag,
                         const MatrixXd& X, const MatrixXd& Xl, const MatrixXd& Xp,
                         const RowVectorXd& yl, const RowVectorXd& yp, const VectorXd& dbar,
                         const vector<size_t>& idxLib, const vector<size_t>& idxPrd) {
    if (method == 0) {  // Perform smapSVD with LMD
        for (size_t i = 0; i < len; ++i) {
            MatrixXd nexTheta = mutTheta(diag);

            MatrixXd wmat = matSubset(getWMat(X, dbar, nexTheta), idxPrd, idxLib);
            MatrixXd C = smapSVD(Xl, yl, wmat);
            double nexEval = SolveSA::evalCost(criterion, C, Xp, yp);

            if (rUnif() <= probSwap(curEval, nexEval)) {
                curTheta = nexTheta;
                curEval = nexEval;
            }
        }
    } else {            // Perform smapRidge with LMD
        for (size_t i = 0; i < len; ++i) {
            MatrixXd nexTheta = mutTheta(diag);
            double nexLambda = mutLambda();

            MatrixXd wmat = matSubset(getWMat(X, dbar, nexTheta), idxPrd, idxLib);
            MatrixXd C = smapRidge(Xl, yl, wmat, nexLambda);
            double nexEval = SolveSA::evalCost(criterion, C, Xp, yp);

            if (rUnif() <= probSwap(curEval, nexEval)) {
                curTheta = nexTheta;
                curLambda = nexLambda;
                curEval = nexEval;
            }
        }
    }
}
