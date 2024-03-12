#ifndef ANDS_NTSA_H
#define ANDS_NTSA_H

#include <thread>
#include "ands_common.h"

// Forward declarations.
size_t getOptHT(const Eigen::VectorXd& sv, size_t m, size_t n);

std::vector<int> genSeeds(size_t threadNo, int gSeed);

Eigen::VectorXd getDBar(const Eigen::MatrixXd& X, const std::vector<size_t>& idx);

Eigen::MatrixXd getWMat(const Eigen::MatrixXd& X, const Eigen::VectorXd& dbar, const Eigen::MatrixXd& theta);

Eigen::MatrixXd smapSVD(const Eigen::MatrixXd& X, const Eigen::RowVectorXd& y, const Eigen::MatrixXd& wmat);

Eigen::MatrixXd smapRidge(const Eigen::MatrixXd& X, const Eigen::RowVectorXd& y, const Eigen::MatrixXd& wmat, double lambda);

Eigen::VectorXd lmdSmplx(const Eigen::MatrixXd& X, const Eigen::RowVectorXd& y, const Eigen::MatrixXd& dmat, const std::vector<size_t>& idxLib, const std::vector<size_t>& idxPrd, size_t nns);

Eigen::MatrixXd lmdSMap(const Eigen::MatrixXd& X, const Eigen::RowVectorXd& y, const std::vector<size_t>& idxLib, const std::vector<size_t>& idxPrd, int method, const Eigen::MatrixXd& theta, double lambda);

Rcpp::List funcTPSA(int gSeed, size_t threadNo, size_t iterLim, size_t tsLength, int criterion, int method, int diag,
                    const std::vector<double>& sigmas, const std::vector<double>& temps, const Eigen::MatrixXd& X,
                    const Eigen::RowVectorXd& y, const std::vector<size_t>& idxLib, const std::vector<size_t>& idxPrd);

class SolveSA
{
public:
    // Constructor
    SolveSA(int seed, const std::vector<double>& sigmas, const Eigen::MatrixXd& iniTheta, double iniLambda, double iniTemp, double iniEval);

    // Accessor
    void setTheta(const Eigen::MatrixXd& theta) { curTheta = theta; }
    void setLambda(double lambda) { curLambda = lambda; }
    void setEval(double eval) { curEval = eval; }

    Eigen::MatrixXd getTheta() { return curTheta; }
    double getLambda() { return curLambda; }
    double getTemp() { return curTemp; }
    double getEval() { return curEval; }

    // Calculate exchange probability b/w current & mutate parameters at temperature 'curTemp'.
    double probSwap(double curE, double nexE) {
        double deltaE = nexE - curE;
        double prob = deltaE < 0.0 ? 1.0 : exp(-(deltaE) / curTemp);
        return prob;
    }

    // Mutation function
    Eigen::MatrixXd mutTheta(int diag);
    double mutLambda();

    // Cost evaluation
    static double evalCost(int criterion, const Eigen::MatrixXd& C, const Eigen::MatrixXd& Xp, const Eigen::RowVectorXd& yp);

    // Exchange of parameters with other SolveSA object
    void paramSwap(SolveSA& sa);

    // Implementation of SA process
    void operator()(size_t len, int criterion, int method, int diag,
                    const Eigen::MatrixXd& X, const Eigen::MatrixXd& Xl, const Eigen::MatrixXd& Xp,
                    const Eigen::RowVectorXd& yl, const Eigen::RowVectorXd& yp, const Eigen::VectorXd& dbar,
                    const std::vector<size_t>& idxLib, const std::vector<size_t>& idxPrd);

private:
    // Arguments
    Eigen::MatrixXd curTheta;
    double curLambda;
    double curTemp;
    double curEval;

    // RNG objects
    RngUnif rUnif;
    RngNorm rNorm0;
    RngNorm rNorm1;
};

#endif
