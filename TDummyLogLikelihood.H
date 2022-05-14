#ifndef TDummyLogLikelihood_H_seen
#define TDummyLogLikelihood_H_seen

#include <TMatrixD.h>
#include <TVectorD.h>
#include <cmath>

// A dummy log likelihood for testing.  This is an example, but don't
// slavishly copy it (or you will be sorry!).  This is a quadratic.  The
// correlations between the different dimensions are set during Init().  See
// the "#define" values there for more information.
class TDummyLogLikelihood {
public:
    // Determine the number of dimensions.  This is where the dimensions are
    // defined, and everything else uses it.
    std::size_t GetDim() const {return 100;}

    // Calculate the log(likelihood).  The dummy likelihood is a Gaussian
    // (with covariance) centered at zero.  The covariance is set in Init()
    // (below).
    double operator()(const Vector& point)  const {
        double logLikelihood = 0.0;

        for (std::size_t i = 0; i<GetDim(); ++i) {
            for (std::size_t j = 0; j<GetDim(); ++j) {
                logLikelihood -= 0.5*point[i]*Error(j,i)*point[j];
            }
        }

        return logLikelihood;
    }

    // Note that this needs to be the grad(log(Likelihood)).
    bool operator() (Vector& g, const Vector& p) {
        for (std::size_t i=0; i<p.size(); ++i) {
            g[i] = 0.0;
            for (std::size_t j=0; j<p.size(); ++j) {
                g[i] -= Error(i,j)*p[j];
            }
        }
        return true;
    }

    void Init() {
        Covariance.ResizeTo(GetDim(),GetDim());
        Error.ResizeTo(GetDim(),GetDim());

        // Set the sigma for each variable.
        for (std::size_t i = 0; i<GetDim(); ++i) {
            double sigma = 1.0;
            // double sigma = 1.0*i + 1.0;
            Covariance(i,i) = sigma*sigma;
        }

        for (std::size_t i = 0; i<GetDim(); ++i) {
            for (std::size_t j = i+1; j<GetDim(); ++j) {
                double sig1 = std::sqrt(Covariance(i,i));
                double sig2 = std::sqrt(Covariance(j,j));
                // Now give some correlations to the likelihood.  (Uncomment
                // the one you want to try).

                // Choose a random correlation
// #define RANDOM_CORRELATION 0.9
#ifdef RANDOM_CORRELATION
#define SET_CORRELATION
#warning Using random correlations with RANDOM_CORRELATION
                {
                    double c = 0.3*RANDOM_CORRELATION*gRandom->Gaus();
                    if (c > RANDOM_CORRELATION) c = RANDOM_CORRELATION;
                    if (c < -RANDOM_CORRELATION) c = -RANDOM_CORRELATION;
                    Covariance(i,j) = c*sig1*sig2;
                }
#endif

                // Choose a correlation based on the variables.  The matrix
                // forms and "X" with large correlation along the crossing
                // diagonal.
#define VERY_CORRELATED 0.99
#ifdef VERY_CORRELATED
#define SET_CORRELATION
#warning Using very correlated correlations with VERY_CORRELATED
                if ((i+j) == GetDim()-1) {
                    int d = i-j;
                    if (std::abs(d) < GetDim()-1) continue;
                    Covariance(i,j)
                        = VERY_CORRELATED*sig1*sig2*(j - i)/(GetDim()-1.0);
                }
#endif

                // Make a tri-diagonal correlation.
// #define TRIDIAGONAL_CORRELATED 0.5
#ifdef TRIDIAGONAL_CORRELATED
#define SET_CORRELATION
#warning Using tridiagonal correlations with TRIDIAGONAL_CORRELATED
                if (j == i+4) {
                    Covariance(i,j) = TRIDIAGONAL_CORRELATED*sig1*sig2;
                }
#endif

                // Choose no correlation
#ifndef SET_CORRELATION
#warning Using no correlations.
                Covariance(i,j) = 0.0;
#endif

                // Make sure it's symmetric
                Covariance(j,i) = Covariance(i,j);
            }
        }

        if (GetDim() < 6) Covariance.Print();

        Error = Covariance;

        // Make sure the covariance is positive definite.
        double reduction = 1.0;
        do {
            TVectorD eigenValues(GetDim());
            Error.EigenVectors(eigenValues);
            bool positiveDefinite = true;
            for (std::size_t i = 0; i<GetDim(); ++i) {
                if (eigenValues(i)<0.0) {
                    positiveDefinite = false;
                }
            }
            if (positiveDefinite) break;
            std::cout << "Not positive definite" << std::endl;
            reduction *= 0.9;
            for (std::size_t i = 0; i<GetDim(); ++i) {
                for (std::size_t j = i+1; j<GetDim(); ++j) {
                    Error(i,j) = reduction*Covariance(i,j);
                    Error(j,i) = reduction*Covariance(i,j);
                }
            }
        } while (true);
        std::cout << "Likelihood correlation reduction of " << reduction
                  << std::endl;

        Error.Invert();
    }

    static TMatrixD Covariance;
    static TMatrixD Error;
};
TMatrixD TDummyLogLikelihood::Covariance;
TMatrixD TDummyLogLikelihood::Error;
#endif
