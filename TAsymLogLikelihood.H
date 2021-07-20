#ifndef TASymLogLikelihood_H_seen
#define TASymLogLikelihood_H_seen

#include <TMatrixD.h>
#include <TVectorD.h>

// A very asymmetric log likelihood for testing.  This is an example, but
// don't slavishly copy it (or you will be sorry!).
class TASymLogLikelihood {
public:
    // Determine the number of dimensions.  This is where the dimensions are
    // defined, and everything else uses it.  The hard likelihood is only
    // defined for two or more dimensions.
    std::size_t GetDim() const {return 100;}

    const double positiveSlope = -1.0;
    const double negativeSlope = 100.0;

    // Calculate the log(likelihood).
    double operator()(const Vector& point)  const {
        double logLikelihood = 0.0;

        for (std::size_t i = 0; i<GetDim(); ++i) {
            double a = point[i];
            if (a<0.0) a *= negativeSlope;
            else a *= positiveSlope;
            logLikelihood += a;
        }

        return logLikelihood;
    }

    // Note that this needs to be the grad(log(Likelihood)).
    bool operator() (Vector& g, const Vector& p) {
        throw;
    }

    // Add a mentod to mimic the TDummyLikelyhood.H behavior.  This isn't
    // doing anything!
    void Init() {}

    // Here to match TDummyLikelyhood.H.  They aren't used.
    static TMatrixD Covariance;
    static TMatrixD Error;
};
TMatrixD TASymLogLikelihood::Covariance;
TMatrixD TASymLogLikelihood::Error;
#endif
