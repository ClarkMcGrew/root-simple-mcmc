#ifndef THardLogLikelihood_H_seen
#define THardLogLikelihood_H_seen

#include <TMatrixD.h>
#include <TVectorD.h>

// A ridiculously "wrong" log likelihood for testing.  This is an example, but
// don't slavishly copy it (or you will be sorry!).  This calculates a
// multi-dimensional likelihood where the variables are completely correlated
// and non-Gaussian.
//
// The unnormalized likelihood calculated is:
//
// L = exp(- 0.5*(x+y+x+...)**2 / sigma**2)
//
// which dones not have a global maximum (and is degenerate for "x+y+z+..."
// equal to zero).
class THorrificLogLikelihood {
public:
    // Determine the number of dimensions.  This is where the dimensions are
    // defined, and everything else uses it.  This gets very slow with more
    // than 100 dimensions.
    std::size_t GetDim() const {return 75;}

    // Calculate the log(likelihood)
    double operator()(const sMCMC::Vector& point)  const {
        const double sigma = 0.01;
        double logLikelihood = 0.0;

        for (std::size_t i = 0; i<GetDim(); ++i) {
            if (std::abs(point[i]) > 1.0) return -1E+30;
            logLikelihood += point[i];
        }
        double naturalSigma = std::sqrt(GetDim()*4.0/12.0);
        logLikelihood /= naturalSigma;
        logLikelihood = - 0.5*logLikelihood*logLikelihood/sigma/sigma;
        return logLikelihood;
    }

    // Note that this should be the grad(log(Likelihood)).
    bool operator() (sMCMC::Vector& g, const sMCMC::Vector& p) {
        return false;
    }

    // Add a method to mimic the TDummyLikelyhood.H behavior.  This isn't
    // doing anything!
    void Init() {}

    // Here to match TDummyLikelyhood.H.  They aren't used.
    static TMatrixD Covariance;
    static TMatrixD Error;
};
TMatrixD THorrificLogLikelihood::Covariance;
TMatrixD THorrificLogLikelihood::Error;
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// End:
