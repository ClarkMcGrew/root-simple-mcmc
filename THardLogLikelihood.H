#ifndef THardLogLikelihood_H_seen
#define THardLogLikelihood_H_seen

#include <TMatrixD.h>
#include <TVectorD.h>

// A hard log likelihood for testing.  This is an example, but don't slavishly
// copy it (or you will be sorry!).  This calculates a multi-dimensional
// variant of the Rosenbrock function with the caveat that it returns the
// opposite of the usual function.  The maximum value is zero, and the
// function is always zero or less than zero.
//
// The Rosenbrock function is designed to test optimization problemsn and is
// non-covariant with long curved valleys.  In two dimensions the Rosenbrock
// function is
//
// f = (R_a - x)**2 + R_b*(y - x**2)**2
//
// which has the global minimum at (R_a, R_a**2).  Both R_a and R_b must be
// positive.  The value of R_a is fixed to 1.0.  The value of R_b is set using
// ROSEN_B below.
//
// There are several forms for the multi-dimensional Rosenbrock function.
// This particular version is suggested by SciPy.  It has a minimum at x =
// [1,1,...,1], and, possibly, other local minima (e.g. near [-1,1,1,...,1]).
// The global minimum might be degenerate.  For D=2 (and D=3), it has one global
// minimum at [1,1] (and [1,1,1]).  For D = 4 to 7 there is a proof that it
// has one global at [1,1,..,1,1] , and one local minimum near [-1,1,...,1].
//
// COMMENT: The Rosenbrock function is a good example of one of the
// idiosyncracies of Bayesian analysis vs parameter optimization.  At high
// dimensionaities, there is a global best fit point at (1,1,...1,1,1) and the
// probability dropping rapidly toward (1,1,...,1,1,0).  Surprisingly, when
// you marginalize to find the most probable value of the last variable,
// because of the behavior away from best fit, the most probable marginalized
// value is very close to zero (NOT ONE, WHICH IS THE MOST PROBABLE POINT).
class THardLogLikelihood {
public:
    // Determine the number of dimensions.  This is where the dimensions are
    // defined, and everything else uses it.  The hard likelihood is only
    // defined for two or more dimensions.
    std::size_t GetDim() const {return 6;}

    // Set the next definition to a positive value to change the curvature of
    // the Rosenbrock function near the minimum.  The TRADITIONAL value is
    // ROSEN_B equal to 100, but the function is OK for any positive value.
    // Larger values make valleys in the function to be deeper and more
    // narrow.  Smaller values make the function easier to manage since it's
    // smaller second derivatives.
#ifdef ROSEN_B
#undef ROSEN_B
#endif
#define ROSEN_B 100.0

    // Calculate the log(likelihood).  The hard likelihood is the opposite of
    // a Rosenbrock function and has a maxima of 0 at X=(1,1,1,...,1).
    double operator()(const sMCMC::Vector& point)  const {
        double logLikelihood = 0.0;

        for (std::size_t i = 0; i<GetDim()-1; ++i) {
            double a = (1.0-point[i]);
            double b = point[i+1] - point[i]*point[i];
            logLikelihood -= a*a + ROSEN_B*b*b;
        }

        return logLikelihood;
    }

    // Note that this needs to be the grad(log(Likelihood)).
    bool operator() (sMCMC::Vector& g, const sMCMC::Vector& p) {

        // Find gradient for first element.
        g[0] = - 2.0*(1.0-p[0]) - 4.0*ROSEN_B*p[0]*(p[1]-p[0]*p[0]);

        // Find gradients for middle elements.
        for (std::size_t i = 1; i<GetDim()-1; ++i) {
            g[i] = 2.0*ROSEN_B*(p[i] - p[i-1]*p[i-1]);
            g[i] += - 2.0*(1.0-p[i]);
            g[i] += - 4.0*ROSEN_B*p[i]*(p[i+1] - p[i]*p[i]);
         }

        // Find gradient for last element.
        std::size_t i = GetDim()-1;
        g[i] = + 2.0*ROSEN_B*(p[i]-p[i-1]*p[i-1]);

        // Take the opposite since the returned value above is the opposite of
        // the Rosenbrock function.
        for (std::size_t i = 0; i<GetDim(); ++i) g[i] = -g[i];

        return true;
    }

    // Add a mentod to mimic the TDummyLikelyhood.H behavior.  This isn't
    // doing anything!
    void Init() {}

    // Here to match TDummyLikelyhood.H.  They aren't used.
    static TMatrixD Covariance;
    static TMatrixD Error;
};
TMatrixD THardLogLikelihood::Covariance;
TMatrixD THardLogLikelihood::Error;
#endif
