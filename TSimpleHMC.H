#ifndef TSimpleHMC_H_SEEN
#define TSimpleHMC_H_SEEN

#include <iostream>
#include <vector>
#include <cmath>

#include <TRandom.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TVectorD.h>

// Define the amount of debugging when running the chain.
#ifndef HMC_DEBUG_LEVEL
#define HMC_DEBUG_LEVEL 2
#endif

// A macro to make outputting debug messages easy.  This can be used as
//
// HMC_DEBUG(1) << "This is debug output" << std::endl;
//
// The message will be printed if the value is less that HMC_DEBUG_LEVEL.
// Generically, the recommendad levels are "0" for once per run output, "1"
// for intermittent updates, "2" for more verbose intermittent updates, "3"
// for once per step, "4" for extra once per step verbosity, "5" for once per
// gradient, and "6" or more for super verbosity.
#ifndef HMC_DEBUG
#define HMC_DEBUG(level) if (level <= (HMC_DEBUG_LEVEL)) std::cout
#endif

// A macro to make outputting error messages easy.  It can be used like
// "std::cout".
#ifndef HMC_ERROR
#define HMC_ERROR (std::cout <<__FILE__<<":: " << __LINE__ << ": " )
#endif

namespace sMCMC {
    // Make a typedef for the type used for the function parameter.  The
    // parameter type will usually be a double, but for some problems that
    // might not be a good idea (e.g. uses too much memory, or it's not
    // needed).  This makes it easy to change between the two).
    typedef double Parameter;

    // Make a typedef for a vector of the parameters.  This is for the same
    // reason as "Parameter".
    typedef std::vector<Parameter> Vector;

    template <typename U, typename G> class TSimpleHMC;

};

namespace {
    // A place holder to flag that the user didn't provide a gradient.  It
    // will return "false" to indicate that the gradient was not calculated.
    struct SimpleHMCInvalidGradient {
        bool operator() (sMCMC::Vector&, const sMCMC::Vector&) {return false;}
    };
};

/// A class to run a Hamiltonian (or Hybrid) Monte Carlo.  This is implemented
/// based on Chapter 5 of the Handbook of Markov Chain Monte Carlo -- MCMC
/// Using Hamiltonian Dynamics by Radford M. Neal, arxiv:1206.1901.  The
/// UserLikelihood template argument must be a class (or struct) which
/// provides a method declared as:
///
/// The UserLikelihood template argument must be a class (or struct) which
/// provides a method declared as:
///
///\code
/// struct ExampleLogLikelihood {
///    double operator() (const Vector& point);
/// }
///\endcode
///
/// Note: the log likelihood should be returned, which means the value will be
/// negative.  Do not return the chi^2 (i.e. -2*log(likelihood).
///
/// The optional UserGradient will should calculate the gradient of the user
/// log(likelihood), and return true if the gradient is calculated.  if this
/// argument is available, they you should be using a pure HMC technique (and
/// not HMC).  If provided, the template argument must be a class (or struct)
/// which provides a method declared as:
///
///\code
/// struct ExampleGradient {
///    bool operator() (Vector& grad, const Vector& point);
/// }
///\endcode
///
/// If the gradient is not available, then you can use an approximate HMC
/// method.
//
/// BACKGROUND: The HMC technique is can very efficiently approximate the
/// posterior probability using a minimum number of "leapfrog" steps (See
/// "Chapter 5 of the Handbook of Markov Chain Monte Carlo -- MCMC Using
/// Hamiltonian Dynamics" by Radford M. Neal, arxiv:1206.1901).
/// Unfortunately, it's efficiency depends on the ability to calculate the
/// gradient at multiple points along the leapfrog step, and calculating the
/// gradient can be much more expensive than calculating the likelihood.
/// Fortunately, the correctness of the HMC technique depends on the
/// reversibility of the leapfrog step, and not the correctness of the
/// gradient calculation.  That means that an approximate gradient can be used
/// in place of the true gradient.  The cost of using an approximate gradient
/// is that the results of the HMC will not converge to the correct posterior
/// quite as efficiently, but the decrease in efficiency will often be
/// compensated for by the faster calculation of the approximate gradient
/// (relative to the true gradient).
///
/// The AHMC technique is based on the realization that in many cases, the
/// prior is approximately a multidimensional Gaussian, and the resulting
/// posterior is also approximately Gaussian.  Because of that, the gradient
/// of the approximate covariance of the posterior is often a pretty good
/// representation of the true posterior gradient.  The estimated posterior
/// will require more steps to approximate the true posterior, but the
/// increased number of steps is may be compenstated by the decreased time
/// required to calculate the gradient, and a more efficient coverage of the
/// posterior.
template <typename UserParameter,
          typename OptionalGradient = SimpleHMCInvalidGradient>
class sMCMC::TSimpleHMC {
public:
    /// Make the user likelihood class available as TSimpleHMC::LogLikelihood.
    typedef UserParameter LogLikelihood;

    /// Make optional gradient class available as TSimpleHMC::UserGradient.
    /// This is mostly here for debugging purposes.
    typedef OptionalGradient UserGradient;

    TSimpleHMC(TTree* tree = NULL, bool saveStep = false)
        : fTree(tree), fStepCount(0),
          fPotentialCount(0), fPotentialGradientCount(0),
          fLeapFrogSteps(10), fAlpha(0.0),
          fCovarianceWindow(1000000) {
        if (fTree) {
            HMC_DEBUG(0) << "TSimpleHMC: Adding branches to "
                         << fTree->GetName()
                         << std::endl;
            fTree->Branch("LogLikelihood",&fAcceptedPotential);
            fTree->Branch("Accepted",&fAccepted);
            fTree->Branch("Trace",&fCurrentCovarianceTrace);
            fTree->Branch("LikelihoodCalls", &fPotentialCount);
            fTree->Branch("Steps", &fStepCount);
            fTree->Branch("Acceptance", &fCurrentAcceptance);
            fTree->Branch("MeanEpsilon", &fMeanEpsilon);
            fTree->Branch("Orbit", &fEstimatedOrbitLength);
            fTree->Branch("Leapfrog", &fLeapFrogSteps);
        }
    }

    /// Get a reference to the likelihood calculation object.  The
    /// LogLikelihood object is constructed by the TSimpleHMC template and
    /// is the class that you handed to the template.  You have write access
    /// to that object with this method, and can access any of your declared
    /// methods.  Because this is a template, TSimpleHMC::LogLikelihood is
    /// actually just a typedef for your class.
    LogLikelihood& GetLogLikelihood() {return fLogLikelihood;}

    /// Get the acceptance rate
    double GetAcceptanceRate() { return fCurrentAcceptance; }

    /// Get a count of the total number of calls to the Potential method.
    int GetPotentialCount() const {
        return fPotentialCount;
    }

    /// Get a count of the total number of calls to the Potential method.
    int GetGradientCount() const {
        return fPotentialGradientCount;
    }

    /// Set the correlation between the last accepted momentum and the new
    /// proposed momentum.  See the ProposeMomentum() method for a description
    /// of Alpha.
    void SetAlpha(double a) {fAlpha = a;}

    /// Get the current value of alpha.
    double GetAlpha() const {return fAlpha;}

    /// Set the current value of epsilon
    void SetMeanEpsilon(double e) {fMeanEpsilon = e;}

    /// Get the current value of epsilon
    double GetMeanEpsilon() const {return fMeanEpsilon;}

    /// Override the automatically calculated number of leapfrog steps.  If
    /// this is zero, then a single approximate step is taken.  If this is
    /// positive, it overrides the automatically calculated number of steps.
    /// Any negative value will automatically calculate the number of steps.
    void SetLeapFrog(int i) {fLeapFrogSteps = -i;}

    /// Get the most recent central point.
    const Vector& GetCentralPoint() const {return fCentralPoint;}
    double GetCentralPotential() const {return fCentralPotential;}

    /// Get the most recently estimated covariance for the likelihood.
    const TMatrixD& GetEstimatedCovariance() const {
        return fEstimatedCovariance;}

    /// Set the last accepted value, and calculate the accepted potential.
    /// This is to allow the user to force a position between steps.
    void SetPosition(const Vector& start) {
        std::copy(start.begin(), start.end(), fAccepted.begin());
        fAcceptedPotential = Potential(start);
    }

    /// Set the starting point for the HMC chain.  If the optional argument is
    /// true, then the point will be saved to the output.  This must be called
    /// in the user code to initialize the HMC chain.
    void Start(const Vector& start, bool save=true) {
        fStepCount = 0;

        // Set the vector size.
        fProposed.resize(start.size());
        fProposedMomentum.resize(start.size());
        fAccepted.resize(start.size());
        fAcceptedMomentum.resize(start.size());
        fCentralPoint.resize(start.size());
        fAveragePoint.resize(start.size());

        SetPosition(start);

        // Use the starting point as the first proposal.
        std::copy(fAccepted.begin(), fAccepted.end(), fProposed.begin());
        fProposedPotential = fAcceptedPotential;

        // Set the value of the step size.  The step size will evolve based on
        // the walk through the log(likeihood).
        fMeanEpsilon = 0.05;
        fReversalLen = 0.0;

        // Set the target value for the step acceptance probability.  The
        // value of 65% comes from the HMC literature.
        fTargetAcceptance = 0.65;
        fCurrentAcceptance = fTargetAcceptance;

        // Set up the running estimation of the position of the minimum
        // potential.
        std::copy(fAccepted.begin(), fAccepted.end(), fCentralPoint.begin());
        fCentralPotential = fAcceptedPotential;

        // Set up the running average of the likelihood.
        std::copy(fAccepted.begin(), fAccepted.end(), fAveragePoint.begin());
        fAveragePointTrials = 0.0;

        // Set up a running calculation of the covariance.  This starts the
        // calculation with all the variances equal to 1, and no correlations
        // between the dimensions.  The initial guess is given a weight of 10
        // (each added point has a weight of 1).
        fEstimatedCovariance.ResizeTo(start.size(), start.size());
        fEXXT.ResizeTo(start.size(), start.size());
        fEstimatedError.ResizeTo(start.size(), start.size());
        fCovarianceTrials = 0;
        for (int i=0; i<start.size(); ++i) {
            for (int j=0; j<start.size(); ++j) {
                if (i == j) fEstimatedCovariance(i,j) = 1.0;
                else fEstimatedCovariance(i,j) = 0.0;
                fEXXT(i,j) = 0.0;
            }
        }
        fEstimatedError = fEstimatedCovariance;
        fEstimatedError.Invert();
        fEstimatedCovarianceTrace = start.size();

        fStepsRemaining = 0;
        fStepsSinceUpdate = 0;

        if (save) SaveStep();
    }

    /// Take a step.  This returns true if a new point has been accepted, and
    /// false if we stay at the old point.  If save is true, then the points
    /// are saved to the output.  This should be called repeatedly by the user
    /// code.  The gradient types are defined in PotentialGradient(), but the
    /// potentially useful values are 0) use the user defined or finite
    /// differences gradient, 2) use the "covariant" approximation for the
    /// gradient, or 5) for the gradient to zero.  The values of 2 and 5 are
    /// used to implement an AHMC (Approximate HMC).
    bool Step(bool save=true, int gradientType = 0) {
        if (fProposed.empty() || fProposedMomentum.empty()
            || fAccepted.empty() || fAcceptedMomentum.empty()) {
            HMC_ERROR << "Must initialize starting point" << std::endl;
            throw;
        }

        ++fStepCount;

        ProposeMomentum(fProposedMomentum,fAcceptedMomentum);

        // Save info needed about the starting point.  The starting position
        // is always fAccepted.
        double initialKinetic = KineticEnergy(fProposedMomentum);

        // Evolve the proposed position and momentum according to the
        // hamiltonial equations.  The step size is epsilon and the number of
        // steps will be fLeapFrogSteps.
        double epsilon = gRandom->Uniform(0.9*std::abs(fMeanEpsilon),
                                          1.1*std::abs(fMeanEpsilon));
        int okLeap = LeapFrog(fProposed,fProposedMomentum,fAccepted,
                               epsilon,std::abs(fLeapFrogSteps),gradientType);

        if (fLeapFrogSteps > 0) {
            if (okLeap != 2) {
                if (fMeanEpsilon > 0 && fReversalLen > fMeanEpsilon) {
                    double targetEpsilon = fReversalLen/8.0;
                    double deltaEpsilon = targetEpsilon - fMeanEpsilon;
                    if (deltaEpsilon > 0.0) fMeanEpsilon += 0.1*deltaEpsilon;
                }
                if (fLeapFrogSteps < 50) fLeapFrogSteps += 1;
            }
            else if (okLeap == 2) {
                if (fReversalLen < fMeanEpsilon) {
                    fReversalLen = std::abs(fLeapFrogSteps*epsilon);
                }
                else {
                    fReversalLen = 0.95*fReversalLen;
                    fReversalLen += 0.05*std::abs(fLeapFrogSteps*epsilon);
                }
                if (fLeapFrogSteps > 3) fLeapFrogSteps -= 1;
                if (fMeanEpsilon > 0) fMeanEpsilon *= 0.99;

            }
        }

        // Find the proposed kinetic and potential energy
        double proposedKinetic = KineticEnergy(fProposedMomentum);
        fProposedPotential = Potential(fProposed);

        // Find the change in energy between the proposed and accepted states.
        // If the leapfrog step was done with perfect accuracy, we would have
        // energy conservation, and the hamiltonian would remain constant
        // (i.e. delta = 0.0).
        double proposedHamiltonian = fProposedPotential + proposedKinetic;
        double acceptedHamiltonian = fAcceptedPotential + initialKinetic;

        if (okLeap && std::isfinite(fProposedPotential)) {
            // Update the running estimate of the covariance.
            UpdateCovariance(fAccepted, fAcceptedPotential,
                             fProposed, fProposedPotential);
            UpdateErrorMatrix();
        }
        else {
            if (fMeanEpsilon > 0) fMeanEpsilon = 0.3*fMeanEpsilon;
        }

        double delta = proposedHamiltonian - acceptedHamiltonian;
        double trial = - std::log(gRandom->Uniform());
        if (delta > trial || !std::isfinite(delta)) {
            HMC_DEBUG(3) << "Reject " << delta
                         << " " << proposedHamiltonian
                         << "(" << fProposedPotential
                         << ":" << proposedKinetic << ")"
                         << " <-- "
                         << acceptedHamiltonian
                         << "(" << fAcceptedPotential
                         << ":" << initialKinetic << ")"
                         << std::endl;
            // The proposed hamiltonian is more than the previously accepted
            // hamiltonian.  Reject it if it's not finite, or it isn't
            // accepted as a random "bad" step.  This depends on IEEE error
            // handling so that - std::log(0.0) is inf which is always more
            // than delta.  The Step failed so reverse direction.  This is
            // often overwritten by the next proposed step.
            for (std::size_t i=0; i<fAcceptedMomentum.size(); ++i) {
                fAcceptedMomentum[i] = -fAcceptedMomentum[i];
            }
            fCurrentAcceptance = (fCurrentAcceptance*4999.0)/5000.0;
        }
        else {
            HMC_DEBUG(3) << "Accept " << delta
                         << " " << proposedHamiltonian
                         << "(" << fProposedPotential
                         << ":" << proposedKinetic << ")"
                         << " <-- "
                         << acceptedHamiltonian
                         << "(" << fAcceptedPotential
                         << ":" << initialKinetic << ")"
                         << std::endl;
            // We're keeping a new step.
            for (int i=0; i<fProposed.size(); ++i) {
                fAccepted[i] = fProposed[i];
                fAcceptedMomentum[i] = fProposedMomentum[i];
            }
            fAcceptedPotential = fProposedPotential;
            // Update the current acceptance.
            fCurrentAcceptance = (fCurrentAcceptance*4999.0 + 1.0)/5000.0;
        }

        // Update the estimate of the central value.  This is simply the
        // position of the point with the minimum potential.  This is done
        // after updating the covariance since it's a convenient place seen
        // after every step.
        if (fAcceptedPotential < fCentralPotential) {
            MoveCentralPoint(fAccepted,fAcceptedPotential);
        }

        // Save the information to the output tree.
        if (save) SaveStep();
        return true;

    }

private:

    /// A wrapper to call the user LogLikelihood.  This should be used
    /// internally since it will count the number of calls. NOTICE THIS IS THE
    /// OPPOSITE OF THE LIKELIHOOD.  This is used because the HMC is mostly
    /// written in terms of a potential U = - log(likelihood).  In this
    /// formalism, the hamiltonian H = U + K where K is the analog of the
    /// kinetic energy.
    double Potential(const Vector& point) {
        ++fPotentialCount;
        return - fLogLikelihood(point);
    }

    /// Calculate the gradient of Potential using finite differences.
    void FiniteDifferenceGradient(Vector& grad, const Vector& point) {
        if (grad.size() != point.size()) {
            HMC_DEBUG(1) << "Resizing gradient to " << point.size()
                         << std::endl;
            grad.resize(point.size());
        }

        HMC_DEBUG(5) << "Finite";

        // FIXME/WARNING!!!! This is a very simple estimate right now!!!
        Vector work(point.size());
        for (int i=0; i<point.size(); ++i) {
            for (int j=0; j<point.size(); ++j) work[j] = point[j];
            // FIXME/WARNING!!! The step for each dimension should be based on
            // the curvature of the function in that dimension!!!  The
            // curvature can be estimated on the fly and du should probably be
            // different for each dimension.
            double du = 0.01;
            work[i] -= du;
            double u1 = Potential(work);
            work[i] += 2.0*du;
            double u2 = Potential(work);
            grad[i] = 0.5*(u2-u1)/du;
            HMC_DEBUG(6) << " " << grad[i];
            HMC_DEBUG(7) << "(" << u1 << ":" << u2 << ")";
        }
        HMC_DEBUG(5) << std::endl;
    }

    /// Use the running estimate of the covariance to estimate the gradient!
    void CovariantGradient(Vector& grad, const Vector& point) {
        for (int i=0; i<point.size(); ++i) {
            grad[i] = 0.0;
            for (int j=0; j<point.size(); ++j) {
                grad[i] += fEstimatedError(i,j)*(point[j]-fAveragePoint[j]);
            }
        }
    }

    /// Calculate the gradient of the potential.  This is a generalized
    /// wrapper for the optional gradient and will calculate the gradient
    /// using finite differences (if required).  It can also make estimates of
    /// the gradient using "tricks".  The type of gradient can be forced using
    /// the optional parameter: 0) use the prefered method to calculate the
    /// gradient, where "preferred" means *my* choice, 1) use a direct
    /// calculation (either the user gradient, or the finite differences), 2)
    /// use an approximate, and fast, gradient, 3) force the use of finite
    /// differences, 4) force the use of the user gradient, and 5) force tne
    /// gradient to zero. Note: The return value is for debugging, and should
    /// be ignored (it's not part of the "API").
    int PotentialGradient(Vector& grad, const Vector& point, int type) {
        Vector tgrad(grad);
        ++fPotentialGradientCount;
        switch (type) {
        default:
        case 0:
            // This is the default gradient calculation.  It will generally be
            // the same as one of the other options, but is handled
            // explicitly.  This needs to take to opposite of the user
            // gradient because the user gradient is the for log(likelihood),
            // while the gradient we need is for -log(likelihood).
            if (fUserGradient(grad,point)) {
                HMC_DEBUG(5) << "User @";
                for (int i=0; i<point.size(); ++i) {
                    HMC_DEBUG(5) << " " << point[i];
                }
                HMC_DEBUG(5) << std::endl;
                HMC_DEBUG(5) << "User  ";
                for (int i=0; i<grad.size(); ++i) {
                    grad[i] = -grad[i];
                    HMC_DEBUG(5) << " " << grad[i];
                }
                HMC_DEBUG(5) << std::endl;
                type = 4;
                break;
            }
            FiniteDifferenceGradient(grad,point);
            type = 3;
            break;
        case 1:
            // Return the best estimate of the true gradient.  This uses the
            // users gradient if it is available, or the finite difference
            // gradient if it's not available.
            if (fUserGradient(grad,point)) {
                for (int i=0; i<grad.size(); ++i) grad[i] = -grad[i];
                type = 4;
                break;
            }
            FiniteDifferenceGradient(grad,point);
            type = 3;
            break;
        case 2:
            // Use the fast, approximate gradient.
            CovariantGradient(grad,point);
            type = 2;
            break;
        case 3:
            // Force the use of a finite difference calculation, even if the
            // user gradient is available.
            FiniteDifferenceGradient(grad,point);
            type = 3;
            break;
        case 4:
            // Force the use of the users gradient or crash.
            if (!fUserGradient(grad,point)) throw;
            for (int i=0; i<grad.size(); ++i) grad[i] = -grad[i];
            type = 4;
            break;
        case 5:
            // Force a flat gradient.
            for (int i=0; i<grad.size(); ++i) grad[i] = 0.0;
            type = 5;
            break;
        }
        return type;
    }

    /// Calculate the analog of the kinetic energy from a momentum vector.
    double KineticEnergy(const Vector& momentum) {
        double ke = 0.0;
        for (int i=0; i<momentum.size(); ++i) {
            double p = momentum[i];
            ke += p*p/2.0;
        }
        return ke;
    }

    /// Propose a new momentum.  This implements a generalized HMC where the
    /// momentum can be correlated from step to step.  This closely resembles
    /// Langevin step.  It is controlled by the parameter Alpha which says how
    /// much correlation there is with the previous momentum.  If Alpha is
    /// 0.0, then this is a normal HMC step with no correlation between the
    /// new and old momentum.  If Alpha is greater than one, then the new
    /// momentum is assigned to the old momentum (but with a damping factor of
    /// Alpha).  This can be useful during the burn-in phase since this turns
    /// the HMC into a version of minimization using steepest decent.  The
    /// parameter Alpha is set using the SetAlpha() method.
    void ProposeMomentum(Vector &pNew, const Vector& momentum) {
        if (fAlpha >= 1.0) {
            // An alpha greater or equal to 1.0 means that the original
            // momentum should be returned (damped by the factor fAlpha).
            fAlpha = std::max(1.0,fAlpha);
            for (int i=0; i<momentum.size(); ++i) {
                pNew[i] = momentum[i]/fAlpha;
            }
            return;
        }
        // Randomize the momentum proposal.
        if (fAlpha < 0.0) fAlpha = 0.0;
        for (int i=0; i<momentum.size(); ++i) {
            pNew[i] = fAlpha*momentum[i]
                + std::sqrt(1.0-fAlpha*fAlpha)*gRandom->Gaus(0.0,1.0);
        }
    }

    /// Do a "Leap Frog" step.  The leap frog step will be made up of a total
    /// of "steps", each of length "epsilon".  This only evolves position and
    /// momentum.  Specifically, it does not resample the momentum.  The total
    /// step length will be "epsilon*steps".  On return, "qNew" will hold the
    /// new position, and "pNew" will hold the new momentum.  The optional
    /// "type" parameter controls how the gradient is calculated.  See the
    /// Gradient method for documentation.
    ///
    /// If the input value of "steps" is zero, then this will take one step of
    /// length "epsilon".
    int LeapFrog(Vector &qNew, Vector& pNew,
                  const Vector &position,
                  double epsilon, int steps,
                  int type=0) {
        qNew = position;
        Vector momentum = pNew;
        Vector grad(momentum.size());
        int leapStatus = 1;

        // First a SPECIAL cheaters short cut.  If "steps" is not positive, do
        // ONE full step and return.  This is NOT normally done.  To use it,
        // set "steps" to zero.  The input value of "epsilon" will be the
        // total step length.  NOTE: Numerically, this is not stable, and the
        // resulting step is not reversible, but can be useful (e.g. during
        // the initial burn-in).  NOTE 2: For HMC to work, the step has to be
        // reversible, so using this breaks the algorithm.
        if (steps<1) {
#ifdef INCLUDE_GRADIENT_IN_FORCED_STEP
            // Find the gradient.
            PotentialGradient(grad,qNew,type);
            std::cout << "forced step" << std::endl;
            for (std::size_t j=0; j<position.size(); ++j) {
                 pNew[j] = pNew[j] - epsilon*grad[j];
            }
#endif
            for (std::size_t j=0; j<position.size(); ++j) {
                qNew[j] = qNew[j] + epsilon*(momentum[j]+pNew[j])/2.0;
            }
            return leapStatus;
        }

        // Find the gradient at the first point.  This was already calculated,
        // so it could have been cached!
        PotentialGradient(grad,qNew,type);

        // Do the first half step for pNew (the momentum)
        for (std::size_t j=0; j<position.size(); ++j) {
            pNew[j] = pNew[j] - epsilon*grad[j]/2.0;
        }

        // Do everything but the step for qNew and last half step for pNew.
        for (int i = 0; i<steps-1; ++i) {
            for (std::size_t j=0; j<position.size(); ++j) {
                qNew[j] = qNew[j] + epsilon*pNew[j];
            }
            PotentialGradient(grad,qNew,type);
            for (std::size_t j=0; j<position.size(); ++j) {
                pNew[j] = pNew[j] - epsilon*grad[j];
            }
            // Check to see if the direction has reversed (reached the outer
            // limit of the orbit).
            double inner = 0.0;
            for (std::size_t j=0; j<position.size(); ++j) {
                inner += pNew[j]*momentum[j];
            }
            if (inner >= 0.0) continue;
            leapStatus = 2;
        }
        // Do the last step for qNew
        for (std::size_t j=0; j<position.size(); ++j) {
            qNew[j] = qNew[j] + epsilon*pNew[j];
        }
        // Do the last half step for pNew
        PotentialGradient(grad,qNew,type);
        for (std::size_t j=0; j<position.size(); ++j) {
            pNew[j] = pNew[j] - epsilon*grad[j]/2.0;
        }

        return leapStatus;
    }

    /// Move the position of the central point, and update the estimated
    /// covariance to account for the change.
    void MoveCentralPoint(const Vector& newCenter, Parameter newPotential) {
        std::copy(newCenter.begin(), newCenter.end(), fCentralPoint.begin());
        fCentralPotential = newPotential;
    }

    /// Keep a running tabulation of the estimated covariance and central
    /// value of the posterior.  The resulting covariance can be used to help
    /// control the evolution of the chain.  The accepted is the last step
    /// that was accepted, the proposed in the next step that is being
    /// checked.  The proposed step may eventually become the accepted step.
    void UpdateCovariance(const Vector& accepted, Parameter acceptedPotential,
                          const Vector& proposed, Parameter proposedPotential) {
        ++fStepsSinceUpdate;
        --fStepsRemaining;

        // Update the running average of the likelihood weighted position.
        for (std::size_t i=0; i<accepted.size(); ++i) {
            double v = fAveragePoint[i];
            v *= fAveragePointTrials;
            v += accepted[i];
            v /= fAveragePointTrials + 1.0;
            fAveragePoint[i] = v;
        }
        fAveragePointTrials = std::min(fCovarianceWindow,
                                       fAveragePointTrials+1.0);

        for (std::size_t i=0; i<accepted.size(); ++i) {
            for (std::size_t j=0; j<i+1; ++j) {
                double v = fEXXT(i,j);
                v *= fCovarianceTrials;
                v += accepted[i]*accepted[j];
                v /= fCovarianceTrials + 1.0;
                fEXXT(i,j) = fEXXT(j,i) = v;
                fEstimatedCovariance(i,j) = fEstimatedCovariance(j,i)
                    = fEXXT(i,j) - fAveragePoint[i]*fAveragePoint[j];
            }
        }
        fCovarianceTrials = std::min(fCovarianceWindow,
                                     fCovarianceTrials+1.0);

    }

    /// Update the parameters estimated from the covariance.  The covariance
    /// is updated with every step, but the error matrix, eigenvalues and
    /// other parameters are too expensive to recalculate.  This method can be
    /// called with every step, but only updates the calculated parameters
    /// when the covariance has changed significantly.  It will also force an
    /// update periodically.
    void UpdateErrorMatrix() {
        if (!fLeapFrogSteps) return;
        if (fCovarianceTrials < 2*fEstimatedCovariance.GetNrows()) return;

        // Find the current trace...
        fCurrentCovarianceTrace = 0.0;
        for (int i=0; i<fEstimatedCovariance.GetNrows(); ++i) {
            fCurrentCovarianceTrace += std::abs(fEstimatedCovariance(i,i));
        }
        double change = std::abs(fCurrentCovarianceTrace
                                 -fEstimatedCovarianceTrace);

        bool doIt = false;
        if (fStepsRemaining < 0) doIt = true;
        if (fStepsSinceUpdate > 2.0*fCentralPoint.size()
            && change > 0.01*fEstimatedCovarianceTrace) doIt = true;
        if (!doIt) return;

        HMC_DEBUG(1) << "Update(" << fStepCount << ") "
                     << " Old Trace: " << fEstimatedCovarianceTrace
                     << " New Trace: " << fCurrentCovarianceTrace
                     << " Change: " << change/fCurrentCovarianceTrace
                     << std::endl;

        // Check if the average point is better than the current central
        // point.  If it is, then the average point becomes the central point.
        double aPot = Potential(fAveragePoint);
        HMC_DEBUG(1) << "   Central Potential " << fCentralPotential
                     << " Average Potential " << aPot
                     << " Accepted " << GetAcceptanceRate()
                     << std::endl;
        if (aPot < fCentralPotential) {
            // NOTE: THIS IS BEING SET TO THE AVERAGE POINT, NOT THE MINIMUM
            // OF THE POTENTIAL.  THE MINIMUM OFTEN ISN'T AT THE AVERAGE.
            MoveCentralPoint(fAveragePoint,aPot);
        }

        double cDist = 0.0;
        for (int i=0; i<fCentralPoint.size(); ++i) {
            cDist += fCentralPoint[i]*fCentralPoint[i];
        }
        cDist = std::sqrt(cDist);
        HMC_DEBUG(1) << "   Central Distance --"
                     << " " << cDist
                     << " " << fCentralPotential
                     << ":";
        for (int i=0; i<fCentralPoint.size(); ++i) {
            if (i > 2) {
                HMC_DEBUG(1) << "..";
                break;
            }
            HMC_DEBUG(1) << " " << fCentralPoint[i];
        }
        HMC_DEBUG(1) << std::endl;

        fStepsRemaining = 2*fCentralPoint.size() + fStepCount;
        fStepsSinceUpdate = 0;

        // Make sure the estimated covariance is positive definite.
        TVectorD eigenValues;
        double maxScale = 0.0;
        double minScale = 1E+20;
        do {
            fEstimatedCovariance.EigenVectors(eigenValues);
            bool positiveDefinite = true;
            double avgEigen = 0.0;
            double maxEigen = 0.0;
            double minEigen = 1E+10;
            for (std::size_t i = 0; i<fEstimatedCovariance.GetNrows(); ++i) {
                double eigen = eigenValues(i);
                avgEigen += std::abs(eigen);
                if (maxEigen < eigen) maxEigen = eigen;
                if (maxScale < std::abs(eigen)) maxScale = std::abs(eigen);
                if (minEigen > eigen) minEigen = eigen;
                if (minScale > std::abs(eigen)) minScale = std::abs(eigen);
                if (eigen < 0) {
                    positiveDefinite = false;
                }
            }
            if (positiveDefinite) break;
            HMC_DEBUG(1) << "Eigenvectors --"
                         << " average: "
                         << avgEigen/fEstimatedCovariance.GetNrows()
                         << " min: " << minEigen
                         << " max: " << maxEigen
                         << " scale: " << std::sqrt(minScale)
                         << " -- " << std::sqrt(maxScale)
                         << std::endl;
            HMC_DEBUG(1) << "Covariance is not positive definite" << std::endl;
            for (std::size_t i = 0; i<fEstimatedCovariance.GetNrows(); ++i) {
                // Make sure the variance for each element isn't too small.
                double r = fEstimatedCovarianceTrace*1E-6;
                r /= fEstimatedCovariance.GetNrows();
                r = std::abs(r);
                if (fEstimatedCovariance(i,i) < r) {
                    fEstimatedCovariance(i,i) = r;
                }
                // Make sure the correlations aren't too large.
                for (std::size_t j = i+1;
                     j<fEstimatedCovariance.GetNrows(); ++j) {
                    fEstimatedCovariance(i,j) = 0.0;
                    fEstimatedCovariance(j,i) = fEstimatedCovariance(i,j);
                }
            }
        } while (true);

        // Recalculate the trace of the estimated covariance.  This could be
        // copied from old value of fCurrentCovarianceTrace, but is redone in
        // case the adjustments to keep the matrix positive definite changed
        // the trace.
        fCurrentCovarianceTrace = 0.0;
        for (int i=0; i<fEstimatedCovariance.GetNrows(); ++i) {
            fCurrentCovarianceTrace += std::abs(fEstimatedCovariance(i,i));
        }
        fEstimatedCovarianceTrace = fCurrentCovarianceTrace;

        // Estimate the scale of the largest dimension.
        maxScale = std::sqrt(maxScale);
        if (maxScale < 0.1) maxScale = 0.1;

        // Estimate the scale of the smallest dimension.
        minScale = std::sqrt(minScale);
        if (minScale < 0.01) minScale = 0.01;

        // Estimate the circumference of the largest "great circle"
        fEstimatedOrbitLength = 2.0*3.14*maxScale;

        // If the value of fMeanEpsilon is positive, then that means the new
        // value should be estimated based on what is known about the
        // covariance.
        if (fMeanEpsilon > 0) {
            fMeanEpsilon = 0.2*maxScale;
            if (fMeanEpsilon > 0.5*minScale) fMeanEpsilon = 0.5*minScale;
            if (fMeanEpsilon < 0.05*maxScale) fMeanEpsilon = 0.05*maxScale;
        }

        if (fLeapFrogSteps>0) {
            double targetLength = 0.4*fEstimatedOrbitLength;
            fLeapFrogSteps = targetLength/std::abs(fMeanEpsilon);
            fLeapFrogSteps = 2*(fLeapFrogSteps/2 + 1);
            if (fLeapFrogSteps > 3*fCentralPoint.size()) {
                fLeapFrogSteps = 3*fCentralPoint.size();
            }
            if (fMeanEpsilon>0) fMeanEpsilon = targetLength/fLeapFrogSteps;
        }

        fEstimatedError = fEstimatedCovariance;
        fEstimatedError.Invert();

        HMC_DEBUG(1) << "   Orbit: " << fEstimatedOrbitLength
                     << " Scale: [" << minScale
                     << ", " << maxScale << "]"
                     << " Epsilon: " << fMeanEpsilon
                     << " Steps: " << std::abs(fLeapFrogSteps)
                     << std::endl;
    }

    /// If possible, save the step.
    void SaveStep() {if (fTree) fTree->Fill();}

    /// A TTree to save the accepted points.
    TTree* fTree;

    /// The loglikelihood being explored.
    LogLikelihood fLogLikelihood;

    /// The users gradient (may be a dummy function).
    UserGradient fUserGradient;

    /// A count of the total steps taken.
    int fStepCount;

    /// A count of the total number of calls to the Potential method.
    int fPotentialCount;

    /// A count of the total number of calls to the PotentialGradient method.
    int fPotentialGradientCount;

    /// The current acceptance of the recent history of the chain.
    double fCurrentAcceptance;

    /// The target acceptance for the chain.
    double fTargetAcceptance;

    /// The approximate mean value for the epsilon
    double fMeanEpsilon;

    /// The number of leap frog steps to take at each iteration
    int fLeapFrogSteps;

    /// The recent distance traveled before the path reverses.
    double fReversalLen;

    /// The correlation of the momentum from step to step.  If this is zero,
    /// then you get the original HMC.  If it's large, it become similar to an
    /// Langevin MC.
    double fAlpha;

    /// The last accepted point.  This will be the same as the proposed point
    /// if the last step was accepted.
    Vector fAccepted;

    /// The last accepted momentum.
    Vector fAcceptedMomentum;

    /// The likelihood at the last accepted piont.
    double fAcceptedPotential;

    /// The proposed point for the most recent step (This may be the same as
    /// the accepted point).
    Vector fProposed;

    /// The last proposed momentum.
    Vector fProposedMomentum;

    /// The likelihood at the last proposed point.
    double fProposedPotential;

    /// The current (running) estimate of the minimum potential value.
    double fCentralPotential;

    /// The current (running) estimate for the central value of each parameter.
    Vector fCentralPoint;

    /// The current (running) average of the posterior.  This is an estimate of
    /// the central point.
    Vector fAveragePoint;

    // The trials being used for the current average point.  This
    // will be a value between one and fCovarianceWindow.
    double fAveragePointTrials;

    /// The estimated covariance of the posterior.
    TMatrixD fEstimatedCovariance;

    // The sum of the "square" terms that are E[X X^T]
    TMatrixD fEXXT;

    // The trials being used for the current estimated covariance.  This
    // will be a value between one and fCovarianceWindow.
    double fCovarianceTrials;

    /// The estimated error matrix of the posterior.
    TMatrixD fEstimatedError;

    /// The estimated orbit length for the posterior.
    double fEstimatedOrbitLength;

    // The trace of the estimated covariance used to calculate the estimated
    // error matrix.
    double fEstimatedCovarianceTrace;

    // The trace of the current estimated covariance
    double fCurrentCovarianceTrace;

    // The number of steps until an update of the estimated error will be
    // forced.
    int fStepsRemaining;

    // the number of steps since the last update.
    int fStepsSinceUpdate;

    // The window to average the covariance and central value over.  The
    // covariance window should be several times the autocorrelation period.
    // It should usually be very large so that all to points used to probe the
    // posterior are used.  Setting it to a smaller value means that the local
    // covariance of the posterior will be used.  This might be important if
    // the posterior is extremely non-Gaussian (e.g. it's a "banana
    // posterior").
    double fCovarianceWindow;
};
#endif
