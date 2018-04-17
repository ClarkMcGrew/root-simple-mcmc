#ifndef TSimpleHMC_H_SEEN
#define TSimpleHMC_H_SEEN

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

#include <TRandom.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TVectorD.h>

// Make a typedef for the type used for the function parameter.  The
// parameter type will usually be a double, but for some problems that might
// not be a good idea (e.g. uses too much memory, or it's not needed).  This
// makes it easy to change between the two).
typedef double Parameter;

// Make a typedef for a vector of the parameters.  This is for the same reason
// as "Parameter".
typedef std::vector<Parameter> Vector;

// Define the amount of debugging when running the chain.
#ifndef HMC_DEBUG_LEVEL
#define HMC_DEBUG_LEVEL 2
#endif

// A macro to make outputting debug messages easy.  This can be used as
//
// HMC_DEBUG(1) << "This is debug output" << std::endl;
//
// The message will be printed if the value is less that HMC_DEBUG_LEVEL.
#ifndef HMC_DEBUG
#define HMC_DEBUG(level) if (level <= (HMC_DEBUG_LEVEL)) std::cout
#endif

// A macro to make outputting error messages easy.  It can be used like
// "std::cout".
#define HMC_ERROR (std::cout <<__FILE__<<":: " << __LINE__ << ": " )

namespace {
    // A place holder to flag that the user didn't provide a gradient.  It
    // will return "false" to indicate that the gradient was not calculated.
    struct SimpleHMCInvalidGradient {
        bool operator() (Vector&, const Vector&) {return false;}
    };
}

/// A class to run a Hamiltonian (or Hybrid) Monte Carlo.  This is implemented
/// based on Chapter 5 of the Handbook of Markov Chain Monte Carlo -- MCMC
/// Using Hamiltonian Dynamics by Radford M. Neal..  The UserLikelihood
/// template argument must be a class (or struct) which provides a method
/// declared as:
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
/// log(likelihood), and return true if the gradient is calculated.  If
/// provided, the template argument must be a class (or struct) which provides
/// a method declared as:
///
///\code
/// struct ExampleGradient {
///    bool operator() (Vector& grad, const Vector& point);
/// }
///\endcode
///
template <typename UserParameter,
          typename OptionalGradient = SimpleHMCInvalidGradient>
class TSimpleHMC {
public:
    /// Make the user likelihood class available as TSimpleHMC::LogLikelihood.
    typedef UserParameter LogLikelihood;

    /// Make optional gradient class available as TSimpleHMC::UserGradient.
    /// This is mostly here for debugging purposes.
    typedef OptionalGradient UserGradient;

    TSimpleHMC(TTree* tree = NULL, bool saveStep = false)
        : fTree(tree), fStepCount(0), fPotentialCount(0),
          fLeapFrogSteps(100), fAlpha(0.0),
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

    /// Set the starting point for the HMC chain.  If the optional argument is
    /// true, then the point will be saved to the output.  This must be called
    /// in the user code to initialize the HMC chain.
    void Start(Vector start, bool save=true) {
        fStepCount = 0;

        // Set the vector size.
        fProposed.resize(start.size());
        fProposedMomentum.resize(start.size());
        fAccepted.resize(start.size());
        fAcceptedMomentum.resize(start.size());

        // Use the starting point as the first proposal.
        std::copy(start.begin(), start.end(), fProposed.begin());
        fProposedPotential = Potential(fProposed);

        // The first step is always accepted.
        std::copy(start.begin(), start.end(), fAccepted.begin());
        fAcceptedPotential = fProposedPotential;

        // Set the value of the step size.  The step size will evolve based on
        // the walk through the log(likeihood).
        fMeanEpsilon = 0.05;

        // Set the target value for the step acceptance probability.  The
        // value of 65% comes from the HMC literature.
        fTargetAcceptance = 0.65;
        fCurrentAcceptance = 0.50;

        // Set up a running calculation of the covariance.  This starts the
        // calculation with all the variances equal to 1, and no correlations
        // between the dimensions.  The initial guess is given a weight of 10
        // (each added point has a weight of 1).
        fEstimatedCovariance.ResizeTo(start.size(), start.size());
        fEstimatedError.ResizeTo(start.size(), start.size());
        fCentralPoint.resize(start.size());
        fCovarianceTrials = 10;
        for (int i=0; i<start.size(); ++i) {
            fCentralPoint[i] = start[i];
            for (int j=0; j<start.size(); ++j) {
                if (i == j) fEstimatedCovariance(i,j) = 1.0;
                else fEstimatedCovariance(i,j) = 0.0;
            }
        }
        fEstimatedError = fEstimatedCovariance;
        fEstimatedError.Invert();
        fEstimatedCovarianceTrace = start.size();
        fStepsRemaining = start.size();

        if (save) SaveStep();
    }
    
    /// Take a step.  This returns true if a new point has been accepted, and
    /// false if we stay at the old point.  If save is true, then the points
    /// are saved to the output.  This should be called repeatedly by the user
    /// code.
    bool Step(bool save=true) {
        if (fProposed.empty() || fProposedMomentum.empty()
            || fAccepted.empty() || fAcceptedMomentum.empty()) {
            HMC_ERROR << "Must initialize starting point" << std::endl;
            throw;
        }

        ++fStepCount;
        
        ProposeMomentum(fProposedMomentum,fAcceptedMomentum);

        //Tentatively accept the starting point of the step.
        fAcceptedMomentum = fProposedMomentum;
        fAcceptedKinetic = KineticEnergy(fAcceptedMomentum);

        double epsilon = gRandom->Uniform(0.9*fMeanEpsilon,1.1*fMeanEpsilon);
        LeapFrog(fProposed,fProposedMomentum,fAccepted,
                 epsilon,std::abs(fLeapFrogSteps));
        
        // Find the proposed kinetic and potential energy
        fProposedKinetic = KineticEnergy(fProposedMomentum);
        fProposedPotential = Potential(fProposed);

        // Find the change in energy between the proposed and accepted states.
        // If the leapfrog step was done with perfect accuracy, we would have
        // energy conservation, and the hamiltonian would remain constant
        // (i.e. delta = 0.0).
        double proposedHamiltonian = fProposedPotential + fProposedKinetic;
        double acceptedHamiltonian = fAcceptedPotential + fAcceptedKinetic;
        double delta = proposedHamiltonian - acceptedHamiltonian;

        double trial = - std::log(gRandom->Uniform());
        if (delta > trial) {
            // The proposed hamiltonian is more than the previously accepted
            // hamiltonian, so see if it should be rejected.  This depends on
            // IEEE error handling so that - std::log(0.0) is inf which is
            // always more than delta.  If the step is rejected, then all that
            // needs to be done is to update the current acceptance.
            fCurrentAcceptance = (fCurrentAcceptance*999.0)/1000.0;
        }
        else {
            // We're keeping a new step.
            for (int i=0; i<fProposed.size(); ++i) {
                fAccepted[i] = fProposed[i];
                fAcceptedMomentum[i] = fProposedMomentum[i];
            }
            fAcceptedPotential = fProposedPotential;
            // Update the current acceptance.
            fCurrentAcceptance = (fCurrentAcceptance*999.0 + 1.0)/1000.0;
        }

        // Update the best guess at the best step length.
        fMeanEpsilon *= 1.0
            + 0.005*(fCurrentAcceptance-fTargetAcceptance)/fTargetAcceptance;

        // Update the running estimate of the covariance.
        UpdateCovariance(fAccepted);
        UpdateErrorMatrix();
        
        // Save the information to the output tree.
        if (save) SaveStep();
        return true;
    }

    /// Get a count of the total number of calls to the Potential method.
    int GetPotentialCount() const {
        return fPotentialCount;
    }

    /// Set the correlation between the last accepted momentum and the new
    /// proposed momentum.  See the ProposeMomentum() method for a description
    /// of Alpha.
    void SetAlpha(double a) {fAlpha = a;}

    /// Get the current value of alpha.
    double GetAlpha() const {return fAlpha;}

    /// Override the automatically calculated number of leapfrog steps.
    void SetLeapFrog(int i) {fLeapFrogSteps = -i;}
    
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
        }
    }

    /// Use the running estimate of the covariance to estimate the gradient!
    /// This is here for testing, but doesn't give the right answer, so DO NOT
    /// USE IT!.
    void CovariantGradient(Vector& grad, const Vector& point) {
        for (int i=0; i<point.size(); ++i) {
            grad[i] = 0.0;
            for (int j=0; j<point.size(); ++j) {
                grad[i] += fEstimatedError(i,j)*(point[j]-fCentralPoint[j]);
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
    /// differences, 4) force the use of the user gradient. Note: The return
    /// value is for debugging, and should be ignored (it's not part of the
    /// "API").
    int PotentialGradient(Vector& grad, const Vector& point, int type=0) {
        switch (type) {
        default:
        case 0:
            // This is the default gradient calculation.  It will generally be
            // the same as one of the other options, but is handled
            // explicitly.  This needs to take to opposite of the user
            // gradient because the user gradient is the for log(likelihood),
            // while the gradient we need is for -log(likelihood).
            if (fUserGradient(grad,point)) {
                for (int i=0; i<grad.size(); ++i) {
                    grad[i] = -grad[i];
                }
                type = 4;
                break;
            }
            FiniteDifferenceGradient(grad,point);
            type = 3;
            break;
        case 1:
            // This is for when we want the users gradient if it is available,
            // or the finite difference gradient.
            if (fUserGradient(grad,point)) {
                for (int i=0; i<grad.size(); ++i) grad[i] = -grad[i];
                type = 4;
                break;
            }
            FiniteDifferenceGradient(grad,point);
            type = 3;
            break;
        case 2:
            // A test form using the covariance.  This is not a good estimate
            // of the gradient.
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
            // Force the use of the users gradient, or crash.
            if (!fUserGradient(grad,point)) throw;
            for (int i=0; i<grad.size(); ++i) grad[i] = -grad[i];
            type = 4;
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
        if (fAlpha > 0.9999) {
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
    void LeapFrog(Vector &qNew, Vector& pNew, const Vector &position, 
                  double epsilon, int steps, int type=0) {
        qNew = position;
        Vector momentum = pNew;
        Vector grad(momentum.size());

        // Find the gradient at the first point.  This was already calculated,
        // so it could have been cached!
        PotentialGradient(grad,qNew,type);

        // First a SPECIAL cheaters short cut.  If steps is not positive, do
        // ONE full step and return.  Numerically, this is not particularly
        // stable, but can be useful...  This is NOT normally done.  To use
        // it, set "steps" to zero.  The input value of "epsilon" will be the
        // total step length.
        if (steps<1) {
            for (std::size_t j=0; j<position.size(); ++j) {
                pNew[j] = pNew[j] - epsilon*grad[j];
            }
            for (std::size_t j=0; j<position.size(); ++j) {
                qNew[j] = qNew[j] + epsilon*(momentum[j]+pNew[j])/2.0;
            }            
            return;
        }
        
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
    }

    /// Keep a running tabulation of the estimated covariance and central
    /// value of the posterior.  The resulting covariance can be used to help
    /// control the evolution of the chain.
    void UpdateCovariance(const Vector& accepted) {
        // Update the estimate of the central value.  This is simply a running
        // average of the position of the points in the posterior.
        for (std::size_t i=0; i<accepted.size(); ++i) {
            fCentralPoint[i] *= fCovarianceTrials;
            fCentralPoint[i] += accepted[i];
            fCentralPoint[i] /= fCovarianceTrials + 1;
        }
        
        // Update the estimate of the covariance.  This is a running
        // calculation of the covariance.
        for (std::size_t i=0; i<accepted.size(); ++i) {
            for (std::size_t j=0; j<i+1; ++j) {
                double v = fEstimatedCovariance(i,j);
                double r = (accepted[i]-fCentralPoint[i])
                    *(accepted[j]-fCentralPoint[j]);
                v *= fCovarianceTrials;
                v += r;
                v /= fCovarianceTrials + 1.0;
                if (i == j) fEstimatedCovariance(i,j) = v;
                else fEstimatedCovariance(i,j) = fEstimatedCovariance(j,i) = v;
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
        // Find the current trace...
        fCurrentCovarianceTrace = 0.0;
        for (int i=0; i<fEstimatedCovariance.GetNrows(); ++i) {
            fCurrentCovarianceTrace += fEstimatedCovariance(i,i);
        }
        double change = std::abs(fCurrentCovarianceTrace
                                 -fEstimatedCovarianceTrace);
        
        if (0 < --fStepsRemaining 
            && change < 0.01*fEstimatedCovarianceTrace) return;

        HMC_DEBUG(1) << "Update estimates -- "
                     << " Old Trace: " << fEstimatedCovarianceTrace
                     << " New Trace: " << fCurrentCovarianceTrace
                     << " Change: " << change/fCurrentCovarianceTrace
                     << std::endl;

        fStepsRemaining = fStepCount;
        fEstimatedCovarianceTrace = fCurrentCovarianceTrace;

        fEstimatedError = fEstimatedCovariance;
        fEstimatedError.Invert();

        // Estimate the circumference of the largest "great circle"
        TVectorD eigenValues;
        fEstimatedCovariance.EigenVectors(eigenValues);
        fEstimatedOrbitLength = 3.14*std::sqrt(eigenValues(0));
        
        // Estimate the scale of the smallest dimension.
        double minScale = eigenValues(fEstimatedCovariance.GetNrows()-1);
        minScale = std::sqrt(minScale);

        // Use the size of the smallest dimension to update the step
        // size.
        fMeanEpsilon = std::min(minScale,fMeanEpsilon);
        
        if (fLeapFrogSteps>0) {
            fLeapFrogSteps = 0.4*fEstimatedOrbitLength/fMeanEpsilon;
            fLeapFrogSteps = 2*(fLeapFrogSteps/2 + 1);
        }

        HMC_DEBUG(1) << "     Estimated orbit: " << fEstimatedOrbitLength
                     << " Min scale: " << minScale
                     << " Epsilon: " << fMeanEpsilon
                     << " Leapfrog steps: " << std::abs(fLeapFrogSteps)
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

    /// The current acceptance of the recent history of the chain.
    double fCurrentAcceptance;

    /// The target acceptance for the chain.
    double fTargetAcceptance;
    
    /// The approximate mean value for the epsilon
    double fMeanEpsilon;

    /// The number of leap frog steps to take at each iteration
    int fLeapFrogSteps;

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

    /// The analog to kinetic energy at the last momentum
    double fAcceptedKinetic;

    /// The proposed point for the most recent step (This may be the same as
    /// the accepted point).
    Vector fProposed;

    /// The last proposed momentum.
    Vector fProposedMomentum;
    
    /// The likelihood at the last proposed point.
    double fProposedPotential;

    /// The analog to kinetic energy at the proposed momentum
    double fProposedKinetic;

    // The current (running) estimate for the central value of each parameter.
    Vector fCentralPoint;

    /// The estimated covariance of the posterior.
    TMatrixD fEstimatedCovariance;

    /// The estimated error matrix of the posterior.
    TMatrixD fEstimatedError;

    /// The estimated orbit length for the posterior.
    double fEstimatedOrbitLength;
    
    // The trials being used for the current estimated covariance.  This
    // will be a value between one and fCovarianceWindow.
    double fCovarianceTrials;

    // The trace of the estimated covariance used to calculate the estimated
    // error matrix.
    double fEstimatedCovarianceTrace;

    // The trace of the current estimated covariance
    double fCurrentCovarianceTrace;
    
    // The number of steps until an update of the estimated error will be
    // forced.
    int fStepsRemaining;

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
