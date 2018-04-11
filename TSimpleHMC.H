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

/// WIP!!! A class to run a Hamiltonian (or Hybrid) Monte Carlo. The
/// UserLikelihood template argument must be a class (or struct) which
/// provides a method declared as:
///
///\code
/// struct ExampleLogLikelihood {
///    double operator() (const Vector& point);
/// }
///\endcode
///
/// The optional UserGradient will should calculate the gradient of the user
/// likelihood, and return true if the gradient is calculated.  If provided,
/// the template argument must be a class (or struct) which provides a method
/// declared as:
///
///\code
/// struct ExampleGradient {
///    bool operator() (const Vector& point, Vector& grad);
/// }
///\endcode
///
template <typename UserParameter,
          typename OptionalGradient = SimpleHMCInvalidGradient>
class TSimpleHMC {
public:
    
    TSimpleHMC(TTree* tree = NULL, bool saveStep = false)
        : fTree(tree), fPotentialCount(0) {
        if (fTree) {
            HMC_DEBUG(0) << "TSimpleHMC: Adding branches to "
                          << fTree->GetName()
                          << std::endl;
            fTree->Branch("LogLikelihood",&fAcceptedPotential);
            fTree->Branch("Accepted",&fAccepted);
        }
    }
    
    /// Make the user likelihood class available as TSimpleHMC::LogLikelihood.
    typedef UserParameter LogLikelihood;

    /// Make optional gradient class available as TSimpleHMC::UserGradient.
    /// This is mostly here for debugging purposes.
    typedef OptionalGradient UserGradient;

    /// The loglikelihood being explored.
    LogLikelihood fLogLikelihood;

    /// The users gradient (may be a dummy function).
    UserGradient fUserGradient;

    /// Get a reference to the likelihood calculation object.  The
    /// LogLikelihood object is constructed by the TSimpleHMC template and
    /// is the class that you handed to the template.  You have write access
    /// to that object with this method, and can access any of your declared
    /// methods.  Because this is a template, TSimpleHMC::LogLikelihood is
    /// actually just a typedef for your class.
    LogLikelihood& GetLogLikelihood() {return fLogLikelihood;}

    /// A wrapper to call the user LogLikelihood.  This should be used
    /// internally since it will count the number of calls. NOTICE THIS IS THE
    /// OPPOSITE OF THE LIKELIHOOD.
    double Potential(const Vector& point) {
        ++fPotentialCount;
        return - fLogLikelihood(point);
    }

    /// Get a count of the total number of calls to the Potential method.
    int GetPotentialCount() const {
        return fPotentialCount;
    }
    
    /// A count of the total number of calls to the Potential method.
    int fPotentialCount;
    
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
    
    /// Calculate the gradient.  This is a generalized wrapper for the
    /// optional gradient and will calculate the gradient using finite
    /// differences (if required).  It can also make estimates of the gradient
    /// using "tricks".  The type of gradient can be forced using the optional
    /// parameter: 0) use the prefered method to calculate the gradient, where
    /// "preferred" means *my* choice, 1) use a direct calculation (either the
    /// user gradient, or the finite differences), 2) use an approximate, and
    /// fast, gradient, 3) force the use of finite differences, 4) force the
    /// use of the user gradient. Note: The return value is for debugging, and
    /// should be ignored (it's not part of the "API").
    int Gradient(Vector& grad, const Vector& point, int type=0) {
        switch (type) {
        case 0:
        default:
            if (fUserGradient(grad,point)) {
                type = 4;
                break;
            }
            FiniteDifferenceGradient(grad,point);
            type = 3;
            break;
        }
        return type;
    }

    double KineticEnergy(const Vector& momentum) {
        double ke = 0.0;
        for (int i=0; i<momentum.size(); ++i) {
            double p = momentum[i];
            ke += p*p/2.0;
        }
        return ke;
    }
    
    /// Propose a new momentum.
    void ProposeMomentum(Vector &pNew, const Vector& momentum) {
        double alpha = 0.9;
        if (alpha > 0.9999) return;
        for (int i=0; i<momentum.size(); ++i) {
            pNew[i] = alpha*momentum[i]
                + std::sqrt(1.0-alpha*alpha)*gRandom->Gaus(0.0,1.0);
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
        Gradient(grad,qNew,type);

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
            Gradient(grad,qNew,type);
            for (std::size_t j=0; j<position.size(); ++j) {
                pNew[j] = pNew[j] - epsilon*grad[j];
            }
        }
        // Do the last step for qNew
        for (std::size_t j=0; j<position.size(); ++j) {
            qNew[j] = qNew[j] + epsilon*pNew[j];
        }
        // Do the last half step for pNew
        Gradient(grad,qNew,type);
        for (std::size_t j=0; j<position.size(); ++j) {
            pNew[j] = pNew[j] - epsilon*grad[j]/2.0;
        }
    }

    /// Set the starting point for the mcmc.  If the optional argument is
    /// true, then the point will be saved to the output.
    void Start(Vector start, bool save=true) {
        fProposed.resize(start.size());
        std::copy(start.begin(), start.end(), fProposed.begin());
        
        fProposedMomentum.resize(start.size());

        fAccepted.resize(start.size());
        std::copy(start.begin(), start.end(), fAccepted.begin());

        fAcceptedMomentum.resize(start.size());

        fProposedPotential = Potential(fProposed);
        fAcceptedPotential = fProposedPotential;

        fMeanEpsilon = 0.1;
        
        if (save) SaveStep();
    }
    
    /// Take a step.  This returns true if a new point has been accepted, and
    /// false if we stay at the old point.  If save is true, then the points
    /// are saved to the output.
    bool Step(bool save=true) {
        if (fProposed.empty() || fProposedMomentum.empty()
            || fAccepted.empty() || fAcceptedMomentum.empty()) {
            HMC_ERROR << "Must initialize starting point" << std::endl;
            throw;
        }

        ProposeMomentum(fProposedMomentum,fAcceptedMomentum);
        
        fAcceptedMomentum = fProposedMomentum;
        fAcceptedKinetic = KineticEnergy(fAcceptedMomentum);

        double epsilon = gRandom->Uniform(0.95*fMeanEpsilon,1.05*fMeanEpsilon);
        LeapFrog(fProposed,fProposedMomentum,fAccepted,epsilon,1);
        
        // Find the proposed kinetic and potential energy
        fProposedKinetic = KineticEnergy(fProposedMomentum);
        fProposedPotential = Potential(fProposed);

        // Find the change in energy between the proposed and accepted states.
        double proposedHamiltonian = fProposedPotential + fProposedKinetic;
        double acceptedHamiltonian = fAcceptedPotential + fAcceptedKinetic;
        double delta = proposedHamiltonian - acceptedHamiltonian;
        if (delta > 0.0 ) {
            // The proposed hamiltonian is more than the previously accepted
            // hamiltonian, so see if it should be rejected.
            double trial = - std::log(gRandom->Uniform());
            if (delta > trial) {
                // The new step should be rejected, so save the old step.
                // This depends on IEEE error handling so that - std::log(0.0)
                // is inf which is always more than delta.
                if (save) SaveStep();
                fMeanEpsilon *= 0.99;
                return false;
            }
        }
        
        fMeanEpsilon *= 1.01;
        // We're keeping a new step.
        for (int i=0; i<fProposed.size(); ++i) {
            fAccepted[i] = fProposed[i];
            fAcceptedMomentum[i] = fProposedMomentum[i];
        }
        fAcceptedPotential = fProposedPotential;

        // Save the information to the output tree.
        if (save) SaveStep();
        return true;
    }

    /// If possible, save the step.
    void SaveStep() {if (fTree) fTree->Fill();}
 
    /// A TTree to save the accepted points.
    TTree* fTree;

    // The approximate mean value for the epsilon
    double fMeanEpsilon;
    
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
};
#endif
