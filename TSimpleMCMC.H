#ifndef TSimpleMCMC_H_SEEN
#define TSimpleMCMC_H_SEEN

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

#include <TRandom.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TDecompChol.h>

typedef double Parameter;
typedef std::vector<Parameter> Vector;

#ifndef MCMC_DEBUG_LEVEL
#define MCMC_DEBUG_LEVEL 2
#endif

#ifndef MCMC_DEBUG
#define MCMC_DEBUG(level) if (level <= (MCMC_DEBUG_LEVEL)) std::cout
#endif

#define MCMC_ERROR (std::cout <<__FILE__<<":: " << __LINE__ << ": " )

class TProposeAdaptiveStep;

/// A ROOT based include-file-only templated class to run an MCMC.  The
/// resulting MCMC normally uses the Metropolis-Hastings algorithm with an
/// adaptive proposal function.  The UserLikelihood template argument must be
/// a class (or struct) which provides a method declared as:
///
///\code
/// struct ExampleLogLikelihood {
///    double operator() (const std::vector<double>& point);
/// }
///\endcode
///
/// The optional UserProposal template argument must be a struct (or class)
/// which provides a method declared as:
///
///\code
/// struct ExampleProposal {
///    void operator() (std::vector<double>& proposed,
///                       const std::vector<double>& previous,
///                       const double previousValue);
/// }
///\endcode
///
/// where proposed is the new proposed point, previous is the last
/// accepted point, and previousValue is the log Likelihood at the
/// previous point.
///
/// This can be used in your root macros:
///
/// \code
/// class TDummyLogLikelihood {
/// public:
///     double operator()(const Vector& point)  const {
///         double logLikelihood = 0.0;
///         for (std::size_t i = 0; i < point.size(); ++i) {
///             logLikelihood += - 0.5*point[i]*point[i];
///         }
///         return logLikelihood;
///     }
/// };
///
/// void SimpleMCMC() {
///
///     TFile *outputFile = new TFile("simple-mcmc.root","recreate");
///     TTree *tree = new TTree("SimpleMCMC",
///                             "Tree of accepted points");
///
///     TSimpleMCMC<TDummyLogLikelihood> mcmc(tree);
///
///     // The next three lines are for example only and should almost never
///     // be used.  This is to show the syntax for controlling
///     // TProposeAdaptiveStep, don't just copy this blindly!
///     mcmc.GetProposeStep().SetDim(5);          // Not needed!
///     mcmc.GetProposeStep().SetGaussian(3,2.0); // Not recommended!
///     mcmc.GetProposeStep().SetUniform(4,-5,5); // If needed for special case.
///
///     Vector point(5);
///     mcmc.Start(point);  // Set initial conditions
///
///     for (int i=0; i<1000000; ++i) mcmc.Step(false);  // Burn-in the chain.
///     for (int i=0; i<1000000; ++i) mcmc.Step();       // Run the chain.
///
///     tree->Write();
///     delete outputFile;
/// }
/// \endcode
///
/// If this macro is in a file named "SimpleMCMC.C", then you can run it using
///
/// \code
/// root -l -q SimpleMCMC.C+
/// \code
///
/// The macro must be compiled using ACLIC.
///
/// The default class for UserProposal is TProposeAdaptiveStep which
/// implements an adaptive Metropolis-Hastings step.  It has several methods
/// that can be accessed using the GetProposeStep() method.  See above for an
/// example.
///
/// \note Copyright 2017 Clark McGrew (details at the end of the source file).
/// The full distribution can be found at
/// https://github.com/ClarkMcGrew/root-simple-mcmc
///
template <typename UserLikelihood,
          typename UserProposal = TProposeAdaptiveStep>
class TSimpleMCMC {
public:

    /// Make the likelihood class available as TSimpleMCMC::LogLikelihood.
    typedef UserLikelihood LogLikelihood;

    /// Make the step proposal cclass available as TSimpleMCMC::ProposeStep.
    typedef UserProposal ProposeStep;

    /// Declare an object to run an MCMC.  The resulting MCMC normally uses
    /// the Metropolis-Hastings algorithm with an adaptive proposal function.
    /// This takes an optional pointer to a tree to save the accepted steps.
    /// If it is provided the constructor will add branch to the tree for the
    /// log likelihood and the parameters at the accepted step.  If the second
    /// optional parameter is true, then the proposed steps will also be added
    /// to the tree.
    TSimpleMCMC(TTree* tree = NULL, bool saveStep = false) : fTree(tree) {
        if (fTree) {
            MCMC_DEBUG(0) << "TSimpleMCMC: Adding branches to "
                          << fTree->GetName()
                          << std::endl;
            fTree->Branch("LogLikelihood",&fAcceptedLogLikelihood);
            fTree->Branch("Accepted",&fAccepted);
            if (saveStep) {
                MCMC_DEBUG(0) << "TSimpleMCMC: Saving the trial steps."
                              << std::endl;
                fTree->Branch("Step",&fTrialStep);
            }
        }
        fLogLikelihoodCount = 0;
    }

    /// Get a reference to the object that will propose the step.  The
    /// ProposeStep object is constructed by the TSimpleMCMC template.  You
    /// can get write access to that object with this method, and can access
    /// any of your declared methods.  Because this is a template, ProposeStep
    /// is actually just a typedef for your class.
    ProposeStep& GetProposeStep() {return fProposeStep;}

    /// Get a reference to the likelihood calculation object.  The
    /// LogLikelihood object is constructed by the TSimpleMCMC template and
    /// is the class that you handed to the template.  You have write access
    /// to that object with this method, and can access any of your declared
    /// methods.  Because this is a template, TSimpleMCMC::LogLikelihood is
    /// actually just a typedef for your class.
    LogLikelihood& GetLogLikelihood() {return fLogLikelihood;}

    /// Get the number of times the log likelihood has been called.
    int GetLogLikelihoodCount() {return fLogLikelihoodCount;}

    /// Set the starting point for the mcmc.  If the optional argument is
    /// true, then the point will be saved to the output.
    void Start(Vector start, bool save=true) {
        fProposed.resize(start.size());
        std::copy(start.begin(), start.end(), fProposed.begin());
        
        fAccepted.resize(start.size());
        std::copy(start.begin(), start.end(), fAccepted.begin());

        fTrialStep.resize(start.size());
        
        fProposedLogLikelihood = GetLogLikelihoodValue(fProposed);
        fAcceptedLogLikelihood = fProposedLogLikelihood;

        if (save) SaveStep();
    }
    
    /// Take a step.  This returns true if a new point has been accepted, and
    /// false if we stay at the old point.  If save is true, then the points
    /// are saved to the output.
    bool Step(bool save=true) {
        if (fProposed.empty() || fAccepted.empty()) {
            MCMC_ERROR << "Must initialize starting point" << std::endl;
            throw;
        }

        fProposeStep(fProposed,fAccepted,fAcceptedLogLikelihood);

        // Only cache the trial step when tree is being saved.
        if (save) {
            for (std::size_t i = 0; i < fProposed.size(); ++i) {
                fTrialStep[i] = fProposed[i] - fAccepted[i];
            }
        }

        // Find the likelihood at the new step.  The old likelihood has been
        // cached.
        fProposedLogLikelihood = GetLogLikelihoodValue(fProposed);
        double delta = fProposedLogLikelihood - fAcceptedLogLikelihood;
        if (delta < 0.0 ) {
            // The proposed likelihood is less than the previously accepted
            // likelihood, so see if it should be rejected.
            double trial = std::log(gRandom->Uniform());
            if (delta < trial) {
                // The new step should be rejected, so save the old step.
                // This depends on IEEE error handling so that std::log(0.0)
                // is -inf which is always less than delta.
                if (save) SaveStep();
                return false;
            }
        }

        // We're keeping a new step.
        std::copy(fProposed.begin(), fProposed.end(), fAccepted.begin());
        fAcceptedLogLikelihood = fProposedLogLikelihood;

        // Save the information to the output tree.
        if (save) SaveStep();
        return true;
    }

    /// Get the likelihood at the most recently accepted point.
    double GetAcceptedLogLikelihood() const {return fAcceptedLogLikelihood;}

    /// Get the most recently accepted point.
    const Vector& GetAccepted() const {return fAccepted;}

    /// Get the likelihood at the most recently proposed point.
    double GetProposedLogLikelihood() const {return fProposedLogLikelihood;}

    /// Get the most recently proposed point.
    const Vector& GetProposed() const {return fProposed;}
    
protected:

    /// If possible, save the step.
    void SaveStep() {if (fTree) fTree->Fill();}

    double GetLogLikelihoodValue(const Vector& point) {
        ++fLogLikelihoodCount;
        return fLogLikelihood(point);
    }
    
    /// A class (called as a functor) to calculate the likelhood.
    LogLikelihood fLogLikelihood;

    /// A class (called as a functor) to propose the next step.
    ProposeStep fProposeStep;

    /// A TTree to save the accepted points.
    TTree* fTree;

    /// The number of times the likelihood has been calculated.
    int fLogLikelihoodCount;
    
    /// The last accepted point.  This will be the same as the proposed point
    /// if the last step was accepted.
    Vector fAccepted;

    /// The likelihood at the last accepted pont.
    double fAcceptedLogLikelihood;
    
    /// The trial step.  The difference between the previous point and the new
    /// proposed point.
    Vector fTrialStep;

    /// The proposed point for the most recent step (This may be the same as
    /// the accepted point).
    Vector fProposed;

    /// The likelihood at the last proposed point.
    double fProposedLogLikelihood;
};

// This is a very simple example of a step proposal class.  It's not actually
// used for anything, but implements a Metropolis-Hastings step.  The proposal
// is the new point to be tried.  The current is the last successful step, and
// the value is the log likelihood of the last successful step.  The step is
// taken relative to the current and the value is ignored.
struct TProposeSimpleStep {
    TProposeSimpleStep(): fSigma(-1.0) {}
    
    void operator ()(Vector& proposal,
                     const Vector& current,
                     const double value) const {
        double sigma = fSigma;
        
        // No width was provided, so make a bogus guess at a reasonable width;
        if (sigma < 0.0) sigma = std::sqrt(1.0/proposal.size());

        // Make the proposal.
        for (std::size_t i = 0; i < proposal.size(); ++i) {
            proposal[i] = current[i] + gRandom->Gaus(0.0,sigma);
        }
    }

    double fSigma;
};

/// A default for the class to propose the next step.  This implements an
/// adaptive Metropolis-Hastings proposal.  It starts with a guess at the
/// covariance of the posterior, and the updates the estimated posterior.  At
/// user defined intervals (normally about the dimensionality squared), the
/// estimated covariance is updated using the current state of the Markov
/// chain. (Notice that this means it's not really a Markov Chain!  You need
/// to check the ergodcity, but it's almost always OK.)  It works OK as long
/// as the posterior is more or less Gaussian.  If the posterior is not
/// Gaussian, then this probably won't fail, but it can become less efficient.
class TProposeAdaptiveStep {
public:
    TProposeAdaptiveStep() :
        fLastValue(0.0), fTrials(0), fSuccesses(0), fAcceptanceWindow(-1),
        fCovarianceWindow(-1), fNextUpdate(-1), fAcceptance(0.0), fSigma(0.0),
        fStateInitialized(false) {
        // Set a default value for the target acceptance rate.  For sum
        // reason, the magic value in the literature is 44%.
        fTargetAcceptance = 0.44;
    }

    /// Take a proposed trial point to fill, the current point, and the
    /// likelihood at the current point.
    void operator ()(Vector& proposal,
                     const Vector& current,
                     const double value) {
        if (proposal.size() != current.size()) {
            // Apply a sanity check.  This MUST be true.
            MCMC_ERROR << "Proposal and current vectors must be same size."
                       << std::endl;
            throw;
        }

        UpdateState(current,value);

        std::copy(current.begin(), current.end(), proposal.begin());
        
        // Make the proposal.
        for (std::size_t i = 0; i < proposal.size(); ++i) {
            if (fProposalType[i].type == 1) {
                // Make a uniform proposal.
                proposal[i] = gRandom->Uniform(fProposalType[i].param1,
                                               fProposalType[i].param2);
                continue;
            }
            // Make a Gaussian Proposal (with the latest estimate of the
            // covariance).
            double r = gRandom->Gaus(0.0,1.0);
            for (std::size_t j = 0; j < proposal.size(); ++j) {
                if (fProposalType[j].type == 1) continue;
                proposal[j] += fSigma*r*fDecomposition(i,j);
            }
        }
    }

    const Vector& GetEstimatedCenter() const {return fCentralPoint;}
    
    /// Set the number of dimensions in the proposal.  This must match the
    /// dimensionality of the likelihood being use.
    void SetDim(int dim) {
        if (!fLastPoint.empty()) {
            // Apply a sanity check.  This MUST be true.
            MCMC_ERROR << "Dimensionality has already been set."
                       << std::endl;
            return;
        }
        fLastPoint.resize(dim);
        fProposalType.resize(dim);
    }

    /// Set the proposal function for a particular dimension to be uniform.
    void SetUniform(int dim, double minimum, double maximum) {
        if (dim < 0 || (std::size_t) dim >= fProposalType.size()) {
            MCMC_ERROR << "Dimension " << dim << " is out of range."
                       << " 0 to " << fProposalType.size()
                       << std::endl;
            return;
        }
        MCMC_DEBUG(0) << "Overriding proposal for dimension " << dim
                      << " to be uniform between "
                      << "[" << minimum
                      << ", " << maximum << "]."
                      << std::endl;
        fProposalType[dim].type = 1;
        fProposalType[dim].param1 = minimum;
        fProposalType[dim].param2 = maximum;
    }

    /// Set the proposal function for a particular dimension to be Gaussian.
    /// This is the default, so it's only needed if you need to give a hint
    /// about the width of the proposal (for instance the dimension is very
    /// wide, or very narrow).
    void SetGaussian(int dim, double sigma) {
        if (dim < 0 || (std::size_t) dim >= fProposalType.size()) {
            MCMC_ERROR << "Dimension " << dim << " is out of range."
                       << std::endl;
            return;
        }
        MCMC_DEBUG(0) << "Overriding proposal for dimension " << dim
                      << " to be Gaussian with "
                      << sigma << " sigma."
                      << std::endl;
        fProposalType[dim].type = 0;
        fProposalType[dim].param1 = sigma*sigma;
    }
    
    /// Set the window over which to estimate the covariance.  This can
    /// normally be very large (and there is a reasonable default), but it
    /// might need to be smaller for some pathelogical distributions.
    void SetCovarianceWindow(int w) {fCovarianceWindow = w;}
    
    /// The proposal steps are chosen based on an estimate of the covariance
    /// of the posterior.  This method forces the covariance of the proposal
    /// to be updated.  It's automatically called during the run, so (while it
    /// can be) it probably doesn't need to be called in user code.
    void UpdateProposal() {
        MCMC_DEBUG(1) << "Update after "
                      << fSuccesses << "/" << fTrials << " successes"
                      << " (Accepting: " << fAcceptance
                      << " w/ width: " << fSigma << ")"
                      << std::endl;

        MCMC_DEBUG(1) << " Covariance estimated with window of"
                      << " " << fCovarianceTrials
                      << std::endl;
        if (MCMC_DEBUG_LEVEL>1 && fLastPoint.size() < 5) {
            fCurrentCov.Print();
        }

        double trace = 0.0;
        for (std::size_t i=0; i<fLastPoint.size(); ++i) {
            trace += fCurrentCov(i,i);
        }
        MCMC_DEBUG(1) << " Covariance Trace: " << trace
                      << std::endl
                      << "        = ";
        for (std::size_t i=0; i<fLastPoint.size(); ++i) {
            MCMC_DEBUG(1) << fCurrentCov(i,i);
            if (i<fLastPoint.size()-1) MCMC_DEBUG(1) << " + ";
            if (i%6 == 5) MCMC_DEBUG(1) << std::endl << "           ";
        }
        MCMC_DEBUG(1) << std::endl;
        
        // The proposal is being updated, so deweight the current covariance
        // so that the new information weighs more.
        fCovarianceTrials = std::max(1000.0,0.1*fCovarianceTrials);
        fCovarianceTrials = std::min(fCovarianceTrials,0.1*fCovarianceWindow);
        
        // The proposal is being updated, so deweight the current acceptance 
        // so that the new information weighs more.
        fAcceptanceTrials = std::max(1000.0,0.1*fAcceptanceTrials);
        fAcceptanceTrials = std::min(fAcceptanceTrials,0.1*fAcceptanceWindow);
        
        TDecompChol chol(fCurrentCov);
        if (chol.Decompose()) {
            fDecomposition = chol.GetU();
            return;
        }

        // The Cholesky decomposition has failed.  That usually means that the
        // current estimate of the covariance has a one or more pairs of
        // variables that are too correlated, or that the variance of one of
        // the variables has become to small.  This can happen when the chain
        // gets "stuck" at a point.  The solution tried here is to look
        // variables that have a very small variance and to increase it.  It
        // also looks for pairs of variables that have a large correlation and
        // "manually" reduces the correlation.

        // Check for very small variances.
        for (std::size_t i=0; i<fLastPoint.size(); ++i) {
            // The prechain estimate of the variance for the variable.
            double expectedVariance = 1.0;
            if (fProposalType[i].type == 0 && fProposalType[i].param1>0) {
                // A Gaussian proposal with a variance other than one...
                expectedVariance = fProposalType[i].param1;
            }
            else if (fProposalType[i].type == 1) {
                // A uniform proposal...
                expectedVariance = fProposalType[i].param2;
                expectedVariance -= fProposalType[i].param1;
                expectedVariance = expectedVariance*expectedVariance/12.0;
            }
            else {
                MCMC_ERROR << "Illegal proposal type for dimension " << i
                           << ": Type is " << fProposalType[i].type
                           << std::endl;
                throw;
            }
            // The minimum scale of the covariance relative to the prechain
            // estimate of the variance.  This should be picked to make sure
            // that there aren't numeric problems...
            double minimum
                = std::sqrt(std::numeric_limits<Parameter>::epsilon());
            if (fCurrentCov(i,i) < minimum*expectedVariance) {
                MCMC_DEBUG(1) << "Variance for dimension " << i
                              << " has been increased from " << fCurrentCov(i,i)
                              << " to " << minimum*expectedVariance
                              << std::endl;
                fCurrentCov(i,i) = minimum*expectedVariance;
            }
        }

        // Check for very large correlations
        for (std::size_t i=0; i<fLastPoint.size(); ++i) {
            for (std::size_t j=i+1; j<fLastPoint.size(); ++j) {
                double correlation = fCurrentCov(i,j);
                correlation /= std::sqrt(fCurrentCov(i,i));
                correlation /= std::sqrt(fCurrentCov(j,j));
                // Don't worry about "small" correlations.
                double maxCorrelation = 0.95;
                if (correlation < maxCorrelation) continue;
                // Oops, the correlation is too large, so reduce it..
                fCurrentCov(i,j) = maxCorrelation*maxCorrelation;
                fCurrentCov(i,j) *= std::sqrt(fCurrentCov(i,i));
                fCurrentCov(i,j) *= std::sqrt(fCurrentCov(j,j));
                fCurrentCov(j,i) = fCurrentCov(i,j);
            }
        }

        // Make another attempt at finding the Cholesky decomposition.
        TDecompChol chol2(fCurrentCov);
        if (chol2.Decompose()) {
            fDecomposition = chol2.GetU();
            return;
        }

        // Something is going very wrong, so reset the Proposal.  This
        // probably won't work, but it pays to try.
        ResetProposal();
    }

    /// Forget information about the covariance, and use the last point as the
    /// new central value.  This can be useful after burnin to completely
    /// forget about the path to stocastic equilibrium since the covariance is
    /// reset back to the initial conditions.
    void ResetProposal() {
        MCMC_DEBUG(2) << "Reset the proposal after "
                      << fSuccesses << " successes "
                      << " in " << fTrials << " trials "
                      << std::endl;
        MCMC_DEBUG(2) << " Recent acceptance rate was " << fAcceptance
                      << " with an adjusted width of " << fSigma
                      << std::endl;
        // Reset the success and trials counts.
        fTrials = 0;
        fSuccesses = 0;
        // Take a wild guess at width to get the right acceptance.
        if (fSigma < 0.01*std::sqrt(1.0/fLastPoint.size())) {
            fSigma = std::sqrt(1.0/fLastPoint.size());
        }
        // Setup the spece for the decomposition.
        fDecomposition.ResizeTo(fLastPoint.size(), fLastPoint.size());
        // Set up the initial estimate of the covariance.
        fCurrentCov.ResizeTo(fLastPoint.size(), fLastPoint.size());
        for (std::size_t i = 0; i < fLastPoint.size(); ++i) {
            for (std::size_t j = i; j < fLastPoint.size(); ++j) {
                if (i == j
                    && fProposalType[i].type == 0
                    && fProposalType[i].param1 > 0) {
                    MCMC_DEBUG(0) << "Overriding covariance for dimension "
                                  << i 
                                  << " from "
                                  << fCurrentCov(i,i)
                                  << " to " 
                                  << fProposalType[i].param1
                                  << std::endl;
                    fCurrentCov(i,i) = fProposalType[i].param1;
                }
                else if (i == j) {
                    fCurrentCov(i,i) = 1.0;
                }
                else fCurrentCov(i,j) = fCurrentCov(j,i) = 0.0;
            }
        }
        // Set a default window to average the covariance over.  This
        // basically averages over everything.
        int minWindow = 1000;
        if (fCovarianceWindow < minWindow) fCovarianceWindow = 10000000;
        // After reseting, erase the acceptance history.
        fAcceptance = fTargetAcceptance;
        fAcceptanceTrials = std::min(10.0, 0.5*fAcceptanceWindow);
        // Initialize the central point.
        fCentralPoint.resize(fLastPoint.size());
        std::copy(fLastPoint.begin(), fLastPoint.end(), fCentralPoint.begin());
        // Start with a non-zero number of trials to represent the prior
        // information coming from the first guess at the central point.
        fCentralPointTrials =  std::min(10.0, 0.1*fCovarianceWindow);
        // This makes sure everything is set properly.
        UpdateProposal();
    }
    
private:

    // Return to a default state.
    void InitializeState(const Vector& current, const double value) {
        if (fStateInitialized) return;
        fStateInitialized = true;
        if (fLastPoint.size() < 1) {
            SetDim(current.size());
        }
        else if (fLastPoint.size() != current.size()) {
            // Sanity check! These must be equal.
            MCMC_ERROR << "Mismatch in the dimensionality."
                       << std::endl;
        }
        fLastValue = value;
        std::copy(current.begin(), current.end(), fLastPoint.begin());
        // Set a default window to average the acceptance over.
        fAcceptanceWindow = std::pow(1.0*fLastPoint.size(),2.0) + 1000;
        // The steps until the next update.
        fNextUpdate = fAcceptanceWindow;
        // Reset the proposal covariance as part of the initialization.
        ResetProposal();
    }
    
    /// This updates the current state.  The new proposals are adjusted based
    /// on the past history of success (or failure).  This helps make the
    /// chain more efficient at exploring the posterior.
    void UpdateState(const Vector& current, const double value) {
        InitializeState(current,value);
        ++fTrials;

        // Make a quick check to see if the point moved.  This gets it right
        // more than 99.9% of the time (i.e. it's Good Enough (tm)).
        bool accepted = false;
        if (value != fLastValue || current[1] != fLastPoint[1]) accepted = true;

        // Track the total number of successes.
        if (accepted) ++fSuccesses;
        
        // Update the acceptance over the last "fAcceptanceWindow" trials.  
        fAcceptance *= fAcceptanceTrials;
        if (accepted) fAcceptance = fAcceptance+1.0;
        fAcceptance /= fAcceptanceTrials + 1.0;
        fAcceptanceTrials = std::min(fAcceptanceWindow,fAcceptanceTrials+1.0);

        // Update the scale for the proposal step.  If the acceptance is to
        // small, then the scale should be reduced.  If the acceptance is to
        // large, then the scale should be increase.  This is set so that if
        // the ratio of measured acceptance to the target acceptance is one,
        // then the value of fSigma will not change.  The exponent is chosen
        // to control how quickly the fSigma value will be changed.w
        fSigma *= std::pow(fAcceptance/fTargetAcceptance,
                           std::min(0.001,0.5/fAcceptanceWindow));
        
        // Update the estimate of the central value.  This is simply a running
        // average of the position of the points in the posterior.
        for (std::size_t i=0; i< current.size(); ++i) {
            fCentralPoint[i] *= fCentralPointTrials;
            fCentralPoint[i] += current[i];
            fCentralPoint[i] /= fCentralPointTrials + 1;
        }
        fCentralPointTrials = std::min(fCovarianceWindow,
                                       fCentralPointTrials+1.0);

        // Update the estimate of the covariance.  This is a running
        // calculation of the covariance.
        for (std::size_t i=0; i<current.size(); ++i) {
            for (std::size_t j=0; j<i+1; ++j) {
                double v = fCurrentCov(i,j);
                double r = (current[i]-fCentralPoint[i])
                    *(current[j]-fCentralPoint[j]);
                v *= fCovarianceTrials;
                v += r;
                v /= fCovarianceTrials + 1.0;
                if (i == j) fCurrentCov(i,j) = v;
                else fCurrentCov(i,j) = fCurrentCov(j,i) = v;
            }
        }
        fCovarianceTrials = std::min(fCovarianceWindow,
                                     fCovarianceTrials+1.0);
        
        // Periodically change how the proposal is updated based on the
        // current estimate of the covariance.
        if (accepted && (--fNextUpdate)<1) {
            fNextUpdate= fAcceptanceWindow + 0.5*fSuccesses;
            UpdateProposal();
        }
        
        // Save the last value and point.
        fLastValue = value;
        std::copy(current.begin(), current.end(), fLastPoint.begin());
    }
    
    // The previous current point.  This is used to (among other things) keep
    // track of when the state has changed.
    Vector fLastPoint;

    // The previous log Likelihood.
    double fLastValue;

    // The current (running) estimate for the central value of each parameter.
    Vector fCentralPoint;

    // The trials being used for the current estimated central value.  This
    // will be a value between one and fCovarianceWindow.
    double fCentralPointTrials;
    
    // The current (running) estimate of the covariance.
    TMatrixD fCurrentCov;

    // The trials being used for the current estimated covariance.  This
    // will be a value between one and fCovarianceWindow.
    double fCovarianceTrials;
    
    // The window to average the covariance and central value over.  The
    // covariance window should be several times the autocorrelation period.
    // It should usually be very large so that all to points used to probe the
    // posterior are used.  Setting it to a smaller value means that the local
    // covariance of the posterior will be used.  This might be important if
    // the posterior is extremely non-Gaussian (e.g. it's a "banana
    // posterior").
    double fCovarianceWindow;

    // The current Cholesky decomposition of the covariance.  This is not
    // updated everytime the fCurrentCov estimate changes.
    TMatrixD fDecomposition;
    
    // Record the type of proposal to use for each dimension
    struct ProposalType {
        ProposalType(): type(0), param1(0), param2(0) {}
        int type; // 0 for Gaussian, 1 for Uniform.
        double param1; // Sigma for Gaussian, Minimum for Uniform
        double param2; // Not used for Gaussian, Maximum for Uniform
    };

    // The type of distribution to draw the propsal for a dimension from.
    std::vector<ProposalType> fProposalType;
    
    // The number of times a step has been proposed since the last reset.
    int fTrials;

    // The total number of successes since the last reset
    int fSuccesses;
    
    // A counter until the next update gets done.
    int fNextUpdate;
    
    // The recent acceptance rate.
    double fAcceptance;

    // The trials being used for the current estimated acceptance rate.  This
    // will be a value between one and fAcceptanceWindow.
    double fAcceptanceTrials;
    
    // The window to average the acceptance over.
    double fAcceptanceWindow;

    // The target acceptance rate for the chain.
    double fTargetAcceptance;
    
    // The current width of the proposal.
    double fSigma;

    // Keep track of whether we've actually been called.
    bool fStateInitialized;
};

// MIT License

// Copyright (c) 2017 Clark McGrew

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#endif
