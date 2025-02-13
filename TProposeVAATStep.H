#ifndef ProposeVAATStep_H_SEEN
#define ProposeVAATStep_H_SEEN

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

#include <TRandom.h>

#ifndef MCMC_DEBUG_LEVEL
#define MCMC_DEBUG_LEVEL 2
#endif

namespace sMCMC {
    class TProposeVAATStep;
};

/// A default for the class to propose the next step.  This implements an
/// adaptive Variable At A Time.
class sMCMC::TProposeVAATStep {
public:
    TProposeVAATStep() :
        fLastValue(0.0), fTrials(0), fSuccesses(0), fAcceptanceWindow(-1),
        fLastIndex(-1), fAcceptanceRigidity(2.0),
        fStateInitialized(false) {
        // Set a default value for the target acceptance rate.  For some
        // reason, the magic value in the literature is 44%.
        fTargetAcceptance = 0.44;
    }

    bool RestoreState(TTree* tree) {return false;}
    bool AttachState(TTree *tree) {return false;}
    bool SaveState(bool fullSave=false) {return false;}
    bool StateSaved() {return false;}

    /// Take a proposed trial point to fill, the current point, and the
    /// likelihood at the current point.
    void operator ()(sMCMC::Vector& proposal,
                     const sMCMC::Vector& current,
                     const double value) {
        if (proposal.size() != current.size()) {
            // Apply a sanity check.  This MUST be true.
            MCMC_ERROR << "Proposal and current vectors must be same size."
                       << std::endl;
            throw;
        }

        UpdateState(current,value);

        std::copy(current.begin(), current.end(), proposal.begin());

        // Make sure we have a valid proposal
        UpdateProposal();

        // Make the proposal.
        fLastIndex = fNextIndex.back();
        fNextIndex.pop_back();
        if (fProposalType[fLastIndex].type == 1) {
            // Make a uniform proposal.
            proposal[fLastIndex]
                = gRandom->Uniform(fProposalType[fLastIndex].param1,
                                   fProposalType[fLastIndex].param2);
            return;
        }

        double expectedVariance = 1.0;
        if (fProposalType[fLastIndex].type == 0
            and fProposalType[fLastIndex].param1>0) {
            // A Gaussian proposal with a variance other than one...
            expectedVariance = fProposalType[fLastIndex].param1;
        }

        // Make a Gaussian Proposal (with the latest estimate of the
        // covariance).
        proposal[fLastIndex] = current[fLastIndex]
            + fSigma[fLastIndex]*gRandom->Gaus(0.0,expectedVariance);

    }

    /// Set the number of dimensions in the proposal.  This must match the
    /// dimensionality of the likelihood being use.
    void SetDim(int dim) {
        if (!fLastPoint.empty()) {
            // Apply a sanity check.  This MUST be true.
            MCMC_ERROR << "Dimensionality has already been set."
                       << std::endl;
            return;
        }
        fProposalType.resize(dim);
        fLastPoint.resize(dim);
        fAcceptance.resize(dim);
        fAcceptanceTrials.resize(dim);
        fNextIndex.reserve(dim);
        fSigma.resize(dim,2.34);
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
        fProposalType[dim].param1 = sigma;
    }

    /// Set (get) the number of trials used to calculate the
    /// acceptance.
    void SetAcceptanceWindow(double a) {fAcceptanceWindow = a;}
    double GetAcceptanceWindow() const {return fAcceptanceWindow;}

    /// Set (get) the acceptance rigidity.  This is the rough exponential
    /// constant needed to relax the acceptance to the target acceptance.  It
    /// is provided as a factor of the acceptance window.  When the rigidity
    /// is more than 100, the step size is fixed and will not be adapted based
    /// on the acceptance rate.  The value can be as small as 0.5, but typical
    /// values should be around 2.0 or greater.  If the acceptance rigidity is
    /// set to a negative value, the step size will be fixed.
    void SetAcceptanceRigidity(double r) {fAcceptanceRigidity = r;}
    double GetAcceptanceRigidity() const {return fAcceptanceRigidity;}

    /// Get the number of successful steps.
    int GetSuccesses() {return fSuccesses;}

    /// Get the total number of steps taken.
    int GetTrials() {return fTrials;}

    /// Get the current average acceptance rate.
    double GetAcceptance() {
        if (fAcceptance.size() < 1) return 0.0;
        double averageAcceptance
            = std::accumulate(fAcceptance.begin(), fAcceptance.end(),
                              (double) 0.0);
        averageAcceptance /= fAcceptance.size();
        return averageAcceptance;
    }

    /// Get the current average proposal length in standard deviations.
    double GetSigma() {
        if (fSigma.size() < 1) return 0.0;
        double averageSigma
            = std::accumulate(fSigma.begin(), fSigma.end(),
                              (double) 0.0);
        averageSigma /= fSigma.size();
        return averageSigma;
    }

    /// The proposal steps are chosen based on the next dimension to vary.
    void UpdateProposal() {
        if (!fNextIndex.empty()) return;

        // Set initial values (in order).  The indices are popped off as they
        // are used.
        fNextIndex.resize(fProposalType.size());
        fLastIndex = -1;
        for (std::size_t i=0; i< fProposalType.size(); ++i) {
            fNextIndex[i] = i;
        }

        // Do a quick shuffle of the order to break spurious parameter
        // correlations.
        for (std::size_t i=0; i< fProposalType.size(); ++i) {
            std::size_t s = fNextIndex.size() * gRandom->Uniform();
            std::swap(fNextIndex[i],fNextIndex[s]);
        }
    }

    // Return to a default state.
    void InitializeState(const sMCMC::Vector& current, const double value) {
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
        fAcceptanceWindow = 100;
    }

private:

    /// This updates the current state.  The new proposals are adjusted based
    /// on the past history of success (or failure).  This helps make the
    /// chain more efficient at exploring the posterior.
    void UpdateState(const sMCMC::Vector& current, const double value) {
        InitializeState(current,value);
        ++fTrials;

        // Make a quick check to see if the point moved.  This gets it right
        // most of the time (which is good enough).
        bool accepted = false;
        if (value != fLastValue) accepted = true;

        // Track the total number of successes.
        if (accepted) ++fSuccesses;

        // Save the last value and point.
        fLastValue = value;
        std::copy(current.begin(), current.end(), fLastPoint.begin());

        if (fLastIndex < 0) return;

        // Update the acceptance over the last "fAcceptanceWindow" trials.
        ++fAcceptanceTrials[fLastIndex];
        fAcceptance[fLastIndex] *= 1.0*std::min(fAcceptanceWindow,
                                                fAcceptanceTrials[fLastIndex]);
        if (accepted) fAcceptance[fLastIndex] += 1.0;
        fAcceptance[fLastIndex] /=
            1.0 + 1.0*std::min(fAcceptanceWindow,fAcceptanceTrials[fLastIndex]);

        if (fAcceptanceTrials[fLastIndex] > 0.1*fAcceptanceWindow
            and fAcceptanceRigidity > 0
            and fAcceptanceRigidity < 100.0) {
            double v = fSigma[fLastIndex];
            v *= std::pow(fAcceptance[fLastIndex]/fTargetAcceptance,
                          std::min(
                              1.0/500.0,
                              1.0/(fAcceptanceRigidity*fAcceptanceWindow)));
            fSigma[fLastIndex] = std::max(v,1.0E-4);
        }
    }

    // The previous current point.  This is used to (among other things) keep
    // track of when the state has changed.
    sMCMC::Vector fLastPoint;

    // The previous log Likelihood.
    double fLastValue;

    // Record the type of proposal to use for each dimension
    struct ProposalType {
        ProposalType(): type(0), param1(0), param2(0) {}
        int type; // 0 for Gaussian, 1 for Uniform.
        double param1; // Sigma for Gaussian, Minimum for Uniform
        double param2; // Not used for Gaussian, Maximum for Uniform
    };

    // The type of distribution to draw the propsal for a dimension from.
    std::vector<ProposalType> fProposalType;

    // The next dimension to vary.
    std::vector<int> fNextIndex;

    // The last dimension that was varied.
    int fLastIndex;

    // The number of times a step has been proposed.
    int fTrials;

    // The total number of successes
    int fSuccesses;

    // The window to average the acceptance over.
    int fAcceptanceWindow;

    // The rigidity of the acceptance relaxation
    double fAcceptanceRigidity;

    // The recent acceptance rate for each index
    std::vector<double> fAcceptance;

    // The number of trials for each index
    std::vector<int> fAcceptanceTrials;

    // The target acceptance rate for the chain.
    double fTargetAcceptance;

    // The fraction of the parameter window to step.
    std::vector<double> fSigma;

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
