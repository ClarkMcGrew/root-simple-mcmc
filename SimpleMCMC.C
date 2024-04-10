#include "TSimpleMCMC.H"

#undef LIKELIHOOD_DEFINED

// #define USE_HORRIFIC_LIKELIHOOD
#ifndef LIKELIHOOD_DEFINED
#ifdef USE_HORRIFIC_LIKELIHOOD
#define LIKELIHOOD_DEFINED
#warning Using HORRIFIC likelihood
#include "THorrificLogLikelihood.H"
typedef THorrificLogLikelihood TDummyLogLikelihood;
#endif
#endif

// #define USE_HARD_LIKELIHOOD
#ifndef LIKELIHOOD_DEFINED
#ifdef USE_HARD_LIKELIHOOD
#define LIKELIHOOD_DEFINED
#warning Using HARD likelihood
#include "THardLogLikelihood.H"
typedef THardLogLikelihood TDummyLogLikelihood;
#endif
#endif

// #define USE_ASYM_LIKELIHOOD
#ifndef LIKELIHOOD_DEFINED
#ifdef USE_ASYM_LIKELIHOOD
#define LIKELIHOOD_DEFINED
#warning Using ASYM likelihood
#include "TAsymLogLikelihood.H"
typedef TASymLogLikelihood TDummyLogLikelihood;
#endif
#endif

#ifndef LIKELIHOOD_DEFINED
#define LIKELIHOOD_DEFINED
#warning Using normal likelihood
#include "TDummyLogLikelihood.H"
#endif

#include <TRandom3.h>

#include <sstream>

void SimpleMCMC(int cycles, int steps,
                const char* outputName,
                const char* restoreName) {
    std::cout << "Simple MCMC Loaded" << std::endl;

    TFile *restoreFile = NULL;
    TTree *restoreTree = NULL;
    if (restoreName) {
        std::cout << "Restore from " << restoreName << std::endl;
        restoreFile = new TFile(restoreName,"old");
        restoreTree = (TTree*) restoreFile->Get("SimpleMCMC");
    }

#ifdef NO_OUTPUT
    TFile *outputFile = NULL;
    TTree *tree = NULL;
#else
    TFile *outputFile = new TFile(outputName,"recreate");
    TTree *tree = new TTree("SimpleMCMC","Tree of accepted points");
    tree->SetDirectory(outputFile);
#endif

    // Initialize the random number generator.  This makes sure that a
    // different sequence is used everytime the code is run.
    gRandom = new TRandom3(0);

#ifdef USE_THIS_PROPOSAL
#undef USE_THIS_PROPOSAL
    // This is the same as the default (below).  It's here for completeness
    // and testing.
    TSimpleMCMC<TDummyLogLikelihood,TProposeAdaptiveStep> mcmc(tree);
#endif

#ifdef USE_THIS_PROPOSAL
#undef USE_THIS_PROPOSAL
    // This is just the simple step.  This is the best choice with SKIP_MCMC
    // is defined (see #define below).  I wouldn't normally suggest it
    // otherwise, but since it's a lot faster than the adaptive step [O(D) vs
    // O(D^3)] it may be worth a try if your likelihood function is extremely
    // fast and you a bazillion steps (I've never seen that in practice).
    TSimpleMCMC<TDummyLogLikelihood,TProposeSimpleStep> mcmc(tree);
    // Should be adjusted for each likelihood.
    mcmc.GetProposeStep().fSigma = 0.001;
#endif

#define USE_THIS_PROPOSAL
#ifdef USE_THIS_PROPOSAL
#undef USE_THIS_PROPOSAL
    // Create the MCMC object and get the likelihood.  Probably best for
    // everything except debugging the likelihood with SKIP_MCMC.
    sMCMC::TSimpleMCMC<TDummyLogLikelihood> mcmc(tree,true);
#endif

    TDummyLogLikelihood& like = mcmc.GetLogLikelihood();

    // Initialize the likelihood (if you need to).  The dummy likelihood
    // setups a covariance to make the PDF more interesting.
    like.Init();

    // Set the number of dimensions for the proposal.
    mcmc.GetProposeStep().SetDim(like.GetDim());

#ifdef IMPOSE_RANDOM_CORRELATIONS
    for (int i=0; i<like.GetDim(); ++i) {
        for (int j=i+1; j<like.GetDim(); ++j) {
            double c = gRandom->Uniform(-0.999,0.999);
            if (i == 2 && j == 3) c = 1.0/0.0;
            mcmc.GetProposeStep().SetCorrelation(i,j,c);
        }
    }
#endif

/// Uncomment this to "minimize" the likelihood
// #define SKIP_MCMC

/// Uncomment this to "scan" the likelihood.  Set the value to the variable to
/// be scanned.
// #define SCAN_MCMC 1

/// Uncomment this to do a separate burnin stage.  Usually you shouldn't use
/// this and then skip the first "N" entries of the chain.  "N" is a determined
/// based on the autocorrelation.
#define BURNIN_CHAIN

    // Set one of the dimensions to use a uniform proposal over a fixed range.
    // mcmc.GetProposeStep().SetUniform(1,-0.5,0.5);

    // Override the target acceptance (this value is for very low dimension).
    // mcmc.GetProposeStep().SetTargetAcceptance(0.44);

#ifdef SCAN_MCMC
    mcmc.GetProposeStep().SetScanDimension(SCAN_MCMC);
#endif

    // The number of dimensions in the point needs to agree with the number of
    // dimensions in the likelihood.  You can either hard code it, or do like
    // I'm doing here and have a likelihood method to return the number of
    // dimensions.
    sMCMC::Vector p(like.GetDim());
    // Randomize the starting point
    for (std::size_t i=0; i<p.size(); ++i) p[i] = gRandom->Uniform(-0.99,0.99);
    // Override p and start near the global minimum for the Rosenbrock function
    for (std::size_t i=0; i<p.size(); ++i) p[i] = gRandom->Uniform(0.50,1.5);
    // Override p and start at zero.
    for (std::size_t i=0; i<p.size(); ++i) p[i] = 0.0;

    mcmc.Start(p,false);

    if (restoreTree) {
        mcmc.Restore(restoreTree);
        delete restoreFile;
        std::cout << "State Restored" << std::endl;
    }

    int verbosity = steps/5;
    int freezeAfter = 0.1*cycles + 1;
    freezeAfter = -1;  // Uncomment to never freeze.

#if defined(BURNIN_CHAIN) && !defined(SKIP_MCMC)
    // This can be useful for debugging, but isn't good practice.  You should
    // be looking at the burn-in steps to make sure the burn-in is complete.
    std::cout << "Start burn-in chain" << std::endl;

    int aWin = 0.1*steps;
    if (aWin > 1000) aWin = 1000;
    if (aWin < 100) aWin = 100;
    mcmc.GetProposeStep().SetAcceptanceWindow(aWin);
    mcmc.GetProposeStep().SetCovarianceWindow(steps);

    int burnin = steps;
    // Burnin the chain (don't save the output)
    for (int i=0; i<burnin; ++i) {
        if (i%verbosity == 0) {
            std::cout << "First burn-in " << i << " / " << burnin << std::endl;
        }
        mcmc.Step(false);
    }
    mcmc.GetProposeStep().ResetProposal();
    std::cout << "Finished first burn-in chain" << std::endl;

    mcmc.GetProposeStep().SetAcceptanceWindow(aWin);
    mcmc.GetProposeStep().SetCovarianceWindow(2*steps);
    mcmc.GetProposeStep().SetCovarianceUpdateDeweighting(0.5);

    for (int cycle=0; cycle < 4; ++cycle) {
        // Burnin the chain (don't save the output)
        for (int i=0; i<burnin; ++i) {
            if (i%verbosity == 0) {
                std::cout << "burn-in " << cycle
                          << " step " << i << " / " << burnin << std::endl;
            }
            mcmc.Step();
        }
        mcmc.GetProposeStep().UpdateProposal();
        std::cout << "Finished burnin chain " << cycle << std::endl;
    }
#endif

    // Freeze the step size.
    mcmc.GetProposeStep().SetAcceptanceWindow(1000);
    mcmc.GetProposeStep().SetAcceptanceRigidity(2.0);
    mcmc.GetProposeStep().SetCovarianceWindow(cycles*steps);
    mcmc.GetProposeStep().SetCovarianceUpdateDeweighting(0.20);
    mcmc.GetProposeStep().SetNextUpdate(1E+9);

    // Run the chain (with output to the tree).
    int trial = 0;
    for (int cycle = 0; cycle < cycles; ++cycle) {
        for (int i=0; i<steps; ++i) {
            ++trial;
            if (trial%verbosity == 0) {
                std::cout << "Trial " << cycle << ": "
                          << trial << " / " << cycles*steps
                          << " A: " << mcmc.GetProposeStep().GetAcceptance()
                          << "/" << mcmc.GetProposeStep().GetSuccesses()
                          << " ("<<mcmc.GetProposeStep().GetAcceptanceWindow()
                          << ":" << mcmc.GetProposeStep().GetNextUpdate() << ")"
                          << " S: " << mcmc.GetProposeStep().GetSigma()
                          << " T: "
                          << mcmc.GetProposeStep().GetCovarianceTrace()
                          << " RMS: " << mcmc.GetStepRMS()
                          << std::endl;
            }

#if !defined(SKIP_MCMC) && !defined(SCAN_MCMC)
            // Make a normal MCMC step.
            mcmc.Step();
#endif
#if defined(SKIP_MCMC)
#warning SKIPPING THE MCMC PART OF THE CHAIN
            // Make a debugging (minimization) step
            mcmc.Step(true,1);
#endif
#if defined(SCAN_MCMC)
#warning SCANNING THE LIKELIHOOD
            // Make a debugging (minimization) step
            mcmc.Step(true,2);
#endif
        }

        std::cout << "UPDATE PROPOSAL" << std::endl;
        mcmc.GetProposeStep().UpdateProposal();
        if (freezeAfter >= 0 && cycle > freezeAfter) {
            mcmc.GetProposeStep().SetAcceptanceRigidity(-1);
            std::cout << "FREEZE STEP SIZE" << std::endl;
        }
        else {
            mcmc.GetProposeStep().SetAcceptanceRigidity(2.0);
        }
        mcmc.GetProposeStep().SetCovarianceUpdateDeweighting(0.0);
        mcmc.GetProposeStep().SetNextUpdate(10*steps);
    }
    std::cout << "Finished with " << mcmc.GetLogLikelihoodCount() << " calls"
              << std::endl;

    // Save the final state.  This is needed to force the proposal to save
    // it's final state.
    mcmc.SaveStep();

    if (tree) {
        std::cout << "Write the tree" << std::endl;
        tree->Write();
    }
    if (outputFile) {
        std::cout << "Close the file" << std::endl;
        delete outputFile;
    }
}

#ifdef MAIN_PROGRAM
// This let's the example compile directly.  To compile it, use the compile.sh
// script and then run it using ./a.out which will produce a file name
// "SimpleAHMC.root"
int main(int argc, char **argv) {
    int cycles = 10;
    int steps = 1000;
    std::string outputName("SimpleMCMC.root");
    char* restoreName = NULL;

    if (argc > 1) {
        std::istringstream input(argv[1]);
        input >> cycles;
    }
    if (argc > 2) {
        std::istringstream input(argv[2]);
        input >> steps;
    }
    if (argc > 3) {
        outputName = argv[3];
    }
    if (argc > 4) {
        restoreName = argv[4];
    }

    SimpleMCMC(cycles,steps,outputName.c_str(),restoreName);
}
#endif
