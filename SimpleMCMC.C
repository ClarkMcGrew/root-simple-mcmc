#include "TSimpleMCMC.H"

#define USE_HARD_LIKELIHOOD
#ifdef USE_HARD_LIKELIHOOD
#include "THardLogLikelihood.H"
typedef THardLogLikelihood TDummyLogLikelihood;
#else
#include "TDummyLogLikelihood.H"
#endif

#include <TRandom3.h>

#include <sstream>

void SimpleMCMC(int trials,
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

    // Create the MCMC object and get the likelihood.
    TSimpleMCMC<TDummyLogLikelihood> mcmc(tree);
    TDummyLogLikelihood& like = mcmc.GetLogLikelihood();

    // Initialize the likelihood (if you need to).  The dummy likelihood
    // setups a covariance to make the PDF more interesting.
    like.Init();

    // Set the number of dimensions for the proposal.
    mcmc.GetProposeStep().SetDim(like.GetDim());

    // Set one of the dimensions to use a uniform proposal over a fixed range.
    // mcmc.GetProposeStep().SetUniform(1,-0.5,0.5);

    // The number of dimensions in the point needs to agree with the number of
    // dimensions in the likelihood.  You can either hard code it, or do like
    // I'm doing here and have a likelihood method to return the number of
    // dimensions.
    Vector p(like.GetDim());
    for (std::size_t i=0; i<p.size(); ++i) p[i] = gRandom->Uniform(-1.0,1.0);
    for (std::size_t i=0; i<p.size(); ++i) p[i] = 1.0;

    mcmc.Start(p,false);

    if (restoreTree) {
        mcmc.Restore(restoreTree);
        delete restoreFile;
        std::cout << "State Restored" << std::endl;
    }

    int verbosity = 100;

#ifndef SKIP_BURNIN
    // This can be useful for debugging, but isn't good practice.  You should
    // be looking at the burn-in steps to make sure the burn-in is complete.
    std::cout << "Start burn-in chain" << std::endl;

    int burnin = 10000+10*p.size();
    // Burnin the chain (don't save the output)
    for (int i=0; i<burnin; ++i) {
        if (i%verbosity == 0) {
            std::cout << "First burn-in " << i << " / " << burnin << std::endl;
        }
        mcmc.Step(false);
    }
    mcmc.GetProposeStep().UpdateProposal();
    std::cout << "Finished first burn-in chain" << std::endl;

    // Burnin the chain (don't save the output)
    burnin = 10000+10.0*p.size()*p.size();
    for (int i=0; i<burnin; ++i) {
        if (i%verbosity == 0) {
            std::cout << "Second burn-in " << i << " / " << burnin << std::endl;
        }
        mcmc.Step(false);
    }
    mcmc.GetProposeStep().UpdateProposal();
    std::cout << "Finished second burn-in chain" << std::endl;
#endif

    // Run the chain (with output to the tree).
    for (int i=0; i<trials; ++i) {
        if (i%verbosity == 0) {
            std::cout << "Trial " << i << " / " << trials << std::endl;
        }
        mcmc.Step();
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
    int trials = 10000;
    std::string outputName("SimpleMCMC.root");
    char* restoreName = NULL;

    if (argc > 1) {
        std::istringstream input(argv[1]);
        input >> trials;
    }
    if (argc > 2) {
        outputName = argv[2];
    }
    if (argc > 3) {
        restoreName = argv[3];
    }

    SimpleMCMC(trials,outputName.c_str(),restoreName);
}
#endif
