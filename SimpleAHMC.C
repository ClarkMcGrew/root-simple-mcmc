#include "TSimpleHMC.H"
#include "TDummyLogLikelihood.H"

#include <TMatrixD.h>
#include <TVectorD.h>

#include <sstream>

void SimpleAHMC(int trials, int maxEvals=-1) {
    std::cout << "Simple AHMC Loaded" << std::endl;

#ifdef NO_OUTPUT
    TFile *outputFile = NULL;
    TTree *tree = NULL;
#else
    TFile *outputFile = new TFile("SimpleAHMC.root","recreate");
    TTree *tree = new TTree("SimpleAHMC","Tree of accepted pqoints");
#endif

#ifdef FORCE_TRUE_GRADIENT
    // Using the full derivative kind of defeats the purpose of an ahmc.
    sMCMC::TSimpleHMC<TDummyLogLikelihood,TDummyLogLikelihood> hmc(tree);
#else
    sMCMC::TSimpleHMC<TDummyLogLikelihood> hmc(tree);
#endif

    TDummyLogLikelihood& like = hmc.GetLogLikelihood();

    // Initialize the likelihood (if you need to).  The dummy likelihood
    // setups a covariance to make the PDF more interesting.
    like.Init();

    // The number of dimensions in the point needs to agree with the number of
    // dimensions in the likelihood.  You can either hard code it, or do like
    // I'm doing here and have a likelihood method to return the number of
    // dimensions.
    sMCMC::Vector p(like.GetDim());
    for (std::size_t i=0; i<p.size(); ++i) p[i] = gRandom->Uniform(-1.0,1.0);

    hmc.Start(p,true);

    int verbose = 500;
    // Burn-in the chain.  This is basically using a Langevin step.
    std::cout << "Start burn-in" << std::endl;
    int burnin = 500 + 2*p.size()*p.size();
    hmc.SetAlpha(0.8);
    hmc.SetMeanEpsilon(-0.1);
    hmc.SetLeapFrog(0);
    for (int i=0; i<burnin; ++i) {
        if (i%verbose == 0) {
            std::cout << "Burn-in: " << i
                      << " Calls: " << hmc.GetPotentialCount()
                      << " Gradients: " << hmc.GetGradientCount()
                      << " Acceptance: " << hmc.GetAcceptanceRate()
                      << std::endl;
        }
        hmc.Step(false,5);
    }

    // Burning the chain a bit more.
    std::cout << "Second burn-in" << std::endl;
    hmc.SetAlpha(0.0);
    hmc.SetMeanEpsilon(-0.05);
    hmc.SetLeapFrog(5);
    burnin = 500 + 2*p.size()*p.size();
    for (int i=0; i<burnin; ++i) {
        if (i%verbose == 0) {
            std::cout << "Burn-in: " << i
                      << " Calls: " << hmc.GetPotentialCount()
                      << " Gradients: " << hmc.GetGradientCount()
                      << " Acceptance: " << hmc.GetAcceptanceRate()
                      << std::endl;
        }
        hmc.Step(false,2);
    }

    // Run the chain
    std::cout << "Run chain" << std::endl;
    hmc.SetAlpha(0.75);
    hmc.SetMeanEpsilon(-0.05);
    hmc.SetLeapFrog(5);
    for (int i=0; i<trials; ++i) {
        if (i%verbose == 0) {
            std::cout << "Trials: " << i
                      << " Calls: " << hmc.GetPotentialCount()
                      << " Gradients: " << hmc.GetGradientCount()
                      << " Acceptance: " << hmc.GetAcceptanceRate()
                      << std::endl;
        }
        // Use the covariant approximation for the derivative.
        hmc.Step(true,2);
    }
    // like.Covariance.Print();
    // hmc.GetEstimatedCovariance().Print();

    std::cout << "Finished " << trials
              << " trials with " << hmc.GetPotentialCount()
              << " calls and " << hmc.GetGradientCount()
              << " gradients "
              << " Accepted: " << hmc.GetAcceptanceRate()
              << std::endl;

    if (tree) tree->Write();
    if (outputFile) delete outputFile;
}

#ifdef MAIN_PROGRAM
// This let's the example compile directly.  To compile it, use the compile.sh
// script and then run it using ./a.out which will produce a file name
// "SimpleAHMC.root"
int main(int argc, char **argv) {
    int trials = 10000;
    int maxEvaluations = -1;
    if (argc > 1) {
        std::istringstream input(argv[1]);
        input >> trials;
    }
    if (argc > 2) {
        std::istringstream input(argv[1]);
        input >> maxEvaluations;
    }
    SimpleAHMC(trials,maxEvaluations);
}
#endif
