#include "TSimpleAHMC.H"
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

//#define FORCE_TRUE_GRADIENT
#ifdef FORCE_TRUE_GRADIENT
    TSimpleAHMC<TDummyLogLikelihood,TDummyLogLikelihood> hmc(tree);
#else
    TSimpleAHMC<TDummyLogLikelihood> hmc(tree);
#endif
    
    TDummyLogLikelihood& like = hmc.GetLogLikelihood();
    
    // Initialize the likelihood (if you need to).  The dummy likelihood
    // setups a covariance to make the PDF more interesting.
    like.Init();
    
    // The number of dimensions in the point needs to agree with the number of
    // dimensions in the likelihood.  You can either hard code it, or do like
    // I'm doing here and have a likelihood method to return the number of
    // dimensions.
    Vector p(like.GetDim());
    for (std::size_t i=0; i<p.size(); ++i) p[i] = gRandom->Uniform(-1.0,1.0);
    for (std::size_t i=0; i<p.size(); ++i) p[i] = 0.0;
    
    hmc.Start(p,true);
    
#define BURNIN
#ifdef BURNIN
    // Burn-in the chain
    int burnin = p.size();
    hmc.SetAlpha(0.8);
    hmc.SetMeanEpsilon(-0.1);
    hmc.SetLeapFrog(0);
    for (int i=0; i<burnin; ++i) {
        for (int j=0; j<p.size(); ++j) {
            p[j] = hmc.GetCentralPoint()[j];
        }
        // hmc.SetPosition(p);
        for (int j = 0; j<2*burnin; ++j) {
            hmc.Step(false,5);
        }
    }
    // Burning the chain a bit more.
    hmc.SetAlpha(0.0);
    hmc.SetMeanEpsilon(0.05);
    hmc.SetLeapFrog(-5);
    for (int i=0; i<4*p.size()*p.size(); ++i) {
        hmc.Step(false);
    }
#endif
    
    // Run the chain
    hmc.SetAlpha(0.0);
    hmc.SetMeanEpsilon(0.05);
    hmc.SetLeapFrog(-5);
    for (int i=0; i<trials; ++i) {
        if (i%1000 == 0) {
            std::cout << i << " " << hmc.GetPotentialCount()
                      << std::endl;
        }
        hmc.Step(true);
        if (maxEvals > 0 && hmc.GetPotentialCount() > maxEvals) break;
    }
    hmc.GetEstimatedCovariance().Print();
    
    std::cout << "Finished " << trials
              << " trials with " << hmc.GetPotentialCount()
              << " calls and " << hmc.GetGradientCount()
              << " gradients "
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
