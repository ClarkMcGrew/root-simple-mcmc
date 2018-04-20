#include "TSimpleHMC.H"
#include "TDummyLogLikelihood.H"

#include <TMatrixD.h>
#include <TVectorD.h>

#include <sstream>

void SimpleHMC(int trials, int maxEvals=-1) {
    std::cout << "Simple HMC Loaded" << std::endl;

#ifdef NO_OUTPUT
    TFile *outputFile = NULL;
    TTree *tree = NULL;
#else
    TFile *outputFile = new TFile("SimpleHMC.root","recreate");
    TTree *tree = new TTree("SimpleHMC","Tree of accepted pqoints");
#endif

#define FORCE_TRUE_GRADIENT
#ifdef FORCE_TRUE_GRADIENT
    TSimpleHMC<TDummyLogLikelihood,TDummyLogLikelihood> hmc(tree);
#else
    TSimpleHMC<TDummyLogLikelihood> hmc(tree);
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
    
    // Run the chain
    for (int i=0; i<trials; ++i) {
        if (i%1000 == 0) {
            std::cout << i << " " << hmc.GetPotentialCount()
                      << std::endl;
        }
        hmc.Step(true);
        if (maxEvals > 0 && hmc.GetPotentialCount() > maxEvals) break;
    }
    
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
// "SimpleHMC.root"
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
    SimpleHMC(trials,maxEvaluations);
}
#endif
