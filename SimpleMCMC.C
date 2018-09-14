#include "TSimpleMCMC.H"

#include <sstream>

#include "TDummyLogLikelihood.H"

void SimpleMCMC(int trials, int maxEvaluations) {
    std::cout << "Simple MCMC Loaded" << std::endl;

#ifdef NO_OUTPUT
    TFile *outputFile = NULL;
    TTree *tree = NULL;
#else
    TFile *outputFile = new TFile("SimpleMCMC.root","recreate");
    TTree *tree = new TTree("SimpleMCMC","Tree of accepted points");
#endif
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

    mcmc.Start(p,false);

    // Burnin the chain (don't save the output)
    for (int i=0; i<10000+p.size()*p.size(); ++i) mcmc.Step(false);
    std::cout << "Finished burnin chain" << std::endl;

    // Burnin the chain (don't save the output)
    mcmc.GetProposeStep().UpdateProposal();
    for (int i=0; i<10000+10*p.size()*p.size(); ++i) mcmc.Step(false);
    std::cout << "Finished burnin chain" << std::endl;

    // Run the chain (now with output to the tree).
    mcmc.GetProposeStep().UpdateProposal();
    for (int i=0; i<trials; ++i) mcmc.Step();
    std::cout << "Finished with " << mcmc.GetLogLikelihoodCount() << " calls"
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
    SimpleMCMC(trials,maxEvaluations);
}
#endif
