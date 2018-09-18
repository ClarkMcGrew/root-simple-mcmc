#include "../TSimpleMCMC.H"

#include <sstream>

#include "TConstrainedLikelihood.H"

void Constrained(int trials, int maxEvaluations) {
    std::cout << "Simple MCMC Loaded" << std::endl;

#ifdef NO_OUTPUT
    TFile *outputFile = NULL;
    TTree *tree = NULL;
#else
    TFile *outputFile = new TFile("Constrained.root","recreate");
    TTree *tree = new TTree("Constrained","Tree of accepted points");
#endif
    TSimpleMCMC<TConstrainedLikelihood> mcmc(tree);
    TConstrainedLikelihood& like = mcmc.GetLogLikelihood();

    // Initialize the likelihood (if you need to)
    like.Init();
    
    // Set the number of dimensions for the proposal.
    mcmc.GetProposeStep().SetDim(like.GetDim());
    Vector p(like.GetDim());
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
// This let's the example compile directly.
int main(int argc, char **argv) {
    int trials = 100000;
    int maxEvaluations = -1;
    if (argc > 1) {
        std::istringstream input(argv[1]);
        input >> trials;
    }
    if (argc > 2) {
        std::istringstream input(argv[1]);
        input >> maxEvaluations;
    }
    Constrained(trials,maxEvaluations);
}
#endif
