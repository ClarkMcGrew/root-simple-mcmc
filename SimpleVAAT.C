#include "TSimpleMCMC.H"
#include "TProposeVAATStep.H"

#include <sstream>
#include <string>

#include "TDummyLogLikelihood.H"

void SimpleVAAT(int cycles, int steps,
                 std::string outputName) {

    std::cout << "Simple VAAT Loaded" << std::endl;

    std::unique_ptr<TFile> outputFile;
    std::unique_ptr<TTree> tree;
    if (not outputName.empty()) {
        outputFile = std::make_unique<TFile>(outputName.c_str(),"recreate");
        tree = std::make_unique<TTree>("SimpleVAAT","Tree of accepted points");
    }

    sMCMC::TSimpleMCMC<TDummyLogLikelihood,sMCMC::TProposeVAATStep> mcmc(tree.get());
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
    sMCMC::Vector p(like.GetDim());
    for (std::size_t i=0; i<p.size(); ++i) p[i] = gRandom->Uniform(-1.0,1.0);

    mcmc.Start(p,false);
    mcmc.GetProposeStep().UpdateProposal();

    // Run the chain (with output to the tree).
    int verbosity = std::max(1,steps*cycles/100);
    int trial = 0;
    for (int cycle = 0; cycle < cycles; ++cycle) {
        for (int step=0; step<steps; ++step) {
            ++trial;
            mcmc.Step();
            if (trial%verbosity == 0) {
                std::cout << "Trial " << cycle+1 << ":" << step+1
                          << " Total: " << trial << "/" << cycles*steps
                          << " Acceptance: " << mcmc.GetProposeStep().GetAcceptance()
                          << " (" << mcmc.GetProposeStep().GetSuccesses()
                          << "/" << mcmc.GetProposeStep().GetTrials() << ")"
                          << " Sigma: " << mcmc.GetProposeStep().GetSigma()
                          << std::endl;
            }
        }
    }

    if (tree) tree->Write();

    std::cout << "Exit" << std::endl;
}

#ifdef MAIN_PROGRAM
// This let's the example compile directly.  To compile it, use the compile.sh
// script and then run it using ./a.out which will produce a file name
// "SimpleAHMC.root"
int main(int argc, char **argv) {
    int cycles = 10;
    int steps = 1000;
    std::string outputName("SimpleMCMC.root");

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

    SimpleVAAT(cycles,steps,outputName);
}
#endif
