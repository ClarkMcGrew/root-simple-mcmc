#include "TSimpleMCMC.H"

// A dummy log likelihood for testing.  This is an example, but don't
// slavishly copy it (or you will be sorry!).
class TDummyLogLikelihood {
public:
    // Determine the number of dimensions.  This is where the dimensions are
    // defined, and everything else uses it.
    std::size_t GetDim() const {return 10;}

    // Calculate the likelihood.  The dummy likelihood is a Gaussian (with
    // covariance) centered at zero.  The covariance is set in Init() (below).
    double operator()(const Vector& point)  const {
        double logLikelihood = 0.0;

        for (std::size_t i = 0; i<GetDim(); ++i) {
            for (std::size_t j = 0; j<GetDim(); ++j) {
                logLikelihood -= 0.5*point[i]*Error(j,i)*point[j];
            }
        }

        return logLikelihood;
    }

    void Init() {
        Covariance.ResizeTo(GetDim(),GetDim());
        Error.ResizeTo(GetDim(),GetDim());

        // Set the sigma for each variable.
        for (std::size_t i = 0; i<GetDim(); ++i) {
            // double sigma = 1.0;
            double sigma = 1.0*i + 0.1;
            Covariance(i,i) = sigma*sigma;
        }

        for (std::size_t i = 0; i<GetDim(); ++i) {
            for (std::size_t j = i+1; j<GetDim(); ++j) {
                double sig1 = std::sqrt(Covariance(i,i));
                double sig2 = std::sqrt(Covariance(j,j));
                // Now give some correlations to the likelihood.  (Uncomment
                // the one you want to try).
                Covariance(i,j) = gRandom->Uniform(0.0,0.9)*sig1*sig2;
                // Covariance(i,j) = 0.90*sig1*sig2;
                // Covariance(i,j) = 0.0;
                Covariance(j,i) = Covariance(i,j);
            }
        }

        // Make sure the covariance is positive definite.
        do {
            TVectorD eigenValues(GetDim());
            Covariance.EigenVectors(eigenValues);
            bool positiveDefinite = true;
            for (std::size_t i = 0; i<GetDim(); ++i) {
                if (eigenValues(i)<0.0) {
                    positiveDefinite = false;
                }
            }
            if (positiveDefinite) break;
            for (std::size_t i = 0; i<GetDim(); ++i) {
                for (std::size_t j = i+1; j<GetDim(); ++j) {
                    Covariance(i,j) = 0.9*Covariance(i,j);
                    Covariance(j,i) = Covariance(i,j);
                }
            }
        } while (true);

        Error = Covariance;
        Error.Invert();

        if (GetDim() < 5) {
            Covariance.Print();
            Error.Print();
        }
    }

    TMatrixD Covariance;
    TMatrixD Error;
};

void SimpleMCMC() {
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
    TProposeAdaptiveStep& proposal = mcmc.GetProposeStep();

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

    mcmc.GetProposeStep().ResetProposal();

    // Burnin the chain (don't save the output)
    for (int i=0; i<10000; ++i) mcmc.Step(false);

    std::cout << "Finished burnin " << std::endl;
    
    mcmc.GetProposeStep().ResetProposal();

    // Run the chain (now with output to the tree).
    for (int i=0; i<1000000; ++i) mcmc.Step();

    if (tree) tree->Write();
    if (outputFile) delete outputFile;
}

#ifdef MAIN_PROGRAM
// This let's the example compile directly.  To compile it, use the compile.sh
// script and then run it using ./a.out which will produce a file name
// "SimpleMCMC.root"
int main(int argc, char **argv) {
    SimpleMCMC();
}
#endif
