#include "TSimpleHMC.H"

#include <TMatrixD.h>
#include <TVectorD.h>

#include <sstream>

// A dummy log likelihood for testing.  This is an example, but don't
// slavishly copy it (or you will be sorry!).
class TDummyLogLikelihood {
public:
    // Determine the number of dimensions.  This is where the dimensions are
    // defined, and everything else uses it.
    std::size_t GetDim() const {return 4;}

    // Calculate the log(likelihood).  The dummy likelihood is a Gaussian
    // (with covariance) centered at zero.  The covariance is set in Init()
    // (below).  
    double operator()(const Vector& point)  const {
        double logLikelihood = 0.0;

        for (std::size_t i = 0; i<GetDim(); ++i) {
            for (std::size_t j = 0; j<GetDim(); ++j) {
                logLikelihood -= 0.5*point[i]*Error(j,i)*point[j];
            }
        }

        return logLikelihood;
    }

    // Note that this needs to be the grad(log(Likelihood)).  
    bool operator() (Vector& g, const Vector& p) {
        for (int i=0; i<p.size(); ++i) {
            g[i] = 0.0;
            for (int j=0; j<p.size(); ++j) {
                g[i] -= Error(i,j)*p[j];
            }
        }
        return true;
    }

    void Init() {
        Covariance.ResizeTo(GetDim(),GetDim());
        Error.ResizeTo(GetDim(),GetDim());

        // Set the sigma for each variable.
        for (std::size_t i = 0; i<GetDim(); ++i) {
            double sigma = 1.0;
            // double sigma = 1.0*i + 1.0;
            Covariance(i,i) = sigma*sigma;
        }

        for (std::size_t i = 0; i<GetDim(); ++i) {
            for (std::size_t j = i+1; j<GetDim(); ++j) {
                double sig1 = std::sqrt(Covariance(i,i));
                double sig2 = std::sqrt(Covariance(j,j));
                // Now give some correlations to the likelihood.  (Uncomment
                // the one you want to try).

                // Choose a random correlation
#ifdef RANDOM_CORRELATION
                Covariance(i,j) = gRandom->Uniform(-0.999,0.999)*sig1*sig2;
#endif
                
                // Choose no correlation
#define NO_CORRELATION
#ifdef NO_CORRELATION
                Covariance(i,j) = 0.0;
#endif

                // Choose a correlation based on the variables.  Neighbors are 
                // not correlated, but there is more correlation as the
                // variables are further apart.
#ifdef VERY_CORRELATED
                if (i+j==GetDim()-1) {
                    Covariance(i,j) = 0.999*sig1*sig2*(j - i)/(GetDim()-1.0);
                }
#endif
                
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

        Covariance.Print();
        Error.Print();
    }

    static TMatrixD Covariance;
    static TMatrixD Error;
};
TMatrixD TDummyLogLikelihood::Covariance;
TMatrixD TDummyLogLikelihood::Error;

void SimpleHMC(int maxEvals=-1) {
    std::cout << "Simple HMC Loaded" << std::endl;

#ifdef NO_OUTPUT
    TFile *outputFile = NULL;
    TTree *tree = NULL;
#else
    TFile *outputFile = new TFile("SimpleHMC.root","recreate");
    TTree *tree = new TTree("SimpleHMC","Tree of accepted pqoints");
#endif
    TSimpleHMC<TDummyLogLikelihood,TDummyLogLikelihood> hmc(tree);
    // TSimpleHMC<TDummyLogLikelihood> hmc(tree);
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

    hmc.Start(p,false);

    // Burnin the chain (don't save the output)
    int trials = 10000;
    for (int i=0; i<trials; ++i) {
        if (i%1000 == 0) {
            std::cout << i << " " << hmc.GetPotentialCount()
                      << std::endl;
        }
        hmc.Step(true);
        if (maxEvals > 0 && hmc.GetPotentialCount() > maxEvals) break;
    }
    std::cout << "Finished " << trials
              << " requested trials with calls " << hmc.GetPotentialCount()
              << std::endl;

    if (tree) tree->Write();
    if (outputFile) delete outputFile;
}

#ifdef MAIN_PROGRAM
// This let's the example compile directly.  To compile it, use the compile.sh
// script and then run it using ./a.out which will produce a file name
// "SimpleHMC.root"
int main(int argc, char **argv) {
    int maxEvaluations = -1;
    if (argc > 1) {
        std::istringstream input(argv[1]);
        input >> maxEvaluations;
    }
    SimpleHMC(maxEvaluations);
}
#endif
