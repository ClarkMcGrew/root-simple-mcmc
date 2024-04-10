#include "TSimpleHMC.H"

#include <TMatrixD.h>
#include <TVectorD.h>

#include <sstream>

// A dummy log likelihood for testing.  This is a test bed for HMC where the
// gradient of the logLikelihood is being incorrectly calculated.  NOTICE: THE
// GRADIENT IS INTENTIONALLY WRONG.
class TDummyLogLikelihood {
public:
    // Determine the number of dimensions.  This is where the dimensions are
    // defined, and everything else uses it.
    std::size_t GetDim() const {return 50;}

    // Calculate the log(likelihood).  The dummy likelihood is a Gaussian
    // (with covariance) centered at zero.  The covariance is set in Init()
    // (below).
    double operator()(const sMCMC::Vector& point)  const {
        double logLikelihood = 0.0;

        for (std::size_t i = 0; i<GetDim(); ++i) {
            for (std::size_t j = 0; j<GetDim(); ++j) {
                logLikelihood -= 0.5*point[i]*Error(j,i)*point[j];
            }
        }

        return logLikelihood;
    }

    // Note that this needs to be the grad(log(Likelihood)).
    bool operator() (sMCMC::Vector& g, const sMCMC::Vector& p) {
        for (int i=0; i<p.size(); ++i) {
            g[i] = 0.0;
            for (int j=0; j<p.size(); ++j) {
                g[i] -= GradientError(i,j)*p[j];
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
#ifdef NO_CORRELATION
                Covariance(i,j) = 0.0;
#endif

                // Choose a correlation based on the variables.  Neighbors are
                // not correlated, but there is more correlation as the
                // variables are further apart.
#define VERY_CORRELATED
#ifdef VERY_CORRELATED
                if (i+j==GetDim()-1) {
                    Covariance(i,j) = 0.900*sig1*sig2*(j - i)/(GetDim()-1.0);
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

        // Make a gradient based on the covariance.  The gradient is
        // intentionally WRONG.
        GradientCovariance.ResizeTo(GetDim(),GetDim());
        GradientError.ResizeTo(GetDim(),GetDim());
        for (std::size_t i = 0; i<GetDim(); ++i) {
            for (std::size_t j = i; j<GetDim(); ++j) {
                double sig1 = std::sqrt(Covariance(i,i));
                double sig2 = std::sqrt(Covariance(j,j));
                double r = Covariance(i,j);
                if (i == j) {
                    double s = 0.1;
                    double e = gRandom->Gaus(1.0,s);
                    while (e < 0.3) e = gRandom->Gaus(1.0,s);
                    r = r*e;
                }
                else {
                    r = r + gRandom->Gaus(0.0,0.3)*sig1*sig2;
                }
                GradientCovariance(i,j) = r;
                GradientCovariance(j,i) = r;
            }
        }

        // Make sure the gradient covariance is positive definite.
        do {
            TVectorD eigenValues(GetDim());
            GradientCovariance.EigenVectors(eigenValues);
            bool positiveDefinite = true;
            for (std::size_t i = 0; i<GetDim(); ++i) {
                if (eigenValues(i)<0.0) {
                    positiveDefinite = false;
                }
            }
            if (positiveDefinite) break;
            for (std::size_t i = 0; i<GetDim(); ++i) {
                for (std::size_t j = i+1; j<GetDim(); ++j) {
                    GradientCovariance(i,j) = 0.9*GradientCovariance(i,j);
                    GradientCovariance(j,i) = GradientCovariance(i,j);
                }
            }
        } while (true);

        GradientError = GradientCovariance;
        GradientError.Invert();

        Covariance.Print();
        GradientCovariance.Print();
    }

    static TMatrixD Covariance;
    static TMatrixD Error;
    static TMatrixD GradientCovariance;
    static TMatrixD GradientError;
};
TMatrixD TDummyLogLikelihood::Covariance;
TMatrixD TDummyLogLikelihood::Error;
TMatrixD TDummyLogLikelihood::GradientCovariance;
TMatrixD TDummyLogLikelihood::GradientError;

void BadGrad(int maxEvals=-1) {
    std::cout << "Simple HMC Loaded" << std::endl;

#ifdef NO_OUTPUT
    TFile *outputFile = NULL;
    TTree *tree = NULL;
#else
    TFile *outputFile = new TFile("badGrad.root","recreate");
    TTree *tree = new TTree("BadGrad","Tree of accepted pqoints");
#endif
    sMCMC::TSimpleHMC<TDummyLogLikelihood,TDummyLogLikelihood> hmc(tree);
    // TSimpleHMC<TDummyLogLikelihood> hmc(tree);
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

    hmc.Start(p,false);

    // Burnin the chain (don't save the output)
    int trials = 100000;
    for (int i=0; i<trials; ++i) {
        if (i%1000 == 0) {
            std::cout << i << " " << hmc.GetPotentialCount()
                      << " " << hmc.GetGradientCount()
                      << std::endl;
        }
        hmc.Step(true);
        if (maxEvals > 0 && hmc.GetPotentialCount() > maxEvals) break;
    }
    std::cout << "Finished " << trials
              << " requested trials with calls " << hmc.GetPotentialCount()
              << " + " << hmc.GetGradientCount()
              << " " << 1.0*hmc.GetGradientCount()/hmc.GetPotentialCount()
              << std::endl;

    if (tree) tree->Write();
    if (outputFile) delete outputFile;
}

#ifdef MAIN_PROGRAM
// This let's the example compile directly.  To compile it, use the compile.sh
// script and then run it using ./a.out which will produce a file name
// "badGrad.root"
int main(int argc, char **argv) {
    int maxEvaluations = -1;
    if (argc > 1) {
        std::istringstream input(argv[1]);
        input >> maxEvaluations;
    }
    BadGrad(maxEvaluations);
}
#endif
