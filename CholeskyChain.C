#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TMatrixD.h>
#include <TDecompChol.h>

// Take an input file containing the mean and covariance of a distribution,
// and then use Cholesky decomposition to write a fake MCMC chain.  This
// expects a file produced by MakeCovariance.C, but any file containing TH2
// histogram of covariance values (named "AcceptedCovariance"), and a TH1
// histogram of mean values (named "AcceptedMean") should work.  The output is
// in a file named "Cholesky.root".  This is run:
//
// root covariance.root CholeskyChain.C
//
void CholeskyChain(int trials=100000) {
    // Get the covariance and mean values from the input file.
    TH2 *covarianceHist = (TH2*) gFile->Get("AcceptedCovariance");
    TH1 *meanHist = (TH1*) gFile->Get("AcceptedMean");

    // Open the output file.
    TFile output("Cholesky.root","recreate");

    // Make the mean vector and covariance matrix
    int dim = covarianceHist->GetNbinsX();
    std::cout << "Dimensions: " << dim << std::endl;

    std::vector<double> mean(dim);
    TMatrixD covariance(dim,dim);
    for (int i=0; i<dim; ++i) {
        mean[i] = meanHist->GetBinContent(i+1);
        for (int j=0; j<dim; ++j) {
            covariance(i,j) = covarianceHist->GetBinContent(i+1,j+1);
        }
    }
    
    TMatrixD decomposition(dim,dim);
    TDecompChol cholesky(covariance);
    if (!cholesky.Decompose()) {
        std::cout << "Decomposition of the covariance has failed"
                  << std::endl;
        std::exit(1);
    }
    decomposition = cholesky.GetU();

    // Create the output tree.
    TTree *tree = new TTree("Cholesky","Tree of accepted points");
    std::vector<double> accepted(dim);
    tree->Branch("Accepted",&accepted);

    for (int trial = 0; trial < trials; ++trial) {
        std::copy(mean.begin(), mean.end(), accepted.begin());
        for (int i = 0; i < dim; ++i) {
            double r = gRandom->Gaus(0.0,1.0);
            for (std::size_t j = 0; j < dim; ++j) {
                accepted[j] += r*decomposition(i,j);
            }
        }
        tree->Fill();
    }    

    // Close the tree.
    tree->Write();
    output.Close();
}
