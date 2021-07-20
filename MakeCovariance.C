#include <string>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TList.h>

/////////////////////////////////////////////////////////////////
// Read a tree written by the TSimpleMCMC and calculate the covariance of the
// accepted points.  This is run:
//
//  root input.root MakeCovariance.C
//
// The output is saved in a file named covariance.root which contains
// histograms:
//
// AcceptedCovariance (TH2D) -- A 2 D histogram with the bin values
//     representing the covariance (bin 1 is for Accepted[0]).
//
// AcceptedMean (TProfile) -- A 1 D profile histogram with the bin values
//     representing the means (bin 1 is Accepted[0]).  The errors on the
//     profile points are the marginalized spread for each parameter
//     (calculated using the TProfile "spread" option).
//
/////////////////////////////////////////////////////////////////
void MakeCovariance() {
    // Find the tree in the file.
    TList *list = gFile->GetListOfKeys();
    TIter iter(list->MakeIterator());
    std::string name;
    while (TObject *obj = iter()) {
        TKey *key = (TKey*)obj;
        if (std::string(key->GetClassName()) != "TTree") continue;
        name = key->GetName();
    }

    // Get the tree out of the file.
    TTree *inputTree = (TTree*) gFile->Get(name.c_str());
    std::cout << "Input Tree Name: " << inputTree->GetName() << std::endl;
    int entries = inputTree->GetEntries();
    std::cout << "Entries: " << entries << std::endl;\
    std::vector<double> *accepted = NULL;

    // Get an entry to collect information about the tree.
    inputTree->SetBranchAddress("Accepted",&accepted);
    inputTree->GetEntry(0);
    std::size_t dim = accepted->size();
    std::cout << "Dimensions: " << accepted->size() << std::endl;

    // Create an output histogram for the covariance.
    TFile output("covariance.root","recreate");
    TH2* covariance = new TH2D("AcceptedCovariance",
                               "Covariance of the Accepted Points",
                               dim, 0, dim,
                               dim, 0, dim);
    TProfile* mean = new TProfile("AcceptedMean",
                                  "Mean value of the Accepted Points",
                                  dim, 0, dim,"S");

    // Calculate the average and covariance.
    std::vector<double> avg(dim);
    for(int e=0; e<entries; ++e) {
        inputTree->GetEntry(e);
        if (e%1000 == 0) std::cout << "entry " << e << std::endl;
        for (int i = 0; i<accepted->size(); ++i) {
            mean->Fill(i+0.1,accepted->at(i));
            avg[i] += accepted->at(i);
            for (int j = 0; j<accepted->size(); ++j) {
                covariance->Fill(i+0.1,j+0.1,
                                 accepted->at(i)*accepted->at(j));
            }
        }
    }

    for (int i=0; i<dim; ++i) avg[i] /= entries;
    double minCov = 0.0;
    double maxCov = 0.0;
    for (int i=0; i<dim; ++i) {
        for (int j=0; j<dim; ++j) {
            double v = covariance->GetBinContent(i+1,j+1);
            // This assumes entries is very large...
            v = v/entries - avg[i]*avg[j];
            covariance->SetBinContent(i+1,j+1,v);
            minCov = std::min(minCov,v);
            maxCov = std::max(maxCov,v);
        }
    }

    // Save it to a file.
    covariance->SetMinimum(minCov);
    covariance->SetMaximum(maxCov);
    covariance->SetStats(false);
    covariance->Draw("colz");
    gPad->Print("AcceptedCovariance.png");

    covariance->Write();
    mean->Write();
}
