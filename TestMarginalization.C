// This should be run on a THardLogLikelihood with a smallish number of
// dimensions (or the plots become unreasonable).

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TList.h>

void TestMarginalization() {
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

    if (dim > 10) {
        std::cout << "I think you're making a mistake" << std::endl;
        return;
    }

    // Define the range for the plots... This could be A LOT more intelligent.
    std::vector<std::pair<double,double>> range;

    for (int i= 0; i< dim; ++i) {
        range.push_back(std::make_pair(accepted->at(i),accepted->at(i)));
    }

    for (int entry = 0; entry < entries; ++entry) {
        inputTree->GetEntry(entry);
        for (int i= 0; i< dim; ++i) {
            range[i].first = std::min(range[i].first,accepted->at(i));
            range[i].second = std::max(range[i].second,accepted->at(i));
        }
        entry += 0.001*entries;
    }

    std::cout << "Dimensions: " << dim << std::endl;
    std::vector<TH1F*> marginalized;
    for (int i= 0; i< dim; ++i) {
        std::ostringstream name;
        std::ostringstream title;
        name << "marginalized" << i;
        title << "Marginalization: Dimension " << i;
        TH1F* hist = new TH1F(name.str().c_str(),title.str().c_str(),
                              100, -1.5, 1.5);
        marginalized.push_back(hist);
        std::cout << range[i].first << " " << range[i].second << std::endl;
    }

    std::vector<TH2F*> correlations;
    for (int i= 0; i< dim; ++i) {
        for (int j= 0; j< dim; ++j) {
            std::ostringstream name;
            std::ostringstream title;
            name << "correlation" << i << "_" << j;
            title << "Dimensions " << j << " vs " << i;
            double ri = 0.05*(range[i].second-range[i].first);
            double rj = 0.05*(range[j].second-range[j].first);
            TH2F* hist = new TH2F(name.str().c_str(),title.str().c_str(),
                                  50, range[i].first-ri, range[i].second+ri,
                                  50, range[j].first-rj, range[j].second+rj);
            correlations.push_back(hist);
        }
    }

    for (int entry = 0; entry < entries; ++entry) {
        inputTree->GetEntry(entry);
        for (int i = 0; i<dim; ++i) {
            marginalized[i]->Fill(accepted->at(i));
            for (int j = 0; j<dim; ++j) {
                correlations[dim*i+j]->Fill(accepted->at(i),accepted->at(j));
            }
        }
    }

    marginalized[0]->SetTitle("Marginalized Rosenbrock (all dimensions)");
    marginalized[0]->Draw("C");
    marginalized[0]->SetLineWidth(3.0);
    for (int i = 1; i<dim; ++i) {
        marginalized[i]->Draw("same,C");
    }
    marginalized[dim-1]->SetLineWidth(3.0);
    marginalized[dim-1]->SetLineColor(kRed);
    gPad->Print("marginalized.png");

    TCanvas *zoned = new TCanvas("zoned","multipads",1000,1000);
    gStyle->SetOptStat(0);
    zoned->Divide(dim,dim,0,0);
    for (int i = 0; i<dim; ++i) {
        for (int j = 0; j<dim; ++j) {
            zoned->cd(dim*i+j+1);
            correlations[dim*i+j]->Draw("colz");
        }
    }
    zoned->cd(0);
    gPad->Print("correlations.png");
}
