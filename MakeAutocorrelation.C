/////////////////////////////////////////////////////////////////
// Read a tree written by the MCMC and calculate the covariance of the
// accepted points.  The tree is expected to have a "std::vector<double>"
// branch named "Accepted" that contains all of the accepted points.  There's
// no error checking for the branch, so it will usually crash if the branch
// doesn't exist (watch the ROOT output)!  This is run:
//
//  root input.root MakeAutocorrelation.C
//
// The output is saved in a file named autocorrelation.root which contains
// (lots of) histograms:
//
//   avgAutocorrelation (TProfile) -- A profile with the average
//       autocorrelation for all dimensions vs the lag.
//
//   meanValues (TProfile) -- A profile with the mean value for each
//       dimension.
//
//   autoCorrD<n> (TH1) -- A histogram for each dimension containing the
//       autocorrelation vs the lag for that dimension.
//
/////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TList.h>
#include <TTree.h>
#include <TKey.h>

void MakeAutocorrelation() {
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
    int dim = accepted->size();

    // The history of the points
    const int depth = 30000;
    std::vector<std::array<float,depth>> ringBuffer(dim);
    int nextBuffer = -1;

    // Open the file to save the correlations.  This overwrites any existing
    // file.
    TFile outputFile("autocorrelation.root","recreate");

    // Find the maximum lag to use for this data set.
    int maxLag = entries - std::sqrt(entries);
    maxLag = std::min(maxLag,depth);
    int bins = std::min(100,maxLag);

    // Create the summary histograms.
    TProfile *meanValues = new TProfile("meanValues",
                                        "Mean and RMS",
                                        dim, 0.0, dim,"s");
    TProfile *avgCorr = new TProfile("avgAutocorrelation",
                                     "Average Autocorrelation",
                                     bins, 0.0, maxLag,"S");

    // Create the autocorrelation histograms for each dimension.
    std::vector<std::pair<TH1*,TH1*>> autoCorr(dim);
    for (int i = 0; i<dim; ++i) {
        std::ostringstream name;
        name << "D" << std::setw(3) << std::setfill('0') << i;
        autoCorr[i].first = new TH1F(("autoCorr" + name.str()).c_str(),
                                     ("Autocorrelation: "+name.str()).c_str(),
                                     bins, 0.0, maxLag);
        autoCorr[i].second = new TH1F(("autoCount" + name.str()).c_str(),
                                      "Autocorrelation Count",
                                      bins, 0.0, maxLag);
    }

    // We don't need to calculate the auto correlation kernel for every lag.
    // Figure out which ones we need to calculate (this is for efficiency).
    int lagStep = 0.5*maxLag/bins;
    if (lagStep < 1) lagStep = 1;

    // Loop over the entries and fill the auto correlation kernel.  Only look
    // at the end of the file.  The number of trials is limited so that this
    // goes faster.
    int trials = 2*maxLag + std::sqrt(entries);
    trials = std::min(trials,entries);
    int fills = 0;            // Track total entries added to the ring buffer.
    for (int entry = entries-trials; entry < entries; ++entry) {
        inputTree->GetEntry(entry);
        nextBuffer = (++nextBuffer)%depth;
        ++fills;
        if (entry%1000 == 0) std::cout << "entry " << entry << std::endl;
        for (int i = 0; i<dim; ++i) {
            double val = accepted->at(i);
            meanValues->Fill(i+0.5,val);
            ringBuffer[i][nextBuffer] = val;
            for (int lag=1; lag < maxLag; lag += lagStep) {
                if (fills <= lag) break;
                int lagBuffer = (nextBuffer+depth-lag)%depth;
                double lVal = ringBuffer[i][lagBuffer];
                autoCorr[i].first->Fill(lag+0.5,lVal*val);
                autoCorr[i].second->Fill(lag+0.5,1.0);
            }
        }
    }

    // Turn the autocorrelation kernels into the "Pearson" autocorrelation.
    // For a "Gaussian" like distribution this will be in [-1,+1], but for
    // pathelogical cases (like the ever popular generalized Rosenbrock
    // distribution) it can take on almost any value.
    for (int i = 0; i<dim; ++i) {
        double mean = meanValues->GetBinContent(i+1);
        double err = meanValues->GetBinError(i+1);
        for (int j = 0; j < bins; ++j) {
            double c = autoCorr[i].first->GetBinCenter(j+1);
            double v = autoCorr[i].first->GetBinContent(j+1);
            double ve = autoCorr[i].first->GetBinError(j+1);
            double e = autoCorr[i].second->GetBinContent(j+1);
            double a = (v/e - mean*mean)/(err*err);
            // std::cout << i << " " << j
            //           << " " << v << " " << e
            //           << " " << mean << " " << err
            //           << " " << a << std::endl;
            autoCorr[i].first->SetBinContent(j+1,a);
            autoCorr[i].first->SetBinError(j+1,ve/e);
            avgCorr->Fill(c,a);
        }
    }

    // Make the plots and write the output.
    avgCorr->Draw("e1");
    for (int i=0; i<dim; ++i) {
        autoCorr[i].first->Draw("H,e2,same");
        autoCorr[i].first->Write();
    }
    avgCorr->SetLineWidth(2);
    avgCorr->SetLineColor(kRed);
    avgCorr->Draw("e1,Same");
    avgCorr->Write();
    meanValues->Write();
    gPad->Print("averageAutocorrelation.pdf");
    gPad->Print("averageAutocorrelation.png");
}
