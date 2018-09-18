#include <string>
#include <sstream>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TList.h>

/////////////////////////////////////////////////////////////////
// Read a tree written by the TSimpleMCMC and check that the constraints in
// the likehood are obeyed by the Markov chain.  This produces a profile
// histogram.  The last bin contains the average parameter value of all the
// bins (i.e. sum of the parameters divided by the number of parameters).  The
// values in the histogram should be compared with the values set in
// TConstrainedLikelihood.C::Init().
/////////////////////////////////////////////////////////////////
void ConstrainedCheck() {
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
    std::cout << "Entries: " << entries << std::endl;
    std::vector<double> *accepted = NULL;

    // Get an entry to collect information about the tree.
    inputTree->SetBranchAddress("Accepted",&accepted);
    inputTree->GetEntry(0);
    std::size_t dim = accepted->size();
    std::cout << "Dimensions: " << accepted->size() << std::endl;

    std::ostringstream title;
    title << "Mean value and sum of parameters for the \""
          << name << "\" tree";
    TProfile* mean = new TProfile("ConstrainedCheck",
                                  title.str().c_str(),
                                  dim+1, 0.5, dim+1.5,"S");
    
    // Calculate the averages.
    for(int e=0; e<entries; ++e) {
        inputTree->GetEntry(e);
        double sum = 0.0;
        for (int i = 0; i<accepted->size(); ++i) {
            mean->Fill(i+1.0,accepted->at(i));
            sum += accepted->at(i);
        }
        mean->Fill(accepted->size()+1.0,sum/accepted->size());
    }

    // Draw the histogram.
    mean->Draw();
    std::ostringstream label;
    label << "Parameter Number [bin " << accepted->size()+1
          << " is the #frac{1}{" << accepted->size()
          << "} #sum^{" << accepted->size() << "} P_{i} ]";
    mean->SetXTitle(label.str().c_str());
    mean->SetStats(false);
}
