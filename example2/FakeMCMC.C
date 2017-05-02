#include "../TSimpleMCMC.H"

#include "FakeLikelihood.H"

#include "TH1D.h"
#include "TPad.h"
#include "THStack.h"

#include <string>
#include <sstream>
#include <iomanip>

const int gBurninCycles = 3;
const int gBurninLength = 3000;

const int gCycles = 4;
const int gChainLength = 50000;

void FakeMCMC() {
    std::cout << "Fake Likelihood MCMC Loaded" << std::endl;
    gRandom->SetSeed();

#ifdef NO_OUTPUT
    TFile *outputFile = NULL;
    TTree *tree = NULL;
#else
    TFile *outputFile = new TFile("FakeMCMC.root","recreate");
    TTree *tree = new TTree("MCMC","Tree of accepted points");
#endif
    TSimpleMCMC<FakeLikelihood> mcmc(tree);
    FakeLikelihood& like = mcmc.GetLogLikelihood();
    TProposeAdaptiveStep& proposal = mcmc.GetProposeStep();

    // Initialize the likelihood (if you need to).  The dummy likelihood
    // setups a covariance to make the PDF more interesting.
    like.Init(1000,1000,10.0);

    THStack *dataStack = new THStack("dataStack", "A toy experiment");
    dataStack->Add(like.ToyData.DecayTag);
    dataStack->Add(like.ToyData.Separated);
    dataStack->Add(like.ToyData.Close);

    like.DataClose->Write();
    like.DataSeparated->Write();
    like.DataDecayTag->Write();
    
    like.SimulatedClose->Write();
    like.SimulatedSeparated->Write();
    like.SimulatedDecayTag->Write();
    
    // Find the total number of data events.
    double dataEvents = like.DataDecayTag->Integral();
    dataEvents += like.DataSeparated->Integral();
    dataEvents += like.DataClose->Integral();

    THStack *nominalStack = new THStack("nominalStack", "A toy experiment");
    nominalStack->Add(like.SimulatedDecayTag);
    nominalStack->Add(like.SimulatedSeparated);
    nominalStack->Add(like.SimulatedClose);
    dataStack->Draw();
    nominalStack->Draw("same");
    gPad->Print("FakeMCMC-nominal.png");
    
    // Set the number of dimensions for the proposal.
    proposal.SetDim(like.GetDim());
    
    // The number of dimensions in the point needs to agree with the number of
    // dimensions in the likelihood.  You can either hard code it, or do like
    // I'm doing here and have a likelihood method to return the number of
    // dimensions.
    Vector p(like.GetDim());
    for (std::size_t i=0; i<p.size(); ++i) {
        p[i] = like.MCTrueValues[i] + gRandom->Uniform(-1.0,1.0);
    }
 
    mcmc.Start(p,false);

    like.WriteSimulation(p,"initial");

    THStack *initialStack = new THStack("initialStack", "A toy experiment");
    initialStack->Add(like.SimulatedDecayTag);
    initialStack->Add(like.SimulatedSeparated);
    initialStack->Add(like.SimulatedClose);
    {
        std::stringstream title;
        title << "Initial @" << std::fixed << std::setprecision(2);
        title << " (" << p[0];
        for (std::size_t i=1; i<p.size(); ++i) {
            title << ", " << p[i];
        }
        title << ")";
        dataStack->SetTitle(title.str().c_str());
    }
    dataStack->Draw();
    initialStack->Draw("same");
    gPad->Print("FakeMCMC-initial.png");

    // First burnin of the chain (don't save the output).  This is looking for
    // the best fit point.
    for (int burnin = 0; burnin<gBurninCycles; ++burnin) {
        proposal.ResetProposal();

        std::stringstream burninName;
        burninName << "burnin" << burnin+1;
        for (int i=0; i<gBurninLength; ++i) mcmc.Step(false);
        like.WriteSimulation(proposal.GetEstimatedCenter(),
                             burninName.str().c_str());
    
        THStack *burninStack = new THStack("burninStack", "A toy experiment");
        burninStack->Add(like.SimulatedDecayTag);
        burninStack->Add(like.SimulatedSeparated);
        burninStack->Add(like.SimulatedClose);
        {
            p = proposal.GetEstimatedCenter();
            std::stringstream title;
            title << "Burnin @" << std::fixed << std::setprecision(2);
            title << " (" << p[0];
            for (std::size_t i=1; i<p.size(); ++i) {
                title << ", " << p[i];
            }
            title << ")";
            dataStack->SetTitle(title.str().c_str());
        }
        dataStack->Draw();
        burninStack->Draw("same");
        gPad->Print(("FakeMCMC-" + burninName.str() + ".png").c_str());
    }

    proposal.UpdateProposal();
    
    // Run the chain (now with output to the tree).
    std::cout << "Start the chain" << std::endl;
    for (int i=0; i<gChainLength; ++i) mcmc.Step();
    like.WriteSimulation(proposal.GetEstimatedCenter(),"midway");
    
    THStack *midwayStack = new THStack("midwayStack", "A toy experiment");
    midwayStack->Add(like.SimulatedDecayTag);
    midwayStack->Add(like.SimulatedSeparated);
    midwayStack->Add(like.SimulatedClose);
    {
        p = proposal.GetEstimatedCenter();
        std::stringstream title;
        title << "Midway Chain @" << std::fixed << std::setprecision(2);
        title << " (" << p[0];
        for (std::size_t i=1; i<p.size(); ++i) {
            title << ", " << p[i];
        }
        title << ")";
        dataStack->SetTitle(title.str().c_str());
    }
    dataStack->Draw();
    midwayStack->Draw("same");
    gPad->Print("FakeMCMC-midway.png");
    
    for (int i=0; i<gChainLength; ++i) mcmc.Step();
    like.WriteSimulation(proposal.GetEstimatedCenter(),"final");
    
    THStack *finalStack = new THStack("finalStack", "A toy experiment");
    finalStack->Add(like.SimulatedDecayTag);
    finalStack->Add(like.SimulatedSeparated);
    finalStack->Add(like.SimulatedClose);
    {
        p = proposal.GetEstimatedCenter();
        std::stringstream title;
        title << "Final Chain @" << std::fixed << std::setprecision(2);
        title << " (" << p[0];
        for (std::size_t i=1; i<p.size(); ++i) {
            title << ", " << p[i];
        }
        title << ")";
        dataStack->SetTitle(title.str().c_str());
    }
    dataStack->Draw();
    finalStack->Draw("same");
    gPad->Print("FakeMCMC-final.png");
    
    if (tree) tree->Write();
    if (outputFile) delete outputFile;
}

#ifdef MAIN_PROGRAM
// This let's the example compile directly.  To compile it, use the compile.sh
// script and then run it using ./a.out which will produce a file name
// "FakeMCMC.root"
int main(int argc, char **argv) {
    FakeMCMC();
}
#endif
