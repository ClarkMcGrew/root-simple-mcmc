#include "../TSimpleMCMC.H"

#include "FakeLikelihood.H"

#include "TH1D.h"
#include "TPad.h"
#include "THStack.h"

#include <string>
#include <sstream>
#include <iomanip>

const int gBurninCycles = 5;
const int gBurninLength = 1000;

const int gChainCycles = 5;
const int gChainLength = 10000;

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
    like.Init(10000,100,10.0);

    THStack *dataStack = new THStack("dataStack", "A toy experiment");
    dataStack->Add(like.ToyData.DecayTag);
    dataStack->Add(like.ToyData.Separated);
    dataStack->Add(like.ToyData.Close);
    dataStack->Add(like.ToyData.VeryClose);

    like.DataVeryClose->Write();
    like.DataClose->Write();
    like.DataSeparated->Write();
    like.DataDecayTag->Write();
    
    like.SimulatedVeryClose->Write();
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
    nominalStack->Add(like.SimulatedVeryClose);
    dataStack->Draw();
    nominalStack->Draw("same");
    gPad->Print("FakeMCMC-nominal.png");
    
    // Set the number of dimensions for the proposal.
    proposal.SetDim(like.GetDim());
    
    proposal.SetGaussian(0,std::sqrt(1.0+like.MCTrueValues[0]));
    proposal.SetGaussian(1,std::sqrt(1.0+like.MCTrueValues[1]));

    // The number of dimensions in the point needs to agree with the number of
    // dimensions in the likelihood.  You can either hard code it, or do like
    // I'm doing here and have a likelihood method to return the number of
    // dimensions.
    Vector p(like.GetDim());

    // Set the starting point.
    for (std::size_t i=0; i<p.size(); ++i) {
        p[i] = like.MCTrueValues[i];
    }
    p[0] += gRandom->Gaus(0.0,std::sqrt(p[0]));
    p[1] += gRandom->Gaus(0.0,std::sqrt(p[1]));
    mcmc.Start(p,false);

    like.WriteSimulation(p,"initial");

    THStack *initialStack = new THStack("initialStack", "A toy experiment");
    initialStack->Add(like.SimulatedDecayTag);
    initialStack->Add(like.SimulatedSeparated);
    initialStack->Add(like.SimulatedClose);
    initialStack->Add(like.SimulatedVeryClose);
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
        int length = gBurninLength*(burnin+1)/gBurninCycles;
        std::cout << "Start new burnin phase ("<< length
                  << " steps)" << std::endl;
        for (int i=0; i<length; ++i) mcmc.Step(false);
        like.WriteSimulation(proposal.GetEstimatedCenter(),
                             burninName.str().c_str());
    
        THStack *simStack = new THStack("simStack", "A toy experiment");
        simStack->Add(like.SimulatedDecayTag);
        simStack->Add(like.SimulatedSeparated);
        simStack->Add(like.SimulatedClose);
        simStack->Add(like.SimulatedVeryClose);
        {
            p = proposal.GetEstimatedCenter();
            std::stringstream title;
            title << "Burn-in " << burnin + 1
                  <<" @" << std::fixed << std::setprecision(2);
            title << " (" << p[0];
            for (std::size_t i=1; i<p.size(); ++i) {
                title << ", " << p[i];
            }
            title << ")";
            dataStack->SetTitle(title.str().c_str());
        }
        dataStack->Draw();
        simStack->Draw("same");
        gPad->Print(("FakeMCMC-" + burninName.str() + ".png").c_str());
        like.Corrections.BackgroundShape->GetHistogram()->Draw("");
        gPad->Print(("FakeMCMC-" + burninName.str() + "-bkg.png").c_str());
        like.Corrections.SignalShape->GetHistogram()->Draw("");
        gPad->Print(("FakeMCMC-" + burninName.str() + "-sig.png").c_str());
    }


    for (int chain = 0; chain < gChainCycles; ++chain) {
        proposal.UpdateProposal();
    
        // Run the chain (now with output to the tree).
        std::cout << "Start chain " << chain << std::endl;
        for (int i=0; i<gChainLength; ++i) mcmc.Step();
        like.WriteSimulation(proposal.GetEstimatedCenter(),"midway");
        
        THStack *simStack = new THStack("simStack", "A toy experiment");
        simStack->Add(like.SimulatedDecayTag);
        simStack->Add(like.SimulatedSeparated);
        simStack->Add(like.SimulatedClose);
        simStack->Add(like.SimulatedVeryClose);
        p = proposal.GetEstimatedCenter();
        std::stringstream title;
        title << "Chain " << chain << " @"
              << std::fixed << std::setprecision(2);
        title << " (" << p[0];
        for (std::size_t i=1; i<p.size(); ++i) {
            title << ", " << p[i];
        }
        title << ")";
        dataStack->SetTitle(title.str().c_str());
        dataStack->Draw();
        simStack->Draw("same");

        std::stringstream name;
        name << "FakeMCMC-chain" << chain;
        gPad->Print((name.str() + ".png").c_str());
        like.Corrections.BackgroundShape->GetHistogram()->Draw("");
        gPad->Print((name.str() + "-bkg.png").c_str());
        like.Corrections.SignalShape->GetHistogram()->Draw("");
        gPad->Print((name.str() + "-sig.png").c_str());
    }
    
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
