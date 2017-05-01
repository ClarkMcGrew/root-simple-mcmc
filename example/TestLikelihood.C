#include "FakeLikelihood.H"

#include <TH1D.h>
#include <TPad.h>
#include <THStack.h>

#include <vector>
#include <iostream>

void TestLikelihood() {
    std::cout << "Start test" << std::endl;
    gRandom->SetSeed();
    
    FakeLikelihood likelihood;

    likelihood.Init(10000,10000,1.0);

    THStack *dataStack = new THStack("dataStack", "A toy experiment");
    dataStack->Add(likelihood.DataDecayTag);
    dataStack->Add(likelihood.DataSeparated);
    dataStack->Add(likelihood.DataClose);
    
    std::vector<double> params(likelihood.GetDim());

    likelihood.FillHistograms(params);

    THStack *simStack = new THStack("simStack", "A toy experiment");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();

    gPad->Print("TestCorrrection-nominal.png");
    
    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kSignalWeight] = -1.0;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Signal -1.0");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-signal1.png");

    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kSignalWeight] = +1.0;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Signal +1.0");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-signal2.png");

    ///////////////////////////////////////////////////////////
    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kBackgroundWeight] = -1.0;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Background -1.0");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-background1.png");

    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kBackgroundWeight] = +1.0;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Background +1.0");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-background2.png");

    ///////////////////////////////////////////////////////////
    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kMassScale] = -1.0;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Scale -1.0");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-scale1.png");

    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kMassScale] = +1.0;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Scale +1.0");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-scale2.png");

    ///////////////////////////////////////////////////////////
    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kMassWidth] = -1.0;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Width -1.0");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-width1.png");

    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kMassWidth] = +1.0;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Width +1.0");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-width2.png");

    ///////////////////////////////////////////////////////////
    double skew = 10.0;
    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kMassSkew] = -skew;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Negative Skew");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-skew1.png");

    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kMassSkew] = +skew;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Positive Skew");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-skew2.png");

    ///////////////////////////////////////////////////////////
    double sigsep = 5.0;
    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kSignalSeparationScale] = -sigsep;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Negative Sig. Sep.");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-sigsep1.png");

    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kSignalSeparationScale] = +sigsep;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Positive Sig. Sep.");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-sigsep2.png");

    ///////////////////////////////////////////////////////////
    double bkgsep = 5.0;
    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kBackgroundSeparationScale] = -bkgsep;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Negative Bkg. Sep.");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-bkgsep1.png");

    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kBackgroundSeparationScale] = +bkgsep;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Positive Bkg. Sep.");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-bkgsep2.png");

    ///////////////////////////////////////////////////////////
    double fake = 10.0;
    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kFakeMuDkProb] = -fake;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Negative Fake");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-fake1.png");

    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kFakeMuDkProb] = +fake;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Positive Fake");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-fake2.png");

    ///////////////////////////////////////////////////////////
    double mudk = 10.0;
    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kMuDkEfficiency] = -mudk;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Negative MuDk");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-mudk1.png");

    std::fill(params.begin(), params.end(), 0.0);
    params[SystematicCorrection::kMuDkEfficiency] = +mudk;
    likelihood.FillHistograms(params);
    simStack = new THStack("simStack", "Positive MuDk");
    simStack->Add(likelihood.SimulatedDecayTag);
    simStack->Add(likelihood.SimulatedSeparated);
    simStack->Add(likelihood.SimulatedClose);
    simStack->Draw();
    gPad->Print("TestLikelihood-mudk2.png");

    
}


