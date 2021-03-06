#ifndef FakeData_H_seen
#define FakeData_H_seen

#include <TH1D.h>
#include <TRandom.h>
#include <cmath>
#include <iostream>

// Create a sample of fake data.  This explicitly uses a different toy model
// for the fake data than is used in the simulated sample.  That let's me play
// with fitting the wrong model.
struct FakeData {
    struct Event {
        double Mass;
        int Type;
        double Separation;
        int MuDk;
        double TrueMass;
        double TrueMassSigma;
    };
    
    typedef std::vector< Event > SampleType;

    SampleType DataSample;

    TH1* VeryClose;

    TH1* Close;

    TH1* Separated;
    
    TH1* DecayTag;

    void MakeSample(int signal, int bkg) {
        DataSample.reserve(signal+bkg);
        DataSample.clear();
        TrueSignal = signal;
        TrueBackground = bkg;
        do{TrueScale = gRandom->Gaus(1.0,0.01);}while (TrueScale < 0.80);
        do{TrueResolution=gRandom->Gaus(0.4,0.05);}while (TrueResolution<0.15);
        TrueSignalSeparation = gRandom->Gaus(150.0,10.0);
        TrueBackgroundSeparation = gRandom->Gaus(70.0,10.0);
        do{TrueFakeMuDk = gRandom->Gaus(0.05,0.01);}while(TrueFakeMuDk<0.03);
        do{TrueMuDkFrac = gRandom->Gaus(0.50,0.01);}while(TrueMuDkFrac<0.3);

        SampleMinimum = 0.0;
        SampleMaximum = 500.0;

        std::cout << "True Signal " << TrueSignal << std::endl;
        std::cout << "True Background " << TrueBackground << std::endl;
        std::cout << "True Mass Scale " << TrueScale << std::endl;
        std::cout << "True Mass Resolution " << TrueResolution << std::endl;
        std::cout << "True Signal Separaration "
                  << TrueSignalSeparation
                  << std::endl;
        std::cout << "True Background Separaration "
                  << TrueBackgroundSeparation
                  << std::endl;
        std::cout  << "True Muon Decay Fake Probability " << TrueFakeMuDk
                  << std::endl;
        std::cout  << "True Muon Decay Efficiency " << TrueMuDkFrac
                  << std::endl;
        
        Event row;
        for (int i=0; i<TrueSignal; ++i) {
            MakeSignalEvent(row);
            DataSample.push_back(row);
        }
        for (int i=0; i<TrueBackground; ++i) {
            MakeBackgroundEvent(row);
            DataSample.push_back(row);
        }
    }

    void MakeSignalEvent(Event& row) {
#define SKIP_LOG_GAUSS_RESOLUTION
#ifndef SKIP_LOG_GAUSS_RESOLUTION
        double resWidth = std::log(1.0+TrueResolution);
        double v = gRandom->Gaus(0.0,resWidth);
        double m = TrueScale*135.0*std::exp(v);
        while (m<SampleMinimum || m>SampleMaximum) {
            v = gRandom->Gaus(0.0,resWidth);
            m = TrueScale*135.0*std::exp(v);
        }
        row.Mass = m;
#else
        TrueResolution = 0.3;
        double resWidth = TrueScale*TrueResolution*135.0;
        do {row.Mass = gRandom->Gaus(TrueScale*135.0,resWidth);}
        while (row.Mass<SampleMinimum || row.Mass>SampleMaximum);
#endif
        row.TrueMass = 135.0;
        row.TrueMassSigma = TrueResolution*135.0;
        row.Type = -1; // Signal (don't look at this for data!);
        row.Separation = std::abs(gRandom->Exp(TrueSignalSeparation));
        if (gRandom->Uniform() < TrueFakeMuDk) row.MuDk = 1;
        else row.MuDk = 0;
    }

    void MakeBackgroundEvent(Event& row) {
        row.Mass = gRandom->Uniform(SampleMinimum,SampleMaximum);
        row.TrueMassSigma = TrueResolution*row.TrueMass;
        do {
            // Yes, this looks backwards.  No, it's not.
            row.TrueMass = gRandom->Gaus(row.Mass,row.TrueMassSigma);
        } while (row.TrueMass < 0);
        row.Type = -2; // Background (don't look at this for data!);
        row.Separation = std::abs(gRandom->Gaus(0.0,TrueBackgroundSeparation));
        if (gRandom->Uniform() < TrueMuDkFrac) row.MuDk = 1;
        else row.MuDk = 0;
    }

    void FillData(int signal, int bkg) {
        MakeSample(signal,bkg);

        DecayTag = new TH1D("DataDecayTag", "A Toy Sample w/  Decay Tag",
                            50, SampleMinimum, SampleMaximum);
        DecayTag->SetLineColor(kBlack);
        
        Separated = new TH1D("DataSeparated", "A Toy Sample w/ Separation",
                             50, SampleMinimum, SampleMaximum);
        Separated->SetLineColor(kBlue);

        Close = new TH1D("DataClose", "A Toy Sample w/ Small Separation",
                             50, SampleMinimum, SampleMaximum);
        Close->SetLineColor(kRed);

        VeryClose = new TH1D("DataVeryClose", "A Toy Sample w/o Separation",
                             50, SampleMinimum, SampleMaximum);
        VeryClose->SetLineColor(kGreen);

        for (std::size_t i = 0; i< DataSample.size(); ++i) {
            if (DataSample[i].MuDk>0) {
                DecayTag->Fill(DataSample[i].Mass);
                continue;
            }
            if (DataSample[i].Separation < 50.0) {
                VeryClose->Fill(DataSample[i].Mass);
                continue;
            }
            if (DataSample[i].Separation < 100.0) {
                Close->Fill(DataSample[i].Mass);
                continue;
            }
            Separated->Fill(DataSample[i].Mass);
        }
    }

    // The true scale.  A value of 1.0 means the reconstructed pizero mass is
    // 135 MeV/c^2
    double TrueScale;

    // The true resolution (in 0.2 is twenty percent).
    double TrueResolution;

    // The average separation for signal.
    double TrueSignalSeparation;

    // The average separation for background
    double TrueBackgroundSeparation;

    // The minimum bound of the sample (usually 0)
    double SampleMinimum;

    // The maximum bound of the sample (usually 500)
    double SampleMaximum;

    // The number of true signal events.
    int TrueSignal;

    // The number of true background events.
    int TrueBackground;

    // The true fake rate for muon decays.
    double TrueFakeMuDk;

    // The true fraction of muon decays in background.
    double TrueMuDkFrac;
    
};
#endif
