#ifndef Simulated_h_seen
#define Simulated_h_seen

#include <TRandom.h>

struct Simulated {
    struct Event {
        double Mass;
        int Type;
        double Separation;
        int MuDk;
        double TrueMass;
        double TrueMassSigma;
    };
    typedef std::vector< Event > SampleType;

    void MakeSample(SampleType& sample,int signal, int background) {
        sample.reserve(signal+background);
        sample.clear();
        Event row;
        for (int i=0.0; i<signal; ++i) {
            MakeSignalEvent(row);
            sample.push_back(row);
        }
        for (int i=0.0; i<background; ++i) {
            MakeBackgroundEvent(row);
            sample.push_back(row);
        }
    }
    
    void MakeSignalEvent(Event& row) {
        row.TrueMass = 135.0;
        row.TrueMassSigma = 0.3*row.TrueMass;
        do {
            row.Mass = gRandom->Gaus(row.TrueMass,row.TrueMassSigma);
        } while (row.Mass < 0);
        row.Type = 0; // Signal
        row.Separation = std::abs(gRandom->Exp(100.0));
        if (gRandom->Uniform() < 0.05) row.MuDk = 1;
        else row.MuDk = 0;
    }

    void MakeBackgroundEvent(Event& row) {
        row.TrueMass = gRandom->Uniform(0.0,1000.0);
        row.TrueMassSigma = 0.3*row.TrueMass;
        do {
            row.Mass = gRandom->Gaus(row.TrueMass,row.TrueMassSigma);
        } while (row.Mass < 0);
        row.Type = 1; // Background
        row.Separation = std::abs(gRandom->Gaus(0.0,50.0));
        if (gRandom->Uniform() < 0.5) row.MuDk = 1;
        else row.MuDk = 0;
    }
    
};
#endif
