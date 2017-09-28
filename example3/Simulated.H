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

    void MakeSample(SampleType& sample, int signal, int background) {
        // Always simulate at least 1000 signal.
        if (signal<1000) signal = 1000;
        // Always simulate some background.
        if (background < signal) background = signal;
        std::cout << "Simulate"
                  << " " << signal << " signal events"
                  << " and " << background << " background events"
                  << std::endl;
        sample.reserve(signal+background);
        sample.clear();
        Event row;
        for (int i=0.0; i<signal; ++i) {
            MakeSignalEvent(row);
            sample.push_back(row);
        }
        while (background > 0) {
            MakeBackgroundEvent(row);
            sample.push_back(row);
            if (row.Mass < 500.0) --background;
        }
    }
    
    void MakeSignalEvent(Event& row) {
        row.TrueMass = 135.0;
        row.TrueMassSigma = 0.3*row.TrueMass;
        do {
            row.Mass = gRandom->Gaus(row.TrueMass,row.TrueMassSigma);
        } while (row.Mass < 0);
        row.Type = 0; // Signal
        row.Separation = std::abs(gRandom->Exp(150.0));
        if (gRandom->Uniform() < 0.05) row.MuDk = 1;
        else row.MuDk = 0;
    }

    void MakeBackgroundEvent(Event& row) {
        row.TrueMass = gRandom->Uniform(0.0,1000.0);
        row.TrueMassSigma = 0.4*row.TrueMass;
        do {
            row.Mass = gRandom->Gaus(row.TrueMass,row.TrueMassSigma);
        } while (row.Mass < 0);
        row.Type = 1; // Background
        row.Separation = std::abs(gRandom->Gaus(0.0,70.0));
        if (gRandom->Uniform() < 0.5) row.MuDk = 1;
        else row.MuDk = 0;
    }
    
};
#endif
