#ifndef FakeLikelihood_H_seen
#define FakeLikelihood_H_seen

#include "../TSimpleMCMC.H"

#include "FakeData.H"
#include "Simulated.H"
#include "SystematicCorrection.H"

#include "TH1D.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

// A likelihood similar to what might be used for the pizero analysis.
class FakeLikelihood {
public:
    /// The source for the toy data.
    FakeData ToyData;

    /// The "data" histograms.  These histograms are filled by ToyData and
    /// should never be changed.
    TH1* DataClose;
    TH1* DataSeparated;
    TH1* DataDecayTag;
    
    /// The sample of simulated events.  
    Simulated::SampleType SimulatedSample;

    /// The simulated histograms.  These are refilled during each trial, so
    /// they reflect the most recently tried point.
    TH1* SimulatedClose;
    TH1* SimulatedSeparated;
    TH1* SimulatedDecayTag;
    
    /// The corrections that are applied to each event.
    SystematicCorrection Corrections;

    /// Determine the number of dimensions.  This is where the dimensions are
    /// defined, and everything else uses it.
    std::size_t GetDim() const {return SystematicCorrection::kParamSize;}

    /// Calculate the likelihood.  This does a bin by bin comparision of the
    /// Data and Simulated histograms.
    double operator()(const Vector& point)  {
        ResetHistograms();
        FillHistograms(point);
        
        double logLikelihood = 0.0;

        for (int i=1; i<=DataClose->GetNbinsX(); ++i) {
            double data = DataClose->GetBinContent(i);
            double mc = SimulatedClose->GetBinContent(i);
            if (mc < 0.001) mc = 0.001;
            double v = data - mc;
            if (data > 0.0) v += data*std::log(mc/data);
            logLikelihood += v;
        }

        for (int i=1; i<=DataSeparated->GetNbinsX(); ++i) {
            double data = DataSeparated->GetBinContent(i);
            double mc = SimulatedSeparated->GetBinContent(i);
            if (mc < 0.001) mc = 0.001;
            double v = data - mc;
            if (data > 0.0) v += data*std::log(mc/data);
            logLikelihood += v;
        }

        for (int i=1; i<=DataDecayTag->GetNbinsX(); ++i) {
            double data = DataDecayTag->GetBinContent(i);
            double mc = SimulatedDecayTag->GetBinContent(i);
            if (mc < 0.001) mc = 0.001;
            double v = data - mc;
            if (data > 0.0) v += data*std::log(mc/data);
            logLikelihood += v;
        }

        return logLikelihood;
    }

    /// Initialize the likelihood.  Normally, this would read the data and
    /// simulated samples.  Instead, this randomly generates new toy data and
    /// a new simulated simpple.
    void Init(int dataSignal = 1000, int dataBackground = 1000,
              double mcOversample = 10.0) {
        std::cout << "Start initialization" << std::endl;

        // Make the toy data.  This creates the toy histograms.
        ToyData.FillData(dataSignal,dataBackground);
        DataClose = ToyData.Close;
        DataSeparated = ToyData.Separated;
        DataDecayTag = ToyData.DecayTag;
        
        // Make the simulated data.
        Simulated sim;
        sim.MakeSample(SimulatedSample,
                       mcOversample*dataSignal,
                       2*mcOversample*dataBackground);

        // Create the histograms to be filled using the simulation.
        SimulatedSeparated = (TH1D*)DataSeparated->Clone("simSeparated");
        SimulatedClose = (TH1D*)DataClose->Clone("simClose");
        SimulatedDecayTag = (TH1D*)DataDecayTag->Clone("simDecayTag");
        
        Corrections.ExposureRatio = 1.0;

        ResetHistograms();
        std::vector<double> point(GetDim());
        FillHistograms(point);

        double dataClose = DataClose->Integral();
        double dataSeparated = DataSeparated->Integral();
        double dataTagged = DataDecayTag->Integral();
        double data = dataClose + dataSeparated + dataTagged;

        double mcClose = SimulatedClose->Integral();
        double mcSeparated = SimulatedSeparated->Integral();
        double mcTagged = SimulatedDecayTag->Integral();
        double mc = mcClose + mcSeparated + mcTagged;

        std::cout << "Raw Integrals"
                  << " " << data
                  << " " << mc
                  << std::endl;
        std::cout << "  Data breakdown"
                  << " " << dataClose
                  << " " << dataSeparated
                  << " " << dataTagged
                  << std::endl;
        std::cout << "  MC breakdown"
                  << " " << mcClose
                  << " " << mcSeparated
                  << " " << mcTagged
                  << std::endl;

        Corrections.ExposureRatio = data/mc;
        std::cout << "Exposure Ratio " << Corrections.ExposureRatio
                  << std::endl;
        
        ResetHistograms();
        FillHistograms(point);

        dataClose = DataClose->Integral();
        dataSeparated = DataSeparated->Integral();
        dataTagged = DataDecayTag->Integral();
        data = dataClose + dataSeparated + dataTagged;

        mcClose = SimulatedClose->Integral();
        mcSeparated = SimulatedSeparated->Integral();
        mcTagged = SimulatedDecayTag->Integral();
        mc = mcClose + mcSeparated + mcTagged;


        std::cout << "Corrected Integrals"
                  << " " << data
                  << " " << mc
                  << std::endl;
        std::cout << "  Data breakdown"
                  << " " << dataClose
                  << " " << dataSeparated
                  << " " << dataTagged
                  << std::endl;
        std::cout << "  MC breakdown"
                  << " " << mcClose
                  << " " << mcSeparated
                  << " " << mcTagged
                  << std::endl;

    }

    void WriteSimulation(const Vector& point, std::string name) {
        ResetHistograms();
        FillHistograms(point);
        SimulatedClose->Clone((name+"Close").c_str())->Write();
        SimulatedSeparated->Clone((name+"Separated").c_str())->Write();
        SimulatedDecayTag->Clone((name+"DecayTag").c_str())->Write();
    }
    
    void ResetHistograms() {
        SimulatedClose->Reset();
        SimulatedSeparated->Reset();
        SimulatedDecayTag->Reset();
    }

    // This fills the simulated histogram with corrections applied.
    void FillHistograms(const std::vector<double>& params) {
#ifdef DUMP_FILL_PARAMETERS
        std::cout << "Fill with ";
        for (std::size_t i=0; i<params.size(); ++i) {
            std::cout << " " << params[i];
        }
        std::cout << std::endl;
#endif
        ResetHistograms();
        Simulated::Event corrected;
        for (std::size_t i = 0; i< SimulatedSample.size(); ++i) {
            double weight = Corrections.CorrectEvent(corrected,
                                                     SimulatedSample[i],
                                                     params);
            // Apply the cuts to see if the event passes.
            if (corrected.Mass > 500.0) continue;
            if (corrected.Mass < 0.0) continue;
            if (corrected.Separation < 0.0) continue;
            if (corrected.MuDk > 0) {
                SimulatedDecayTag->Fill(corrected.Mass,weight);
                continue;
            }
            if (corrected.Separation < 100.0) {
                SimulatedClose->Fill(corrected.Mass,weight);
                continue;
            }
            SimulatedSeparated->Fill(corrected.Mass,weight);
        }
    }
    
};
#endif
