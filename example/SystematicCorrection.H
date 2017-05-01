#ifndef SystematicCorrection_H_seen
#define SystematicCorrection_H_seen

#include "Simulated.H"

#include <iostream>

struct SystematicCorrection {
    typedef  std::pair<double, double> Correction;

    enum {
        kSignalWeight = 0,
        kBackgroundWeight,
        kMassScale,
        kMassWidth,
        kMassSkew,
        kSignalSeparationScale,
        kBackgroundSeparationScale,
        kFakeMuDkProb,
        kMuDkEfficiency,
        kParamSize,
    };
    
    // This should be the ratio between the raw data sample and the raw
    // simulated sample.  It keeps the signal weight and background weight
    // close to unity.
    double ExposureRatio;

    bool IsData(const Simulated::Event& evt) const {return (evt.Type < 0);}
    
    bool IsSignal(const Simulated::Event& evt) const {return (evt.Type == 0);}
    
    bool IsBackground(const Simulated::Event& evt) const {return (evt.Type > 0);}
    
    double Separation(const Simulated::Event& evt,
                      const std::vector<double>& param) const {
        if (IsData(evt)) return evt.Separation;

        double scale = 0.0;
        if (IsSignal(evt)) scale += param[kSignalSeparationScale];
        if (IsBackground(evt)) scale += param[kBackgroundSeparationScale];

        scale = std::exp(scale/10.0);
        
        double sep = evt.Separation;

        return sep*scale;
    }
        
    double InvariantMass(const Simulated::Event& evt,
                         const std::vector<double>& param) const {
        double mass = evt.Mass;
        if (IsData(evt)) return mass;
        double nominalMass = evt.TrueMass;
        double nominalSigma = evt.TrueMassSigma;

        double nominalLogMass = std::log(nominalMass);
        double nominalLogSigma = std::log(nominalMass+nominalSigma);
        nominalLogSigma = nominalLogSigma - nominalLogMass;

        double logMass = std::log(mass);
        double logSigma = (logMass-nominalLogMass)/nominalLogSigma;

        double scale = param[kMassScale]/10.0;
        double width = std::exp(param[kMassWidth]/10.0);
        // Limit the skew to a valid range.  The skew function isn't defined
        // for skew values greater than +/- 0.3.
        double skew = 0.3*std::erf(param[kMassSkew]/10.0);
        skew = std::exp(logSigma*skew);

        // The order of the corrections matter.
        logMass = nominalLogMass + (logMass-nominalLogMass)*skew;
        logMass = nominalLogMass + (logMass-nominalLogMass)*width;
        logMass = logMass + scale;

        mass = std::exp(logMass);
        
        return mass;
    }
    
    double EventWeight(const Simulated::Event& evt,
                       const std::vector<double>& param) const {
        double weight = 1.0;
        if (IsData(evt)) return weight;
        
        // Apply weight for signal (or background) normalization.
        if (IsSignal(evt)) weight *= std::exp(param[kSignalWeight]/10.0);
        else weight *= std::exp(param[kBackgroundWeight]/10.0);

        // Apply weight for muon decay fake probability.  Only apply this to
        // the signal since the fake rate for the background can be covered by
        // the efficiency.
        double trueFakes = 0.05; // From Simulated.H
        double correctedFakes = std::tan(M_PI*(trueFakes-0.5));
        correctedFakes += param[kFakeMuDkProb];
        correctedFakes = std::atan(correctedFakes)/M_PI + 0.5;
        
        if (IsSignal(evt)) {
            if (evt.MuDk>0) weight *= correctedFakes/trueFakes;
            else weight *= (1.0-correctedFakes)/(1.0-trueFakes);
        }
        
        // Apply weight for muon decay efficiency.  Only apply this to
        // the background since the signal doesn't have any true muon decays.
        double trueEfficiency = 0.5; // From Simulated.H
        double correctedEfficiency = std::tan(M_PI*(trueEfficiency-0.5));
        correctedEfficiency += param[kMuDkEfficiency];
        correctedEfficiency = std::atan(correctedEfficiency)/M_PI + 0.5;

        if (IsBackground(evt)) {
            if (evt.MuDk>0) weight *= correctedEfficiency/trueEfficiency;
            else weight *= (1.0-correctedEfficiency)/(1.0-trueEfficiency);
        }
        
        weight *= ExposureRatio;
        return weight;
    }

    // Return the corrected event and the event weight.
    double CorrectEvent(Simulated::Event& corrected,
                        const Simulated::Event& evt,
                        const std::vector<double>& param) const {
        corrected.Mass = InvariantMass(evt,param);
        corrected.Type = evt.Type;
        corrected.Separation = Separation(evt,param);
        corrected.MuDk = evt.MuDk;
        return EventWeight(evt,param);
    }
    
};
#endif
