#ifndef SystematicCorrection_H_seen
#define SystematicCorrection_H_seen

#include "../TFakeGP.H"
#include "Simulated.H"
#include <iostream>

struct SystematicCorrection {
    typedef  std::pair<double, double> Correction;

    enum {
        kSignalWeight = 0,
        kBackgroundWeight,
        kSignalSeparationScale,
        kBackgroundSeparationScale,
        kFakeMuDkProb,
        kMuDkEfficiency,
        kMassScale,
        kMassWidth,
        kMassSkew,
        kBackgroundShapeBeg,
        kBackgroundShapeEnd = kBackgroundShapeBeg+10,
        kSignalShapeBeg,
        kSignalShapeEnd = kSignalShapeBeg+10,
        kParamSize
    };
    
    bool IsData(const Simulated::Event& evt) const {return (evt.Type < 0);}
    
    bool IsSignal(const Simulated::Event& evt) const {return (evt.Type == 0);}
    
    bool IsBackground(const Simulated::Event& evt) const {return (evt.Type>0);}
    
    double Separation(const Simulated::Event& evt) const {
        if (IsData(evt)) return evt.Separation;

        double scale = 0.0;
        if (IsSignal(evt)) scale += fParams[kSignalSeparationScale];
        if (IsBackground(evt)) scale += fParams[kBackgroundSeparationScale];

        scale = std::exp(scale/10.0);
        
        double sep = evt.Separation;

        return sep*scale;
    }
        
    double InvariantMass(const Simulated::Event& evt) const {
        double mass = evt.Mass;
        if (IsData(evt)) return mass;
        double nominalMass = evt.TrueMass;
        double nominalSigma = evt.TrueMassSigma;

        double nominalLogMass = std::log(nominalMass);
        double nominalLogSigma = std::log(nominalMass+nominalSigma);
        nominalLogSigma = nominalLogSigma - nominalLogMass;

        double logMass = std::log(mass);
        double logSigma = (logMass-nominalLogMass)/nominalLogSigma;

        double scale = fParams[kMassScale]/10.0;
        double width = std::exp(fParams[kMassWidth]/10.0);
        // Limit the skew to a valid range.  The skew function isn't defined
        // for skew values greater than +/- 0.3.
        double skew = 0.3*std::erf(fParams[kMassSkew]/10.0);
        skew = std::exp(logSigma*skew);

        // The order of the corrections matter.
        logMass = nominalLogMass + (logMass-nominalLogMass)*skew;
        logMass = nominalLogMass + (logMass-nominalLogMass)*width;
        logMass = logMass + scale;

        mass = std::exp(logMass);
        
        return mass;
    }
    
    double EventWeight(const Simulated::Event& evt) const {
        double weight = 1.0;
        if (IsData(evt)) return weight;
        
        // This does not apply a weight for the signal and background.  It's
        // assumed that the signal and background strength will be handled
        // external to the correction.  I've left this here so it's very clear
        // that the weighting is not happening.
#ifdef WEIGHT_SIGNAL_ANYWAY
        if (IsSignal(evt)) weight *= std::exp(fParams[kSignalWeight]/10.0);
        else weight *= std::exp(fParams[kBackgroundWeight]/10.0);
#endif
        // Apply weight for muon decay fake probability.  Only apply this to
        // the signal since the fake rate for the background can be covered by
        // the efficiency.
        double trueFakes = 0.05; // From Simulated.H
        double correctedFakes = std::tan(M_PI*(trueFakes-0.5));
        correctedFakes += fParams[kFakeMuDkProb]/10.0;
        correctedFakes = std::atan(correctedFakes)/M_PI + 0.5;
        
        if (IsSignal(evt)) {
            if (evt.MuDk>0) weight *= correctedFakes/trueFakes;
            else weight *= (1.0-correctedFakes)/(1.0-trueFakes);
            // This is reweighting against the uncorrected reconstructed mass,
            // not the corrected mass since this shape variation is
            // independent of the skew, width and energy scale.
            weight *= std::exp(SignalShape->GetValue(evt.Mass));
        }
        
        // Apply weight for muon decay efficiency.  Only apply this to
        // the background since the signal doesn't have any true muon decays.
        double trueEfficiency = 0.5; // From Simulated.H
        double correctedEfficiency = std::tan(M_PI*(trueEfficiency-0.5));
        correctedEfficiency += fParams[kMuDkEfficiency]/10.0;
        correctedEfficiency = std::atan(correctedEfficiency)/M_PI + 0.5;

        if (IsBackground(evt)) {
            if (evt.MuDk>0) weight *= correctedEfficiency/trueEfficiency;
            else weight *= (1.0-correctedEfficiency)/(1.0-trueEfficiency);
            // This is reweighting against the uncorrected reconstructed mass,
            // not the corrected mass since this shape variation is
            // independent of the skew, width and energy scale.
            weight *= std::exp(BackgroundShape->GetValue(evt.Mass));
        }
        
        return weight;
    }

    // Return the corrected event and the event weight.
    double CorrectEvent(Simulated::Event& corrected,
                        const Simulated::Event& evt) const {
        corrected.Mass = InvariantMass(evt);
        corrected.Type = evt.Type;
        corrected.Separation = Separation(evt);
        corrected.MuDk = evt.MuDk;
        return EventWeight(evt);
    }
    
    void SetParameters(const std::vector<double>& param) {
        if (fParams.size() != param.size()) fParams.resize(param.size());
        std::copy(param.begin(),param.end(),fParams.begin());
        for (int i=kBackgroundShapeBeg; i<=kBackgroundShapeEnd; ++i) {
            BackgroundShape->SetBinValue(i-kBackgroundShapeBeg,fParams[i]/10.0);
            // std::cout << "bkg " << i
            //           << " " << i-kBackgroundShapeBeg
            //           << " " << BackgroundShape->GetBinValue(i-kBackgroundShapeBeg) << std::endl;
        }
        for (int i=kSignalShapeBeg; i<=kSignalShapeEnd; ++i) {
            SignalShape->SetBinValue(i-kSignalShapeBeg+1,fParams[i]/10.0);
        }

    }
    
    SystematicCorrection() {
        BackgroundShape
            = new TFakeGP("backgroundShape", 0.0, 500.0,
                          kBackgroundShapeEnd - kBackgroundShapeBeg + 1);
        BackgroundShape->GaussianKernel(100.0,0.7);
        // The two additional bins are because the end points are fixed to
        // zero and not allowed to vary so there aren't free parameters for
        // those bins.
        SignalShape
            = new TFakeGP("signalShape", 0.0, 250.0,
                          kSignalShapeEnd - kSignalShapeBeg + 1 + 2);
        SignalShape->GaussianKernel(50.0);
    }

    TFakeGP* BackgroundShape;
    TFakeGP* SignalShape;

    std::vector<double> fParams;
};
#endif
