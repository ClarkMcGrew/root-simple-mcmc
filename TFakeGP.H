#ifndef FakeGP_H_SEEN
#define FakeGP_H_SEEN

#include <iostream>

#include <TH1.h>
#include <TMatrixD.h>

#ifndef FakeGP_DEBUG_LEVEL
#define FakeGP_DEBUG_LEVEL 2
#endif

#ifndef FakeGP_DEBUG
#define FakeGP_DEBUG(level) if (level <= (FakeGP_DEBUG_LEVEL)) std::cout
#endif

#define FakeGP_ERROR (std::cout <<__FILE__<<":: " << __LINE__ << ": " )

/// A Gaussian process provides a method for determining the possible
/// variations of a function between measured points, and is a very powerful
/// method for approximating functions.  For a real reference, check Rasmussen
/// & Williams, "Gaussian Processes for Machine Learning" (available on line
/// at www.GaussianProcess.org/gpml).  This class takes the general idea of a
/// Gaussian process (i.e. the being able to set the correlation between
/// points) and simplifies the overall behaviour to have linear interpolation
/// between the control (i.e. measured) points.  This is much faster than
/// dealing with the necessary matrix inversions when handling a real GP, and
/// is good enough for finding a set of interpolated weights.
///
/// The object is created between a lower and upper bound with a fixed number
/// of control points.  It's assumed that the results of this will be a
/// perturbation relative to some other function, so the mean value at each
/// control point is taken to be zero.  The kernel can then be set using the
/// SetKernel method, or if you want to use a Gaussian (or Exponential)
/// kernel, you can use the GaussianKernel (or ExponentialKernel) with a
/// coherence length.
class TFakeGP {
public:

    /// Create a fake Gaussian Process.  The first argument is a name used to
    /// construct an internal histogram, the low and high arguments specify
    /// the full range of the process, and the number of bins is set with the
    /// final argument.
    TFakeGP(const char* name, double low, double high, int bins) {
        fHist = new TH1D(name,name,bins,low,high);
        fKernel.ResizeTo(bins,bins);
    }
    
    ~TFakeGP() {
        delete fHist;
    }

    /// Get the total number of bins in the fake GP.  This is the value
    /// provided in the constructor.
    double GetBinCount() {
        return fHist->GetNbinsX();
    }
    
    /// Get the bin centers.
    double GetBinCenter(int i) {
        return fHist->GetBinCenter(i+1);
    }
    
    /// Set the control value for a bin.
    void SetBinValue(int i, double v) {
        fHist->SetBinContent(i+1,v);
    }

    /// Get the control value for a bin.
    double GetBinValue(int i) {
        return fHist->GetBinContent(i+1);
    }

    /// Get the value for a particular independent variable.  This uses linear
    /// interpolation between the control bins.  Handles the case where the
    /// independent variable is out of range by returning the value from the
    /// closest bin.
    double GetValue(double v) {
        return fHist->Interpolate(v);
    }

    /// Set the kernel between two bins.  This forces a symmetric matrix.
    void SetKernel(int i, int j, double v) {
        fKernel(i,j) = v;
        if (i != j) fKernel(j,i) = v;
    }

    /// Get the kernel value.
    double GetKernel(int i, int j) {
        return fKernel(i,j);
    }

    /// Construct a Gaussian kernel with a set coherence length.  The
    /// variation of the function around a mean value of zero can also be set.
    double GaussianKernel(double coherence, double sigma = 1.0) {
        for (int i=0; i<GetBinCount(); ++i) {
            for (int j=i; j<GetBinCount(); ++j) {
                double r = fHist->GetBinCenter(i+1) - fHist->GetBinCenter(j+1);
                r = r/coherence;
                r = 0.5*r*r;
                if (r > 20) r = 0.0;
                else r = sigma*sigma*std::exp(-r);
                SetKernel(i,j,r);
            }
        }
        fKernelDecomp.ResizeTo(0,0);
        fKernelInv.ResizeTo(0,0);
    }

    /// Construct an exponential kernel with a set coherence length.  The
    /// variation of the function around a mean value of zero can also be set.
    double ExponentialKernel(double coherence, double sigma = 1.0) {
        for (int i=0; i<GetBinCount(); ++i) {
            for (int j=i; j<GetBinCount(); ++j) {
                double r = fHist->GetBinCenter(i+1) - fHist->GetBinCenter(j+1);
                r = r/coherence;
                r = std::abs(r);
                if (r > 20) r = 0.0;
                else r = sigma*sigma*std::exp(-r);
                SetKernel(i,j,r);
            }
        }
        fKernelDecomp.ResizeTo(0,0);
        fKernelInv.ResizeTo(0,0);
    }

    /// Get the penalty for the current set of values.  This is the
    /// chi-squared of the control points relative to the expected value of
    /// zero for each control point.
    double GetPenalty() {
        if (fKernelInv.GetNcols() != fKernel.GetNcols()) {
            fKernelInv.ResizeTo(fKernel.GetNcols(),fKernel.GetNcols());
            fKernelInv = fKernel;
            fKernelInv.Invert();
        }
        double penalty = 0.0;
        for (int i=0; i<GetBinCount(); ++i) {
            for (int j=0; j<GetBinCount(); ++j) {
                double r = fHist->GetBinContent(i+1)*fHist->GetBinContent(j+1)
                    * fKernelInv(i,j);
                penalty += r;
            }
        }
        return penalty;
    }

    /// Generate a set of weights with the right correlations.  This hasn't
    /// been debugged.  It could be used to generate MCMC steps.
    void MakeProposal() {
        if (fKernelDecomp.GetNcols() != fKernel.GetNcols()) {
            fKernelDecomp.ResizeTo(fKernel.GetNcols(),fKernel.GetNcols());
            fKernelDecomp = fKernel;
            TDecompChol chol(fKernelDecomp);
            if (chol.Decompose()) {
                fKernelDecomp = chol.GetU();
            }
        }
        for (std::size_t i = 0; i < GetBinCount(); ++i) {
            fHist->SetBinContent(i+1,0.0);
        }
        for (std::size_t i = 0; i < GetBinCount(); ++i) {
            double r = gRandom->Gaus(0.0,1.0);
            for (std::size_t j = 0; j < GetBinCount(); ++j) {
                double v = fHist->GetBinContent(j+1);
                v += r*fKernelDecomp(i,j);
                fHist->SetBinContent(j+1,v);
            }
        }
    }

    /// Get the internal histogram
    TH1* GetHistogram() const {return fHist;}
    
public:
    TH1* fHist;
    TMatrixD fKernel;
    TMatrixD fKernelInv;
    TMatrixD fKernelDecomp;
};
    
#endif
