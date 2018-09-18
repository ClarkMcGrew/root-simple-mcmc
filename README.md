# root-simple-mcmc

A simple MCMC template for use with ROOT (tested with 5.34+ and 6.06+).  It
can be used as in a macro with ACLiC, or directly in C++ code.

The TSimpleMCMC templated class runs an Markov Chain Monte Carlo using a
user provided likelihood and stepping proposal.  The resulting MCMC
normally uses default stepping proposal which implements an adaptive
proposal function.  

The MCMC object is created with one required and one optional template argument.

```
typedef TSimpleMCMC<UserLikelihood> TUserMCMC;
```

or 
	
```
typedef TSimpleMCMC<UserLikelihood,UserProposal> TUserMCMC;
```

The UserLikelihood template argument must be a class (or struct) which
provides a method declared as:

```
struct ExampleLogLikelihood {
   double operator() (const std::vector<double>& point);
}
```

The optional UserProposal template argument must be a struct (or class)
which provides a method declared as:

```
struct ExampleProposal {
   void operator() (std::vector<double>& proposed,
                      const std::vector<double>& previous,
                      const double previousValue);
}
```

where proposed is the new proposed point, previous is the last
accepted point, and previousValue is the log Likelihood at the
previous point.

This can be used in your root macros, or C++ code as follows.

```
class TDummyLogLikelihood {
public:
    double operator()(const Vector& point)  const {
        double logLikelihood = 0.0;
        for (std::size_t i = 0; i < point.size(); ++i) {
            logLikelihood += - 0.5*point[i]*point[i];
        }
        return logLikelihood;
    }
};

void SimpleMCMC() {

    TFile *outputFile = new TFile("simple-mcmc.root","recreate");
    TTree *tree = new TTree("SimpleMCMC",
                            "Tree of accepted points");

    TSimpleMCMC<TDummyLogLikelihood> mcmc(tree);

    // The next three lines are for example only and should almost never
    // be used.  They are documented in the TProposeAdaptiveStep header file.  
    // This is to show the syntax for controlling TProposeAdaptiveStep, don't 
    // just copy this blindly!
    mcmc.GetProposeStep().SetDim(5);          // Not needed!
    mcmc.GetProposeStep().SetGaussian(3,2.0); // Not recommended!
    mcmc.GetProposeStep().SetUniform(4,-5,5); // Maybe for a special case.

    Vector point(5);
    mcmc.Start(point);  // Set initial conditions

    for (int i=0; i<1000000; ++i) mcmc.Step(false);  // Burn-in the chain.
    for (int i=0; i<1000000; ++i) mcmc.Step();       // Run the chain.

    tree->Write();
    delete outputFile;
}
```

If this macro is in a file named "SimpleMCMC.C", then you can run it using

```
root -l -q SimpleMCMC.C+
```

The macro must be compiled using ACLIC.

The default class for UserProposal is TProposeAdaptiveStep which
implements an adaptive Metropolis-Hastings step.  It has several methods
that can be accessed using the GetProposeStep() method.  See above for an
example.

# MCMC Versions Here

This repo actually contains a few different MCMC examples that I have used
to learn about the different types of behaviors.  None of these examples is
intended as an end user program, and I control a lot of the input
parameters by editing the source and recompiling (hey, it's test code).
However, the associated TSimple<blah>.H classes are fairly well tested.

- TSimpleMCMC.H (and friends) : This is the adaptive MCMC described above.
It's the best tested class, and is my first choice when I'm looking at the
behavior of an MCMC.  An alternative for the proposal is provided by
TProposeGibbsStep.h (the Gibbs step is not adaptive).

- TSimpleHMC.H (and friends) : This is a "pure" Hamiltonian MC.  It handles
the relatively rare special case where you can write down the derivative of
the likelihood, but for the right problem it converges much more quickly.

TSimpleAHMC.H (and friends) -- This is an HMC implementation that uses an
approximate version of the gradient.  The gradient is estimated based on
the accumulated covariance of the posterior.  I avoid this because it's not
faster than TSimpleMCMC, and doesn't seem to be as reliable.  My feeling is
that it makes to many approximations.

BadGrad.C -- This is just a toy to see how accurately the gradient needs to
be calculated.

# Other Tools

A (tiny) collection of other tools to look at the results of the MCMC
chains have been included.  In general, they are in the form of ROOT macros
(i.e. .C files), and have running instructions in the comments at the top
of the file.

MakeCovariance.C -- Read a root tree containing an MCMC chain (for example,
one produced by SimpleMCMC.C), and produce a covariance matrix for the
posterior.  The results are saved in histograms.

CholeskyChain.C -- Get the mean and covariance (as produced by
MakeCovariance.C) from a pair of histograms, and then produce a "chain"
using Cholesky Decomposition.   
