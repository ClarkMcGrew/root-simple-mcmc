# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

A header-only, ROOT-based MCMC template library (`TSimpleMCMC.H` and
friends), plus driver macros, chain-analysis macros, and four example
directories. There is no build system, no library to install, and no test
suite -- a header is "installed" by copying it into another project.

## Build and run

Everything compiles through `root-config`, so ROOT must be on `PATH` (6.36.04 here).

```bash
./mcmc-compile.sh                       # SimpleMCMC.C  -> mcmc.exe  (-O2 -Wall)
./mcmc.exe <cycles> <steps> <out.root> [restore.root]
```

Sibling scripts follow the same pattern: `hmc-compile.sh` -> `hmc.exe`,
`ahmc-compile.sh` -> `ahmc.exe`, `vaat-compile.sh` -> `vaat.exe`, `bg-compile.sh` ->
`bg.exe`, `example*/compile-fake.sh`, `example4/constrained-compile.sh`. Each just runs
`$(root-config --cxx) $(root-config --cflags) -DMAIN_PROGRAM <driver>.C $(root-config --libs)`.

The README's claim that `mcmc.exe` takes three arguments is stale -- `main()` at the
bottom of `SimpleMCMC.C` takes four: cycles, steps, output file, restore-from file.

Any `.C` driver can also be run as a ROOT macro (`root -l -q SimpleMCMC.C+`);
`./cleanTools.sh` removes the `*_C.d`, `*_C.so`, and `*_rdict.pcm` ACLiC leftovers.

## Testing

There are no unit tests. Verification means running a chain against one of the
`T*LogLikelihood.H` test functions and inspecting it with the analysis macros, each of
which documents its inputs and outputs in a header comment:

```bash
root chain.root MakeCovariance.C          # -> covariance.root (AcceptedCovariance, AcceptedMean)
root chain.root MakeAutocorrelation.C     # -> autocorrelation.root (autocorrelation vs lag)
root covariance.root CholeskyChain.C      # -> Cholesky.root (synthetic chain from the covariance)
root chain.root TestMarginalization.C     # for THardLogLikelihood, low dimension only
```

`continue-chain.sh` runs chains in sequence or in parallel on a cluster. It substitutes
`:INPUT:` and `:OUTPUT:` in the command and manages a
`prefix_epoch_parentmd5_childmd5.{open,closed,input,running}.root` naming convention to
hand files from one job to the next. Its header comment is the manual.

## Architecture

Everything is in namespace `sMCMC`, where `Parameter` is `double` and `Vector` is
`std::vector<Parameter>`.

`TSimpleMCMC<UserLikelihood, UserProposal = TProposeAdaptiveStep>` is the
production class. The lifecycle is `Start(point)`, a loop of `Step()`, then
`SaveStep()`.  `Step(save, metropolis)` takes a debugging mode: `0` is a
normal Metropolis step, `1` takes only uphill steps (a very slow
maximizer), `2` accepts everything (a likelihood scan).

The **likelihood contract** is only `double operator()(const
sMCMC::Vector&)`. The drivers additionally call `Init()` and `GetDim()` on
the likelihood, but those are driver conventions, not template
requirements.

The **proposal contract** is six methods, documented in full at the top of
`TSimpleMCMC.H`. `operator()` fills the proposed point and returns
`log(g(new|old)/g(old|new))` -- **0.0 for a symmetric proposal**, which is
the usual case; `TSimpleMCMC::Step` adds this to the log likelihood
difference. The other five -- `InitializeState`, `RestoreState`,
`AttachState`, `SaveState`, `StateSaved` -- may all be no-ops;
`TProposeSimpleStep` in the same header is the minimal implementation.

Proposals available:

- `TProposeAdaptiveStep` (default, in `TSimpleMCMC.H`) -- adaptive
  Metropolis-Hastings.  It accumulates the posterior covariance,
  Cholesky-decomposes it for the step direction, and tunes a scalar
  `fSigma` toward a target acceptance. Reach its knobs through
  `mcmc.GetProposeStep()`: `SetDim`, `SetGaussian`, `SetUniform`,
  `SetCorrelation`, `SetCovarianceWindow`, `SetAcceptanceWindow`,
  `SetAcceptanceRigidity`, `SetTargetAcceptance`, `SetSigma`,
  `UpdateProposal`, `ResetProposal`. Adapting breaks the strict Markov
  property, so update infrequently.
- `TProposeVAATStep.H` -- variable-at-a-time, also usable in production.
- `TSimpleHMC.H` -- a separate template, `TSimpleHMC<Likelihood, Gradient>`,
  for Hamiltonian MC. It needs a gradient (or will approximate one from the
  running covariance) and uses its own `HMC_DEBUG` / `HMC_ERROR` /
  `HMC_DEBUG_LEVEL` macros.

## Chain state and continuation

The output `TTree` gets `LogLikelihood`, `TotalSteps`, `Accepted`, and
`StepRMS` from the MCMC (plus `Step` when the constructor's `saveStep`
argument is true). The proposal then appends its own branches through
`AttachState()` -- `TProposeAdaptiveStep` writes `Adaptive*` branches
holding sigma, acceptance, the estimated central point, and the packed
lower triangle of the covariance.

This is what makes chains continuable: `Restore(tree)` reads those branches
back out of a previous run's tree. Two consequences worth remembering:

- A run **must** end with `mcmc.SaveStep()`. Without it the final proposal
  state is never written and the chain cannot be continued.
- `StateSaved()` zeroes the `Adaptive*` branch variables after each fill to
  save space, so nearly every entry holds zeros -- only the final forced
  save is meaningful.

## Conventions and gotchas

- `*.H` is a header-only class; `*.C` is a driver with an `#ifdef
  MAIN_PROGRAM main()` at the bottom, so the same file works both as a ROOT
  macro and as a standalone program.
- Drivers are configured by editing `#define`s and recompiling --
  `USE_HARD_LIKELIHOOD`, `USE_HORRIFIC_LIKELIHOOD`, `USE_ASYM_LIKELIHOOD`,
  `BURNIN_CHAIN`, `SKIP_MCMC`, `SCAN_MCMC`, `NO_OUTPUT`,
  `USE_THIS_PROPOSAL`. This is deliberate; they are test code.
- Debug output goes through `MCMC_DEBUG(level) << ...`, gated by
  `MCMC_DEBUG_LEVEL` (default 2), and `MCMC_ERROR` for errors.
- `example/`, `example2/`, and `example4/` have been updated for the move
  into the `sMCMC` namespace.  `example3/` has not: it still uses
  unqualified `TSimpleMCMC<>` and `Vector`, so it needs `sMCMC::`
  qualification before it will compile.  It is deliberately left alone
  because it was written to test the `TFakeGP.H` idea, which is not being
  used.
- `TFakeGP.H` does not compile on its own.  `MakeProposal()` uses
  `TDecompChol` and `gRandom` without including `<TDecompChol.h>` or
  `<TRandom.h>`, and `GaussianKernel()` and `ExponentialKernel()` are
  declared to return `double` but return nothing.  It builds inside
  `example3/` only because `TSimpleMCMC.H` is included first and supplies
  the missing declarations.
- Commit messages are short imperative one-liners, e.g. "Add explicit
  override for the 'sigma' step size".
- Limit the character set used on commits, source and documentation to
  UTF-7
