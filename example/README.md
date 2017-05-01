# Example Data/MC likelihood using a binned goodness.

This is a simple test of the idea for using binned expectations for the MC
and binned histograms for the data.  It's modeled off of what would need to
be done for the P0D pizero analysis.  This particular example determines
the data/mc ratio, not the number of data events.

This is implemented using the TSimpleMCMC code in the parent directory.  It
produces a file with a tree of accepted steps.

It can be compile using the compile-fake.sh script and run as ./a.out, or
it can be run directly in ROOT using

root -l -q FakeMCMC.C++ 

I'm assuming that you have ROOT in your path!
