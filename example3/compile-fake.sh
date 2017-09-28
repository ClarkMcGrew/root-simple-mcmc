#!/bin/bash

$(root-config --cxx) $(root-config --cflags) \
		     -DMAIN_PROGRAM FakeMCMC.C \
		     $(root-config --libs) \
		     -o fake-mcmc.exe
