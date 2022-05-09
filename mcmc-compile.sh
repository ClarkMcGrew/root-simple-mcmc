#!/bin/bash

$(root-config --cxx) $(root-config --cflags) -O2 \
		     -o mcmc.exe \
		     -DMAIN_PROGRAM SimpleMCMC.C \
		     $(root-config --libs)
