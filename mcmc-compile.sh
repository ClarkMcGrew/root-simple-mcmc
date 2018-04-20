#!/bin/bash

$(root-config --cxx) $(root-config --cflags) \
		     -o mcmc.exe \
		     -DMAIN_PROGRAM SimpleMCMC.C \
		     $(root-config --libs)
