#!/bin/bash

$(root-config --cxx) $(root-config --cflags) \
		     -o hmc.exe \
		     -DMAIN_PROGRAM SimpleHMC.C \
		     $(root-config --libs)
