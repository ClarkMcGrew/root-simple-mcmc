#!/bin/bash

$(root-config --cxx) $(root-config --cflags) \
		     -DMAIN_PROGRAM SimpleMCMC.C \
		     $(root-config --libs)
