#!/bin/bash

$(root-config --cxx) $(root-config --cflags) \
		     -DMAIN_PROGRAM SimpleHMC.C \
		     $(root-config --libs)
