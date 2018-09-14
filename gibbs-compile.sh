#!/bin/bash

$(root-config --cxx) $(root-config --cflags) \
		     -o gibbs.exe \
		     -DMAIN_PROGRAM SimpleGibbs.C \
		     $(root-config --libs)
