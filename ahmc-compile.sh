#!/bin/bash

$(root-config --cxx) $(root-config --cflags) \
		     -o ahmc.exe \
		     -DMAIN_PROGRAM SimpleAHMC.C \
		     $(root-config --libs)
