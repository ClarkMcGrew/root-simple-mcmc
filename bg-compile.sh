#!/bin/bash

$(root-config --cxx) $(root-config --cflags) \
		     -o bg.exe \
		     -DMAIN_PROGRAM BadGrad.C \
		     $(root-config --libs)
