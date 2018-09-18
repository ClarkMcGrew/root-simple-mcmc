#!/bin/bash

$(root-config --cxx) $(root-config --cflags) \
		     -o constrained.exe \
		     -DMAIN_PROGRAM Constrained.C \
		     $(root-config --libs)
