#!/bin/bash

$(root-config --cxx) $(root-config --cflags) -g \
		     -o vaat.exe \
		     -DMAIN_PROGRAM SimpleVAAT.C \
		     $(root-config --libs)
