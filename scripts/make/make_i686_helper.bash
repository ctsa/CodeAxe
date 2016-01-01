#!/usr/bin/env bash
#
# helper for making a 32-bit p4 build on my x86-64 box
#

. ~/.bash_profile force_i686
make -j 2 gcc-linux-i686

