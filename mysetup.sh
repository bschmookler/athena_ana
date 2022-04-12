#!/bin/sh
export DETECTOR=athena
## May want to use source directory instead...
## ..if you use the 'share' version, you will need to make a softlink to ip6
export DETECTOR_PATH=${ATHENA_PREFIX}/share/athena
export DETECTOR_VERSION=master
export BEAMLINE_CONFIG=ip6
export BEAMLINE_CONFIG_VERSION=master
## note: we will phase out the JUGGLER_* flavor of variables in the future
export JUGGLER_DETECTOR=$DETECTOR
export JUGGLER_DETECTOR_VERSION=$DETECTOR_VERSION
export JUGGLER_DETECTOR_PATH=$DETECTOR_PATH
export JUGGLER_BEAMLINE_CONFIG=$BEAMLINE_CONFIG
export JUGGLER_BEAMLINE_CONFIG_VERSION=$BEAMLINE_CONFIG_VERSION
export JUGGLER_INSTALL_PREFIX=$ATHENA_PREFIX

## Export detector libraries
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${ATHENA_PREFIX}/lib

## modify PS1 for this detector version
export PS1="${PS1:-}"
export PS1="nightly${PS1_SIGIL}>${PS1#*>}"
unset branch
