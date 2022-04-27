#!/bin/sh
export S3_ACCESS_KEY=eicS3read
export S3_SECRET_KEY=eicS3read

## For using 'official' setup, uncomment line below and comment all following lines
#source /opt/detector/setup.sh

export DETECTOR=athena
## If you use the 'share' version, you will need to make a softlink to ip6
#export DETECTOR_PATH=${ATHENA_PREFIX}/share/athena
## May want to use source directory instead
export DETECTOR_PATH=${ATHENA_PREFIX}/../athena

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
export PS1="local${PS1_SIGIL}>${PS1#*>}"
unset branch
