#!/bin/bash

echo "-----------------------------------"
echo "Running Beam Effects Afterburner to DJANGOH Simulation!"
echo "..."
echo ""

abconv outfiles/djangoh.NC.18x275_evt.hepmc -o outfiles/djangoh.NC.18x275_beameffects --plot-off

echo "Done"


