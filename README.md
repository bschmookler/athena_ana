# athena_ana
<br/>

ATHENA Software Container
-------------------------
As described [here](https://eic.phy.anl.gov/tutorials/eic_tutorial/getting-started/quickstart/#step-1-setup-the-eic-software-container-jug_xl), run the following commands to enter the software container:

```
mkdir eic
cd eic
curl https://eicweb.phy.anl.gov/containers/eic_container/-/raw/master/install.sh | bash
./eic-shell
```

You have now entered into the container and can begin work.
<br/>

Using the default detector and reconstruction setup
---------------------------------------------------
In the container, first type
```
source /opt/detector/setup.sh
```

Running the simulation consists of 2 steps. First, generated particles are passed through the detector and an output ROOT file is created; second, the output ROOT file of the detector simulation step is used as the input to the reconstruction software.
<br/>

To throw 100 single muon events through the detector, do the following:
```
npsim --compactFile $DETECTOR_PATH/athena.xml --enableGun --gun.distribution uniform --numberOfEvents 100 --outputFile output.edm4hep.root
```

To instead generated 100 single particle pi+ events, do the following:
```
npsim --compactFile $DETECTOR_PATH/athena.xml --enableGun --gun.distribution uniform --gun.particle pi+ --numberOfEvents 100 --outputFile output.edm4hep.root
```

If you have the output of an event generator that is saved in HEPMC format, then do the following to run events from that generator through the simulation:
```
npsim --compactFile $DETECTOR_PATH/athena.xml --numberOfEvents 25 --inputFiles input.hepmc --outputFile output.edm4hep.root
```

In all these cases, an output file called <i>output.edm4hep.root</i> is created by the detector simulation. This file contains information on the generated and secondary particles, the hit positions in the tracking detectors, the hit positions and energy deposited in the calorimeter cells, and hit information for other detectors.
<br/>

To then run the output of the simulation through the reconstruction framework (which is called Juggler), do the following:

```
export JUGGLER_SIM_FILE=output.edm4hep.root JUGGLER_REC_FILE=rec_output.edm4hep.root JUGGLER_N_EVENTS=10
gaudirun.py /opt/benchmarks/physics_benchmarks/options/reconstruction.py
```

This creates a file called <i>rec_output.edm4hep.root</i>, which contains information on the reconstructed tracks, calorimeter clusters, etc...

If you want to focus on just the calorimeters, for example, you can run this reconstruction script:
```
export JUGGLER_SIM_FILE=output.edm4hep.root JUGGLER_REC_FILE=rec_cal_output.edm4hep.root JUGGLER_N_EVENTS=10
gaudirun.py /opt/benchmarks/reconstruction_benchmarks/benchmarks/clustering/options/full_cal_reco.py
```

Using local builds of the detector and reconstruction software
---------------------------------------------------------------
In order to build local copies of the detector, beamline, and reconstruction software, as discussed in part [here](https://eic.phy.anl.gov/tutorials/eic_tutorial/getting-started/quickstart#step-2-clone-the-repos), first clone the repositories for the beamline, the athena detector, and the reconstruction software: 

```
git clone https://eicweb.phy.anl.gov/EIC/detectors/athena.git
git clone https://eicweb.phy.anl.gov/EIC/detectors/ip6.git
ln -s ../ip6/ip6 athena/ip6
git clone https://eicweb.phy.anl.gov/EIC/juggler.git
```

Now compile each of these one at a time.

For the ATHENA detector:
```
cd athena
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$ATHENA_PREFIX
make
make install
```

For the beamline:
```
cd ip6
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$ATHENA_PREFIX
make
make install
```

For the reconstruction software:
```
cd juggler
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$ATHENA_PREFIX
make
make install
```

Doing the above will install a bunch of things in your <i>$ATHENA_PREFIX</i> directory. In order to use your local install, instead of sourcing the default <i>setup.sh</i> file, source the [mysetup.sh](mysetup.sh) script included with this repository. Then you can simulate events with the same commands as above, but now using your local install of the detector, beamline, and reconstruction software. 
<br/>

Analyzing the simulation output
-------------------------------
Take a look at the [Analysis_examples](Analysis_examples) folder in this repository to see some examples of how the analyze the output of the simulation. Both ROOT-based macros and an Uproot code are included as examples.
<br/>

Standalone detector example
---------------------------
An example of a simple standalone detector can be found [here](Detector_examples/samplinghcal). As this detector is implemented using a standard class in the ATHENA framework (<i>ffi_ZDC_Sampling</i>), we do not need to compile anything. Begin by loading the environment
```
source /opt/detector/setup.sh
```
and then run single particles through the detector as
```
npsim --runType run --enableG4GPS --macroFile gps.mac --compactFile samplinghcal.xml --outputFile sim_out.root
```

Note how this simulation uses the file [gps.mac](Detector_examples/samplinghcal/gps.mac) to define the particle generation and number of events to simulate.

To view the geometry of this standalone dectector, first run this command:
```
dd_web_display --export samplinghcal.xml
```
and then upload the created ROOT file to this [page](https://eic.phy.anl.gov/geoviewer/). You should be able to see the following picture:
![detector_geometry](Detector_examples/samplinghcal/samplinghcal_geometry.png?raw=true)
<br/>

Generation with Beam Effects
----------------------------
To include the EIC beam effects to generated physics events using this [afterburner](https://eicweb.phy.anl.gov/monte_carlo/afterburner), see this [directory](Beam_effects).

General Documentation
---------------------
The official documentation, including various tutorials, can be found [here](https://eic.phy.anl.gov/tutorials/eic_tutorial/).
<br/>
