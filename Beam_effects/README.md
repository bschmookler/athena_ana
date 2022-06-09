The files in this directory will first allow you to run a Djangoh simulation for the 18x275GeV beam energy setting. Djangoh generates the events in the colinear frame with fixed beam energy and interaction vertex. The output of this simulation is then used as input to our beam effects afterburner. The output of the afterburner can then be put through the detector simulation.

Instructions for running:

1) To run the generator, you need to leave the ATHENA container, and set the environment as described in the instruction included in this [repository](https://github.com/cipriangal/eicGenTutorials).

2) Run the generator:
```
./run_ep_norad.sh
```
This will create a HepMC3 file.

3) The beam effects utility can be found [here](https://eicweb.phy.anl.gov/monte_carlo/afterburner). You can download and compile it yourself if you wish, but it is already available on the ATHENA EIC container. So go back into the original container and run the following script:
```
./run_afterburner.sh
```
This will create a new HepMC3 file which can be used as input to the detector simulation.


