# athena_ana
<br/>

ATHENA Software Container
-------------------------
As described [here](https://eic.phy.anl.gov/tutorials/eic_tutorial/getting-started/quickstart/#step-1-setup-the-eic-software-container-jug_xl), run the following commands to enter the software container:

mkdir eic

cd eic

curl https://eicweb.phy.anl.gov/containers/eic_container/-/raw/master/install.sh | bash

./eic-shell

You have now entered into the container and can begin work.
<br/>

Clone Repositories
------------------
As discussed here(https://eic.phy.anl.gov/tutorials/eic_tutorial/getting-started/quickstart#step-2-clone-the-repos) , first clone the repositories for the beamline and the athena detector: 

git clone https://eicweb.phy.anl.gov/EIC/detectors/athena.git

git clone https://eicweb.phy.anl.gov/EIC/detectors/ip6.git

ln -s ../ip6/ip6 athena/ip6
