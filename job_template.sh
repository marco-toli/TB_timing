echo $SHELL
pwd
source /afs/cern.ch/sw/lcg/external/geant4/10.1/x86_64-slc6-gcc49-opt/bin/geant4.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.26/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
source /afs/cern.ch/sw/lcg/external/gcc/4.9.1/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/external/geant4/10.1/x86_64-slc6-gcc49-opt/CMake-setup.sh
unbuffer /afs/cern.ch/work/m/mlucchin/TB_timing/build/TB_timing /afs/cern.ch/work/m/mlucchin/TB_timing/templates/template.cfg out_ID
cp out_ID.root /afs/cern.ch/work/m/mlucchin/TB_timing/jobs/job_ID