# Analysis #####################################################################

This directory contains Mathematica and Matlab files related to
analyzing our simulations and plotting figures. These files generally
deal with simulation files; these are the files written by the clojure
code that is the bulk of this git repository. These files are written
by the clojure function
brainardlab.nben.retina.simulation/write-simulation. If you use the
main function interface or the simulate-some-retinas function (in the
brainardlab.nben.retina.util namespace), these will by default have
names that match the pattern sim[*].bin where the contents of the *
describe the simulation parameters.

The contents of this directory are detailed below:

#### retina-figures.nb

This Mathematica notebook contains code for reading in simulation
files, classifying cones, and plotting various analyses. All of the
figures in the paper can be made using this notebook.

#### readClojureSimFile.m

This Matlab function will read in a clojure simulation file and return
a matlab structure containing its data.

#### writeClojureSimFile.m

This Matlab function will write out a clojure simulation file given a
Matlab structure that contains its data; this structure must have the
same general shape as the structures read in by readClojureSimFile.m.

#### embedAll.m

This Matlab function will use the readClojureSimFile and
writeClojureSimFile to read all the sim*.bin files in a given
directory, construct a 3D embedding out of them using Matlab's mdscale
function, then write them out into the current directory (keeping the
same filenames aside from changing the directory). It safe to use in
the current directory, but we advise against overwriting original
simulation files. Any simulation file may have an embedding already
written into it; this is ignored if so.

## License #####################################################################

Copyright © 2013-2014 Noah C. Benson

Distributed under the Eclipse Public License, the same as Clojure.
