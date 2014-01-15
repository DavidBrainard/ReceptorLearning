# Receptor Learning ############################################################

A Clojure library designed to enable simulation of the responses of
retinal mosaics to natural images. This library was produced as part
of the paper "Unsupervised Learning of Cone Spectral Classes from
Natural Images" by Noah C. Benson, Jeremy R. Manning, and David
H. Brainard. It is intended for use as either a library or a
standalone tool.

## Authors #####################################################################

Primary Author: [Noah C. Benson](mailto:n@nben.net)

Principle Investigator: [David H. Brainard](mailto:brainard@psych.upenn.edu)

## Organization ################################################################

This library is organized into a few directories, detailed here:
 1. analysis - The analysis directory contains all the code needed to
    post-process the simulations and plot the figures from the
    paper. This directory inludes both Matlab code (for embedding the
    correlation matrices) and Mathematica code (for flattening and
    classifying retinas and for generating figures).
 2. doc - The doc directory contains html documentation for the
    clojure library. This html may also be access directly
    [here](http://davidbrainard.github.io/ReceptorLearning/).
    The Matlab and Mathematica code in the analysis
    directory is documented internally.
 3. samples - The samples directory contains examples of natural image
    patches that were drawn from our hyperspectral database for use in
    simulation of a 20 x 20 retina.
 4. src - The src (source) directory contains clojure source code for
    the receptor learning retina simulation library. It is documented
    in the doc directory.
 5. resources - The resources directory contains additional resources
    that are not required for use in the library, but that may be useful.

## Usage #######################################################################

This library is designed for use in simulating the responses of retinas to 
natural images. The expectation is that it will usually be used in
conjunction with [Leiningen](http://leiningen.org), which provides a single
script, lein, that can be downloaded and run in order to resolve and
install all dependencies in this local directory. To do this, simply
follow these instructions:
 1. place the lein script (at the time of the writing of this README, the
    script can be found
    [here](https://raw.github.com/technomancy/leiningen/stable/bin/lein)
    for Mac or Unix and
    [here](https://raw.github.com/technomancy/leiningen/stable/bin/lein.bat)
    for Windows), in your path,
 2. clone this git repository,
 3. in a terminal, switch to the project's root directory (the
    directory containing the project.clj file),
 4. enter 'lein repl'.

This should result in a number of dependencies being resolved followed
by a prompt. This prompt is a clojure repl; to import the library,
simply enter:
<code>
(use '(brainardlab.nben.retina core constants hyperspectral
                               simulation util))
</code>
You may then access functions such as simulate-some-retinas directly.

### Using the main function ####################################################

Alternately, you may instead enter 'lein run' in place of 'lein repl';
this is the interface for calling the main function in the
library. The main function expects 4 arguments: the name of the plan
to simulate (these are listed in util.clj or 
brainardlab.nben.retina.util/simulation-plans, which is
well-documented; note that these are keywords in the code, but this
argument expects a string without the preceding :, as in "tritanopes",
not ":tritanopes"), the hyperspectral cache file, the output directory
(into which simulation files will be written), the number of workers
processing the plan, and the id of this worker. If you are running the
simulations on a single computer, then you would be using 1 worker and
that worker's id would be 0. For distributed calculation of a plan
(such as the rather large :standard plan), if the number of workers is
n and each worker is given an id 0 ... n-1, then the plan will be
divided up and executed as evenly as possible.

### Hyperspectral-cache files ##################################################

The simulations expect to load their images from a hyperspectral
cache. This is effectively a pre-cached set of natural images, all of
which have had a fair amount of pre-processing done so as to optimize
the speed of simulations. Hyperspectral cache files can be built using
the functions in the hyperspectral namespace
(brainardlab.nben.retina.hyperspectral). The relevant functions are
build-hyperspectral-cache (which requires that images be added
sequentially to the cache using the push-to-cache function) and
autobuild-hyperspectral-cache. The latter is recommended, as it
requires only the paths of the databases used in Benson et
al. (2014), and will automatically load the images and build the cache
from them.  For example:

     (autobuild-hyperspectral-cache "data/my-hyperspectral-cache.bin"
                                    :Chakrabarti2011 "hsdbs/Chakrabarti2011/CZ_hsdb"
                                    :Foster2004 "hsdbs/Foster2004"
                                    :Parragas1998 "hsdbs/Parragas1998/brelstaff")

### Multi-threaded computations ################################################

Note that this entire library was written in clojure partially because
clojure's paradigm is optimized for multi-threaded programming. This
library automatically threads across as many processors as possible
when performing simulations, and it is much much more efficient to run
simulations in bulk. That said, there are some perils. Primarily, if
you simulate a very small retina with a very large retina, the small
retina will finish scanning the entire hyperspectral cache before the
large retina has advanced very far. This can result in memory
problems, as the library holds in memory all hyperspectral images in
the hyperspectral cache that have (1) been seen by at least one retina
and (2) have not been seen by all retinas in the simulation. Because
the retinas scan the cache in order, one ideally wants them to be
approximately in sync with each other.

## Documentation ###############################################################

Documentation is available in the doc directory of the
repository. This documentation was generated with
[Codox](https://github.com/weavejester/codox), and provides
descriptions of all the public functions in the package. This
documentation can also be found online
[here](http://davidbrainard.github.io/ReceptorLearning/index.html).

## License #####################################################################

Copyright © 2013-2014 Noah C. Benson

Distributed under the Eclipse Public License, the same as Clojure.
