# retina #######################################################################

A Clojure library designed to enable simulation of the responses of
retinal mosaics to natural images. This library was produced as part
of the paper "Unsupervised Learning of Cone Spectral Classes from
Natural Images" by Noah C. Benson, Jeremy R. Manning, and David
H. Brainard. It is intended for use as either a library or a
standalone tool.

## Authors #####################################################################

Primary Author: [[Noah C. Benson|mailto:n@nben.net]]

Principle Investigator: [[David H. Brainard|mailto:brainard@psych.upenn.edu]]

## Usage #######################################################################

This library is designed for use in simulating the responses of retinas to 
natural images. The expectation is that it will usually be used in
conjunction with Leiningen (leiningen.org), which provides a single
script, lein, that can be downloaded and run in order to resolve and
install all dependencies in this local directory. To do this, simply
follow these instructions:
 1. place the lein script (at the time of the writing of this README, the
    script can be found
    [[here|https://raw.github.com/technomancy/leiningen/stable/bin/lein]]
    for Mac or Unix and
    [[here|https://raw.github.com/technomancy/leiningen/stable/bin/lein.bat]]
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
   

## License #####################################################################

Copyright © 2013 Noah C. Benson

Distributed under the Eclipse Public License, the same as Clojure.
