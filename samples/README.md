# retina samples directory #####################################################

This directory contains sample natural image patches used in our
simulations. These patches were generated with the following code
block, executed from the project root directory in the lein repl:

        (use '(brainardlab.nben.retina constants core hyperspectral simulation util))
        (map (fn [idx hsimage]
               (write-hyperspectral-as-png (format "samples/patch%02d.png" idx)
                                           hsimage))
             (range)
             (draw-patches (read-hyperspectral-cache "data/cache.bin") 100 20 20))

## License #####################################################################

Copyright © 2013-2014 Noah C. Benson

Distributed under the Eclipse Public License, the same as Clojure.
