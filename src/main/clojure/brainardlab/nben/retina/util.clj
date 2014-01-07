;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; util.clj
;; Code for simulating a human retina as it responds to natural images; this code is specific to the
;; paper by Benson, Manning, and Brainard, for running simulations.
;; by Noah C. Benson

(ns brainardlab.nben.retina.util
  (:use brainardlab.nben.retina.constants)
  (:use brainardlab.nben.retina.core)
  (:use brainardlab.nben.retina.hyperspectral)
  (:use brainardlab.nben.retina.simulation)
  (:use clojure.contrib.generic.math-functions))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; These variables are "config" variables; they can be edited to change the parameters of everything
;; that works together below the next line of ;s

;; this is a list of simulation plans; it represents the kinds of simulations that get run based
;; on which simulation plan is requested
(def ^{:dynamic true
       :doc (str "simulation-plans is a variable that describes all the simulations that can be run"
                 " via the simulate-some-retinas function or the command line (via the main"
                 " function. It is marked as dynamic, so you may dynamically rebind it to add"
                 " simulation plans, or you may edit the code by hand (only way to change the"
                 " simulations that can be run from the command line. The variable itself is a map"
                 " of plan-id (a keyword) to a map of plan instructions. The instructions consist"
                 " of several optional keys that map to the values that the particular key should"
                 " take during simulation. When the plan executes, all combinations of the given"
                 " parameters are simulated:\n"
                 "  * :surrounds (default: [[0.25 3.0]]) is the list of parameters to pass to the"
                 " mosaic function when creating the retina. A value of nil indicates that no"
                 " surround should be used; [weight std] indicates that a surround with a total"
                 " given weight and with a standard deviation of std should be used. Additionally,"
                 " [weight std true] specifies that cone-specific surrounds should be used.\n"
                 "  * :L-to-Ms (default: [[1 1]]) is a list of the L:M ratios to use ([L M]).\n"
                 "  * :L-to-A-to-Ms (default: [[1 1 1]]) is the same as :L-to-Ms but for"
                 " simulations of tetrachromatic mosaics.\n"
                 "  * :L-lambda-maxs :M-lambda-maxs, :A-lambda-maxs, :S-lambda-maxs (defaults:"
                 " [558.9], [530.3], [545.0], [420.7], respectively) indicate the lambda-max values"
                 " to simulate for each cone class.\n"
                 "  * :sizes (default: [20]) is the retina size to simulate; this may be given as a"
                 " number, in which case the retina is square, or as [width height].\n"
                 "  * :runs (default: [0]) is the identifiers for the runs of these simulations."
                 " This parameter may be used to force separate runs of simulations that have"
                 " identical parameters.\n"
                 "  * :save-every (default: nil), if an integer, indicates that the simulation"
                 " should save a copy of the simulation so far (including the full correlation"
                 " matrix up to that point) every <save-every> images shown.\n"
                 "  * :S-cone-flags (default: [:human]) indicates flags that can be given to the"
                 " S-cones when creating the retina. If this is :human, then the S-cones are evenly"
                 " spaced across the retina and are held at 6% of the total cones. If this is :none"
                 " then the retina excludes S-cones. If this is a number or a map, the argument is"
                 " passed directly to the label-cones function (see"
                 " (doc brainardlab.nben.retina.core/label-cones) for more information).\n"
                 "  * :samples (default: [2500000]) specifies the number of samples to show the"
                 " retinas.\n"
                 "  * :spectral-indices (default: [nil]) is, if a collection is given, the indices"
                 " (in the hyperspectra) that should be shown one after the other during the"
                 " simulation as if a film that blocked all but that wavelength of light were"
                 " placed over the retina during simulation. This was included to simulate"
                 " conditions reported by Sugita et al (2004).\n"
                 "  * :noises (default: [nil]): if this is a number s between 0 and 1, includes"
                 " noise in the simulation at the level of cone responses such that the noise is"
                 " drawn from a standard deviation with mean 0 and standard deviation of s*u"
                 " where u is the mean of the cone responses.\n"
                 "  * :types (default: [:trichromat]) is the type of retina to assume; this may"
                 " include :trichromat, :dichromat, or :tetrachromat.")}
  simulation-plans
  {:basic     {:surrounds [[0.25 3.0]]
               :L-to-Ms [[16 1] [8 1] [4 1] [2 1] [1 1] [1 2] [1 4] [1 8] [1 16]]
               :M-lambda-maxs [530.3 535 540 545 550 555]
               :L-lambda-maxs [558.9]
               :S-lambda-maxs [420.7]
               :sizes [20]
               :S-cone-flags [:human]
               :samples [2500000]
               :runs [1 2]
               :save-every 100000}
   :standard  {:surrounds [[0.25 3.0] nil]
               :L-to-Ms [[16 1] [8 1] [4 1] [2 1] [1 1] [1 2] [1 4] [1 8] [1 16]]
               :M-lambda-maxs [530.3 535 540 545 550 555]
               :L-lambda-maxs [558.9]
               :S-lambda-maxs [420.7]
               :sizes [20 15 10 5]
               :S-cone-flags [:human]
               :samples [2500000]
               :runs [0]
               :save-every 100000}
   :blur      {:surrounds [[0.25 3.0] nil]
               :L-to-Ms [[16 1] [8 1] [4 1] [2 1] [1 1] [1 2] [1 4] [1 8] [1 16]]
               :M-lambda-maxs [530.3 535 540 545 550 555]
               :L-lambda-maxs [558.9]
               :S-lambda-maxs [420.7]
               :sizes [20]
               :S-cone-flags [:human]
               :samples [2500000]
               :runs ["blur_4"]
               :save-every 100000}
   :noise5    {:surrounds [[0.25 3.0]]
               :L-to-Ms [[16 1] [8 1] [4 1] [2 1] [1 1] [1 2] [1 4] [1 8] [1 16]]
               :M-lambda-maxs [530.3 535 540 545 550 555]
               :L-lambda-maxs [558.9]
               :S-lambda-maxs [420.7]
               :sizes [20]
               :S-cone-flags [:human]
               :samples [2500000]
               :noises [0.05]
               :runs [0]
               :save-every 100000}
   :noise1    {:surrounds [[0.25 3.0]]
               :L-to-Ms [[16 1] [8 1] [4 1] [2 1] [1 1] [1 2] [1 4] [1 8] [1 16]]
               :M-lambda-maxs [530.3 535 540 545 550 555]
               :L-lambda-maxs [558.9]
               :S-lambda-maxs [420.7]
               :sizes [20]
               :S-cone-flags [:human]
               :samples [2500000]
               :noises [0.01]
               :runs [0]
               :save-every 100000}
   :opponent  {:surrounds [[0.25 3.0 true]]
               :L-to-Ms [[16 1] [8 1] [4 1] [2 1] [1 1] [1 2] [1 4] [1 8] [1 16]]
               :M-lambda-maxs [530.3 535 540 545 550 555]
               :L-lambda-maxs [558.9]
               :S-lambda-maxs [420.7]
               :sizes [20]
               :S-cone-flags [:human]
               :samples [2500000]
               :runs ["opponent"]
               :save-every 100000}
   :Sugita    {:surrounds [[0.25 3.0]]
               :L-to-Ms [[1 1]]
               :M-lambda-maxs [530.3]
               :L-lambda-maxs [558.9]
               :S-lambda-maxs [420.7]
               :sizes [20]
               :S-cone-flags [:human]
               :samples [2500000]
               :spectral-indices [[6 12 19 24]]
               :runs ["Sugita2004"]
               :save-every 100000}
   :ditetra   {:surrounds [[0.25 3.0] nil]
               :L-to-A-to-Ms [[4 1 1] [2 2 1] [2 1 2] [1 1 1] [1 2 2] [1 1 4]]
               :M-lambda-maxs [530.3]
               :A-lambda-maxs [535.0 540.0 545.0 550.0 555.0]
               :L-lambda-maxs [558.9]
               :S-lambda-maxs [420.7]
               :sizes [20]
               :S-cone-flags [:human]
               :types [:dichromat :tetrachromat]
               :samples [2500000]
               :runs [0]
               :save-every 100000}
   :periphery {:surrounds [[0.25 3.0]]
               :L-to-Ms [[16 1] [8 1] [4 1] [2 1] [1 1] [1 2] [1 4] [1 8] [1 16]]
               :M-lambda-maxs [535.0 540.0 545.0 550.0 555.0]
               :L-lambda-maxs [558.9]
               :S-lambda-maxs [420.7]
               :sizes [40]
               :spacings [2]
               :S-cone-flags [:human]
               :samples [2500000]
               :runs ["periphery"]
               :save-every 100000}
   :spectral-spacing
              {:surrounds [[0.25 3.0]]
               :L-to-A-to-Ms [[2 2 1] [2 1 2] [1 2 2] [1 1 1] [1 1 2] [1 2 1] [2 1 1]]
               :L-lambda-maxs [558.9]
               :M-lambda-maxs [512.8]
               :A-lambda-maxs [466.8]
               :S-lambda-maxs [420.7]
               :types [:tetrachromat]
               :sizes [20]
               :S-cone-flags [:human]
               :samples [2500000]
               :runs ["spectral_spacing"]
               :save-every 100000}
   :LM-only   {:surrounds [[0.25 3.0]]
               :L-to-Ms [[16 1] [8 1] [4 1] [2 1] [1 1] [1 2] [1 4] [1 8] [1 16]]
               :M-lambda-maxs [530.3]
               :L-lambda-maxs [558.9]
               :S-lambda-maxs [420.7]
               :sizes [20]
               :S-cone-flags [:none]
               :samples [5000000]
               :runs ["tritanope"]
               :save-every 100000}})

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; create a list of retinas from a simulation plan
(defn retinas-from-plan [plan]
  (let [cone-lambda-maxs
        (for [M-lambda-max (get plan :M-lambda-maxs [530.3])
              L-lambda-max (get plan :L-lambda-maxs [558.9])
              S-lambda-max (get plan :S-lambda-maxs [420.7])
              A-lambda-max (get plan :A-lambda-maxs [545.0])]
          {:L L-lambda-max :M M-lambda-max :S S-lambda-max :A A-lambda-max})
        mosaics
        (for [retina-size (get plan :sizes [20])
              spacing (get plan :spacings [1])
              blur (get plan :blurs [nil])
              specindices (get plan :spectral-indices [nil])]
          (mosaic :size retina-size :blur blur :spacing spacing :spectral-indices specindices))]
    (for [retina-type (get plan :types [:trichromat])
          surround (get plan :surrounds [[0.25 3.0]])
          noise (get plan :noises [nil])
          L-to-M (cond (= retina-type :trichromat) (get plan :L-to-Ms [[1 1]])
                       (= retina-type :tetrachromat) (get plan :L-to-A-to-Ms [[1 1 1]])
                       (= retina-type :dichromat) [1]
                       :else (throw (IllegalArgumentException. "invalid type")))
          S-cone-flags (get plan :S-cone-flags [:human])
          mosaic mosaics
          lmaxs (loop [q cone-lambda-maxs, res #{}]
                  (if (nil? q)
                    (if (= S-cone-flags :none)
                      (map #(dissoc % :S) (seq res))
                      (seq res))
                    (recur (next q)
                           (conj res
                                 (cond (= retina-type :trichromat) (dissoc (first q) :A)
                                       (= retina-type :tetrachromat) (first q)
                                       (= retina-type :dichromat) (dissoc (first q) :A :M))))))]
      (retina :type retina-type
              :surround surround
              :noise noise
              :mosaic mosaic
              :cones (let [start-cones (cond (= retina-type :trichromat) [:L :M :S]
                                             (= retina-type :tetrachromat) [:L :A :M :S]
                                             (= retina-type :dichromat) [:L :S])]
                       (if (= :none S-cone-flags)
                         (butlast start-cones)
                         start-cones))
              :cone-makeup (let [allbutS (cond (= retina-type :trichromat)
                                               {:L (first L-to-M)
                                                :M (fnext L-to-M)}
                                               (= retina-type :tetrachromat)
                                               {:L (nth L-to-M 0)
                                                :A (nth L-to-M 1)
                                                :M (nth L-to-M 2)}
                                               :else {:L 1})]
                             (cond (= S-cone-flags :human)
                                   (assoc allbutS :S {:even-spacing true, :fraction 0.06})
                                   (= S-cone-flags :none)
                                   allbutS
                                   (or (map? S-cone-flags)
                                       (and (number? S-cone-flags) (> S-cone-flags 0)))
                                   (assoc allbutS :S S-cone-flags)
                                   :else (throw (IllegalArgumentException. "invalid s-cone flag"))))
              :lambda-max (cond (= S-cone-flags :none) (dissoc lmaxs :S)
                                :else lmaxs)))))

;; yields a filename for writing out a simulation
(defn simulation-filename [r image-count run-id]
  (let [params (:params r)
        Mlm (if-let [m (:M (let [lm (:lambda-max params)]
                             (if (= lm :automatic)
                               human-lambda-max
                               lm)))]
              (Math/round (double m))
              "none")
        cone-makeup (:cone-makeup params)
        L-M (str (:L cone-makeup) ":" (or (:M cone-makeup) 0))
        surr (let [ss (:surround params)
                   s (if (= ss :automatic) [0.25 3.0] ss)]
               (if s
                 (if (number? s)
                   (str s "_" 3.0)
                   (str (first s) "_" (fnext s)))
                 "none"))
        noise (if-let [n (:noise params)] (str ",noise=" n) "")
        mosaic (:mosaic r)
        mparams (:params mosaic)
        sz (:size mosaic)
        mstr (str (if (= :rectilinear (:layout mparams)) "r" "h")
                  (if (number? sz) (str sz "x" sz) (str (first sz) "x" (fnext sz)))
                  "_" (:separation mparams))]
    (str "sim[mosaic=" mstr
         ",L:M=" L-M
         ",Mlm=" Mlm
         ",surround=" surr
         noise
         ",samples=" image-count
         ",run=" run-id
         (let [tt (:type (:params r))]
           (cond (= tt :tetrachromat)
                 (str ",tetrachromat[" (:A cone-makeup)
                      "," (Math/round (double (:A (:lambda-max params)))) "]")
                 (= tt :dichromat) ",dichromat"
                 :else ""))
         "].bin")))

;; used when setting up retina simulations automatically by simulate-some-retinas
(defn- analysis-print-every [[analfn finfn every-k iter dat] sig]
  (if (= 0 (mod iter every-k))
    (println (format "Iteration %8d..." iter)))
  [analfn finfn every-k (inc iter) (analfn dat sig)])
(defn- finish-print-every [[analfn finfn every-k iter dat]]
  (println "Simulation finished after " iter " iterations.")
  (finfn dat))
(defn- analysis-save-every [[analfn finfn R run embed every-k iter outdir dat] sig]
  (if (and (> iter 0) (= 0 (mod iter every-k)))
    (write-simulation (str outdir "/" (simulation-filename R iter run)) R (finfn dat) :embed embed))
  [analfn finfn R run embed every-k (inc iter) outdir (analfn dat sig)])
(defn- finish-save-every [[analfn finfn R run embed every-k iter outdir dat]]
  (finfn dat))
(defn- init-print-every [analfn finfn every-k start-dat]
  [analfn finfn every-k 0 start-dat])
(defn- init-save-every [analfn finfn R run embed every-k outdir start-dat]
  [analfn finfn R run embed every-k 0 outdir start-dat])
(defn- make-analysis [save-every outdir print-every Rs run embed]
  (let [n (count Rs)
        [init1 anal1 fin1]
        (if save-every
          [(map #(init-save-every correlation-analysis correlation-analysis
                                  %1 run embed save-every outdir nil)
                Rs)
           (repeat n analysis-save-every)
           (repeat n finish-save-every)]
          [(repeat n nil) (repeat n correlation-analysis) (repeat n correlation-analysis)])]
    (if print-every
      [(cons (init-print-every (first anal1) (first fin1) print-every (first init1))
             (next init1))
       (cons analysis-print-every (next anal1))
       (cons finish-print-every (next fin1))]
      [init1 anal1 fin1])))

;; a function for simulating some of the retinas; given the total number of jobs and
;; a a job id, this will simulate and write out the results for the appropriate chunk
;; of simulations from the given set. The plan should be :basic or :standard.
(defn simulate-some-retinas [plan-id hs-cache-filename output-dir total-nodes node-id
                             &{:keys [verbose print-every] :or {print-every 10000}}]
  (let [plan (get simulation-plans plan-id)
        plan-retinas (retinas-from-plan plan)
        my-retinas (if (coll? node-id)
                     (apply concat
                            (map #(take-nth total-nodes (nthnext plan-retinas %))
                                 node-id))
                     (take-nth total-nodes (nthnext plan-retinas node-id)))
        runs (get plan :runs [0])
        cache (read-hyperspectral-cache hs-cache-filename)
        samples (get plan :samples [2500000])
        embed (get plan :embed true)
        save-every (get plan :save-every)]
    (if (nil? plan) (throw (IllegalArgumentException. "Plan not found")))
    (doseq [image-count samples
            run runs]
      (if verbose
        (println (format "Simulating %d retinas with %d samples from cache %s (run = %s)..."
                         (count my-retinas) image-count hs-cache-filename (str run))))
      (let [[inits analfns finfns] (make-analysis save-every
                                                  output-dir
                                                  (when verbose print-every)
                                                  my-retinas
                                                  run
                                                  embed)
            corr-mtcs (simulate-retinas my-retinas cache
                                        :image-count image-count
                                        :init-reduce inits
                                        :analysis analfns
                                        :finish finfns)
            flnms (map #(str output-dir "/" (simulation-filename % image-count run)) my-retinas)]
        (doall
         (map #(write-simulation %1 %2 %3 :embed embed)
              flnms my-retinas corr-mtcs))))))

    
        
