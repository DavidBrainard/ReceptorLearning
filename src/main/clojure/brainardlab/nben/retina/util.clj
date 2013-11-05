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
(def simulation-plans
     {:basic    {:surrounds [[0.25 3.0]]
                 :L-to-Ms [[16 1] [8 1] [4 1] [2 1] [1 1] [1 2] [1 4] [1 8] [1 16]]
                 :M-lambda-maxs [530.3 535 540 545 550 555]
                 :L-lambda-maxs [558.9]
                 :S-lambda-maxs [420.7]
                 :sizes [20]
                 :S-cone-flags [:human]
                 :samples [500000]
                 :runs [1 2]
                 :caches ["/Users/iron/nben/data/hyperspectral-cache.bin"]
                 :output-dir "/Users/iron/nben/data/simulations/basic"}
      :standard {:surrounds [[0.25 3.0] nil]
                 :L-to-Ms [[16 1] [8 1] [4 1] [2 1] [1 1] [1 2] [1 4] [1 8] [1 16]]
                 :M-lambda-maxs [530.3 535 540 545 550 555]
                 :L-lambda-maxs [558.9]
                 :S-lambda-maxs [420.7]
                 :sizes [20 15 10 5]
                 :S-cone-flags [:human]
                 :samples [2500000]
                 :runs [0]
                 :caches ["/Users/iron/nben/data/hyperspectral-cache.bin"]
                 :output-dir "/Users/iron/nben/data/simulations/standard"}})

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; create a list of retinas from a simulation plan
(defn retinas-from-plan [plan]
  (for [surround (get plan :surrounds [[0.25 3.0]])
        L-to-M (get plan :L-to-Ms [[1 1]])
        M-lambda-max (get plan :M-lambda-maxs [530.3])
        L-lambda-max (get plan :L-lambda-maxs [558.9])
        S-lambda-max (get plan :L-lambda-maxs [420.7])
        S-cone-flags (get plan :S-cone-flags [:human])
        retina-size (get plan :sizes [20])]
    (retina :surround surround
            :mosaic (mosaic :size retina-size)
            :cones [:L :M :S]
            :cone-makeup (cond (= S-cone-flags :human)
                               {:L (first L-to-M)
                                :M (fnext L-to-M)
                                :S {:even-spacing true, :fraction 0.06}}
                               (= S-cone-flags :none)
                               {:L (first L-to-M)
                                :M (fnext L-to-M)}
                               :else (throw (IllegalArgumentException. "invalid s-cone flag")))
            :lambda-max {:L L-lambda-max, :M M-lambda-max, :S S-lambda-max})))

;; yields a filename for writing out a simulation
(defn simulation-filename [r image-count run-id]
  (let [params (:params r)
        Mlm (Math/round (double (:M (let [lm (:lambda-max params)]
                                      (if (= lm :automatic)
                                        human-lambda-max
                                        lm)))))
        cone-makeup (:cone-makeup params)
        L-M (str (:L cone-makeup) ":" (:M cone-makeup))
        surr (let [ss (:surround params)
                   s (if (= ss :automatic) [0.25 3.0] ss)]
               (if s
                 (if (number? s)
                   (str s "_" 3.0)
                   (str (first s) "_" (fnext s)))
                 "none"))
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
         ",samples=" image-count
         ",run=" run-id
         (let [tt (:type (:params r))]
           (cond (= tt :tetrachromat) ",tetrachromat"
                 (= tt :dichromat) ",dichromat"
                 :else ""))
         "].bin")))

;; a function for simulating some of the retinas; given the total number of jobs and
;; a a job id, this will simulate and write out the results for the appropriate chunk
;; of simulations from the given set. The plan should be :basic or :standard.
(defn simulate-some-retinas [plan-id total-nodes node-id &{:keys [verbose]}]
  (let [plan (get simulation-plans plan-id)
        plan-retinas (retinas-from-plan plan)
        my-retinas (if (coll? node-id)
                     (apply concat
                            (map #(take-nth total-nodes (nthnext plan-retinas %))
                                 node-id))
                     (take-nth total-nodes (nthnext plan-retinas node-id)))
        runs (get plan :runs [0])
        caches (loop [r {}, q (:caches plan)]
                 (if q
                   (recur (assoc r (first q) (read-hyperspectral-cache (first q)))
                          (next q))
                   r))
        samples (get plan :samples [2500000])
        output-dir (get plan :output-dir ".")
        embed (get plan :embed true)]
    (doseq [cache caches
            image-count samples
            run runs]
      (if verbose
        (println (format "Simulating %d retinas with %d samples from cache %s (run = %d)..."
                         (count my-retinas) image-count (key cache) run)))
      (let [corr-mtcs (if verbose
                        (simulate-retinas
                         my-retinas (val cache)
                         :image-count image-count
                         :init-reduce (cons [0 nil] (repeat (dec (count my-retinas)) nil))
                         :analysis (cons (fn [[iter dat] sig]
                                           (if (= 9999 (mod iter 10000))
                                             (println "Iteration " iter "..."))
                                           [(inc iter) (correlation-analysis dat sig)])
                                         (repeat (dec (count my-retinas)) correlation-analysis))
                         :finish (cons (fn [[iter dat]]
                                         (println "Simulation complete after " iter " iterations.")
                                         (correlation-analysis dat))
                                       (repeat (dec (count my-retinas)) correlation-analysis)))
                        (simulate-retinas my-retinas (val cache) :image-count image-count))
            flnms (map #(str output-dir "/" (simulation-filename % image-count run)) my-retinas)]
        (doall
         (map #(write-simulation %1 %2 %3 :embed embed)
              flnms my-retinas corr-mtcs))))))

    
        
