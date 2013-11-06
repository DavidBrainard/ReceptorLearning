;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; main.clj
;; Code for simulating a human retina as it responds to natural images; this code is specific to the
;; paper by Benson, Manning, and Brainard, for running simulations. Basically this file wraps a main
;; function around the simulate-some-retinas function.
;; by Noah C. Benson

(ns brainardlab.nben.retina.main
  (:use brainardlab.nben.retina.constants)
  (:use brainardlab.nben.retina.core)
  (:use brainardlab.nben.retina.hyperspectral)
  (:use brainardlab.nben.retina.simulation)
  (:use brainardlab.nben.retina.util)
  (:gen-class))

(defn- die [& txts] (println (apply str txts)) (System/exit 0))

(defn -main [& args]
  (if (not= (count args) 5)
    (die "Arguments required: plan-id cache-filename output-dir number-of-workers worker-id")
    (let [[plan-id-str cache-flnm output-dir number-of-workers-str worker-id-str] args
          plan-id (clojure.lang.Keyword/intern plan-id-str)
          number-of-workers (Integer/parseInt number-of-workers-str)
          worker-id (Integer/parseInt worker-id-str)]
      (simulate-some-retinas plan-id
                             cache-flnm
                             output-dir
                             number-of-workers
                             worker-id
                             :verbose true)))
  (System/exit 0))

