;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; simulation.clj
;; This file contains code used in simulation of perceptual systems.
;; By Noah C. Benson

(ns brainardlab.nben.retina.simulation
  (:use brainardlab.nben.retina.core)
  (:use brainardlab.nben.retina.hyperspectral)
  (:use clojure.contrib.generic.math-functions)
  (:import mdsj.MDSJ))

(defn- arg-err [& txt]
  (throw (IllegalArgumentException. (apply str txt))))

(defn correlation-analysis
  "This function is used to perform correlation analysis via simulate-retina.  It requires a ref
   and an optional seq of signals from a retina.  If the seq is not provided, the function finalizes
   the calculation, whose temporary data is stored in the ref, and returns the final correlation
   matrix.  If the ref seq is passed, it updates the ref with new correlation information.  The ref
   may initially point to nil."
  ([ca] (ca))
  ([ca sig]
     (let [coranal (or ca (brainardlab.nben.retina.jvm.CorrelationAnalysis. (count sig)))]
       (coranal sig)
       coranal)))

(defn embed-correlation-matrix
  "Given a correlation matrix produced by the simulate-retina function in conjunction with the
   correlation-analysis mode (the default) or represented as a double[][], yields the non-metric
   multidimensional scaling embedding of the negative log of the data squared."
  [data dims]
  (let [n (if (instance? brainardlab.nben.retina.jvm.Matrix2D data)
            (first (data))
            (count data))
        mtx (make-array Double/TYPE n n)]
    (if (instance? brainardlab.nben.retina.jvm.Matrix2D data)
      (doseq [i (range n), j (range n)]
        (aset mtx i j (- (Math/log (double (* (data i j) (data i j)))))))
      (doseq [i (range n), j (range n)]
        (aset mtx i j (- (Math/log (double (* (nth (nth data i) j) (nth (nth data i) j))))))))
    (mdsj.MDSJ/stressMinimization mtx dims)))

(defn simulate-retinas
  "Runs the simulation for the given retinas using the given input database.  The retinas argument
   may be a single retina or a collection of retinas.  The hs-database argument must be an
   hs-cache file name or File object.  The following options may be used:
    :image-count may either be a positive integer or a map of retina to positive integer or a
       collection of positive integers in the same order as retinas.  If such a collection is given
       it may be longer or shorter than the list of retinas, but will cycle upon reaching the end.
       The image count specifies how many images to show to each retina.  The default is 25,000.
    :analysis may be nil or a function that takes a retina and a sequence of signal values from the
       input or a collection of such functions, following the general seq/map form as specified in
       :image-count."
  [retinas hs-database & {:keys [image-count analysis finish full-reduce init-reduce image-filter]
                          :or {image-count 25000 full-reduce true}}]
  (let [rets (or (and (map? retinas) (list retinas))
                 (and (coll? retinas) (every? map? retinas)
                      retinas)
                 (arg-err "retinas must be a retina or a collection of retinas"))
        nrets (count rets)
        cache (or (and (map? hs-database) hs-database)
                  (and (or (string? hs-database) (instance? java.io.File hs-database))
                       (read-hyperspectral-cache hs-database))
                  (arg-err "hs-database must be a cache, filename, or File object"))
        analfns (or (and (nil? analysis) (repeat nrets correlation-analysis))
                    (and (seq? analysis) (= nrets (count analysis)) (every? ifn? analysis)
                         analysis)
                    (and (ifn? analysis) (repeat nrets analysis))
                    (arg-err ":analysis must be a function or seq of functions"))
        finfns ((if (nil? analysis)
                  (fn [ffs] (map (fn [ff] #(ff (correlation-analysis %))) ffs))
                  identity)
                (or (and (nil? finish) (repeat nrets identity))
                    (and (seq? finish) (= nrets (count finish)) (every? ifn? finish)
                         finish)
                    (and (ifn? finish) (repeat nrets finish))
                    (and finish (arg-err ":finish must be a function"))))]
    ;; get signals and feed them to the analysis function
    (doall
     (if (or full-reduce (nil? analysis))
       (pmap (fn [r analfn finfn init-red]
               (finfn
                (reduce (fn [dat sig] (analfn dat sig))
                        init-red
                        (let [[pheight pwidth] (:size (:mosaic r))]
                          (map #(cone-responses r %)
                               (if (ifn? image-filter)
                                 (map image-filter (draw-patches cache image-count pheight pwidth))
                                 (draw-patches cache image-count pheight pwidth)))))))
             rets
             analfns
             finfns
             (or (and (seq? init-reduce) init-reduce) (repeat nrets init-reduce)))
       (map (fn [r analfn finfn]
              (finfn
               (analfn (let [[pheight pwidth] (:size (:mosaic r))]
                         (map #(cone-responses r %)
                              (if (ifn? image-filter)
                                (map image-filter (draw-patches cache image-count pheight pwidth))
                                (draw-patches cache image-count pheight pwidth)))))))
             rets
             analfns
             finfns)))))

;; euclidean distance
(defn- ED [a b] (Math/sqrt (reduce + (let [dif (map - a b)] (map * dif dif)))))
(defn- ED2 [a b] (reduce + (let [dif (map - a b)] (map * dif dif))))
(defn- kohonen-analysis [signals neurons memory nmem iter Rfn afn save-every]
  (if (> (apply max (flatten (map #(or (fnext %) 0) neurons))) 2)
    (throw (Exception. (str "!!! " [iter (apply max (flatten (map fnext neurons)))]))))
  (cond
   (nil? memory)
   ;; we start by making a seq for memory, then we recur
   (let [mem (map #(into-array Double/TYPE %) (apply map list (take nmem signals)))
         var (map (fn [row]
                    (let [[Ex Ex2] (map #(/ % nmem)
                                        (reduce (fn [[ex ex2] v] [(+ ex v) (+ ex (* v v))])
                                                [0 0] row))]
                      (- Ex2 (* Ex Ex))))
                  (apply map (cons list mem)))
         ;; randomly draw weights from near the memories to start
         ws (doall (map #(when (nil? (fnext %))
                           (let [tmp (java.util.Arrays/copyOf (nth mem (rand-int (count mem)))
                                                              nmem)]
                             (doall
                              (map (fn [i v] (aset tmp i (+ (nth tmp i) (- (rand v) (* v 0.5)))))
                                   (range nmem)
                                   var))))
                        neurons))]
     (recur (nthnext signals (dec nmem))
            (doall (map #(if %2 [(first %1) %2] %1) neurons ws))
            mem
            nmem 0 Rfn afn [save-every nil]))
   ;; if the signals are out, we need to return the neurons
   (nil? signals) (if (not= save-every 0)
                    [memory neurons (reverse (fnext save-every))]
                    [memory neurons])
   ;; we add the next signal on and continue
   :else
   (let [k (mod (+ iter (dec nmem)) nmem)
         s1 (first signals)
         s (let [norm (Math/sqrt (reduce #(+ %1 (* %2 %2)) 0 s1))]
             (map #(/ % norm) s1))
         unscmem (map #(do (aset %1 k %2) %1) memory s)
         mem (map (fn [m] (let [norm (Math/sqrt (reduce #(+ %1 (* %2 %2)) 0 m))]
                            (doseq [j (range (count m))]
                              (aset m j (/ (nth m j) norm)))
                            m))
                  unscmem)
         R (Rfn iter)
         const (/ -0.5 (* R R))
         a (afn iter)
         [save-mod saved] save-every
         [new-neurons bmus]
         (loop [m mem, nrns (map #(assoc % 1 nil) neurons), bmus (transient [])]
           (if (nil? m)
             [nrns (persistent! bmus)]
             (let [sig (first m)
                   [bmuPos bmuW] (loop [n (next neurons)
                                        best (first neurons)
                                        sc (ED2 sig (fnext (first neurons)))]
                                   (if (nil? n)
                                     best
                                     (let [top (first n), tmp (ED2 sig (fnext top))]
                                       (if (< tmp sc)
                                         (recur (next n) top tmp)
                                         (recur (next n) best sc)))))]
               (recur (next m)
                      (doall (pmap (fn [[pos deltas] [pos0 w0]]
                                    (let [sc (* a (Math/exp (* const (ED2 pos bmuPos))))
                                          invsc (- 1.0 sc)]
                                      [pos (cons (doall (map #(+ (* %1 sc) (* %2 invsc)) bmuW w0))
                                                 deltas)]))
                                   nrns neurons))
                      (conj! bmus bmuPos)))))]
     (println [a R const (* a (Math/exp const))])
     (recur (next signals)
            (map (fn [[pos deltas]] (let [n (count deltas)]
                                      [pos (map #(/ % n) (reduce #(map + %1 %2) deltas))]))
                 new-neurons)
            mem nmem
            (inc iter)
            Rfn afn
            (if (and (not= save-mod 0) (= 0 (mod iter save-mod)))
              [save-mod (cons [bmus (doall (map fnext neurons))] saved)]
              save-every)))))

(defn simulate-retinas-kohonen
  "Yields a simulation of the retinas (using simulate-retinas), but wraps the call such that all
   analysis is done using Kohonen networks."
  [retinas hs-database & {:keys [grid-size memory image-count radius alpha save-every]
                          :or {grid-size [20 20] memory 20 image-count 500
                               radius nil alpha nil save-every 10}}]
  (let [rets (or (and (map? retinas) (list retinas))
                 (and (coll? retinas) (every? map? retinas)
                      retinas)
                 (arg-err "retinas must be a retina or a collection of retinas"))
        gsz (or (and (number? grid-size) (> grid-size 0) [grid-size grid-size])
                grid-size)
        meandim (* 0.5 (+ (first grid-size) (fnext grid-size)))
        R (or radius
              (fn [iter]
                (* meandim 0.2
                   (Math/exp (* -5.0 (/ iter image-count))))))
        a (or alpha (fn [iter] (- 1.0 (* 0.8 (/ (* iter iter) (* image-count image-count))))))
        ret1 (first rets)
        ncones (count (get ret1 :labels))
        neurons (map (fn [pos] [pos nil])
                     (:coordinates (mosaic :size grid-size)))]
    (simulate-retinas retinas hs-database :image-count image-count
                                          :analysis #(kohonen-analysis %1 neurons nil memory
                                                                       0 R a save-every)
                                          :full-reduce false)))

(defn- hebbian-analysis [signals cones opts]
  (cond
   ;; end of analysis...
   (nil? signals) cones
   ;; start of analysis...
   (nil? (get opts :initialized))
   (let [V0 0, nmem (:memory opts)]
     (recur (nthnext signals nmem)
            (apply map vector
                   (map (fn [xy M0] [(conj xy (* 0.05 (- (rand) 0.5))) V0 M0])
                        cones
                        (map #(into-array Float/TYPE %) (apply map list (take nmem signals)))))
            (assoc opts :k (dec nmem) :initialized true)))
   ;; otherwise, we'll recur eventually...
   :else
   (let [pre-s (first signals) ;; the signal we process
         max-s (apply max pre-s)
         s (map #(/ % max-s) pre-s)
         steps (get opts :steps-per-signal) ;; how many steps we run each signal
         nmem (get opts :memory) ;; the amount of memory for each cone
         k0 (get opts :k)
         k (mod k0 nmem) ;; the place we're putting this signal in the memory
         dfn (get opts :distance-fn) ;; ideal distance between cones by response angle cos
         strength (get opts :spring-strength) ;; strength of a spring
         dampen (- 1.0 (get opts :dampen))
         dt (get opts :timestep)
         dt2 (* dt dt)
         [startX startV M] cones
         step-mon (get opts :step-monitor)]
     ;; update memory...
     (doall (map #(aset %1 k (float %2)) M s))
     ;; make a similarity matrix (ideal distances)
     (let [S (let [norms (map (fn [m] (sqrt (reduce #(+ %1 (* %2 %2)) 0 m))) M)]
               (vec (map (fn [m1 n1]
                           (vec (map (fn [m2 n2]
                                       (dfn (/ (reduce + (map * m1 m2)) (* n1 n2))))
                                     M norms)))
                         M norms)))]
       (if step-mon (step-mon k0 startX startV))
       ;; recur (after running the simulation loops)
       (recur
        (next signals)
        (loop [X0 startX, V0 startV, iter 0]
          (let [;; directions and distances between cones
                D (vec (map (fn [x1]
                              (vec (map (fn [x2]
                                          (let [v (map - x2 x1)
                                                d (sqrt (reduce #(+ %1 (* %2 %2)) 0 v))]
                                            [d (if (or (= x1 x2) (= 0 d))
                                                 0
                                                 (/ (nth v 2) d))]))
                                        X0)))
                            X0))
                ;; get the accelerations
                A (doall (map (fn [drow srow]
                                (reduce + (map (fn [[d u] s] (* strength (- d s) u))
                                               drow srow)))
                              D S))
                ;; now, we update positions...
                X (doall (map (fn [x0 v a]
                                (assoc x0 2 (+ (nth x0 2) (* v dt) (* 0.5 dt2 a))))
                                        ;(vec (map #(+ %1 (* dt %2) (* dt2 0.5 %3))
                                        ;          x0 v a)))
                              X0 V0 A))
                ;; and the velocities
                V (doall (map (fn [v0 a] (+ (* dampen v0) (* dt a)))
                              V0 A))]
            ;; if we've done enough iterations, yield the result, otherwise, continue
            (if (< iter steps)
              (recur X V (inc iter))
              [X V M])))
        (assoc opts :k (inc k)))))))

(defn simulate-retinas-hebbian
  "Yields a simulation of the retinas (using simulate-retinas), but wraps the call such that all
   analysis is done using Hebbian learning implemeted by spring simulation."
  [retinas hs-database & {:keys [memory image-count dampen timestep steps-per-signal
                                 spring-strength distance-fn step-monitor]
                          :or {memory 20 image-count 1000 dampen 0.05 timestep 0.01
                               steps-per-signal 10 spring-strength 1
                               distance-fn (let [c (Math/exp -2.0)]
                                             #(+ 2.0 (* 20.0 (- (Math/exp (* -2.0 %)) c))))}}]
  (let [rets (or (and (map? retinas) (list retinas))
                 (and (coll? retinas) (every? map? retinas)
                      retinas)
                 (arg-err "retinas must be a retina or a collection of retinas"))
        ret1 (first rets)
        X (:coordinates (:mosaic ret1))
        ncones (count (get ret1 :labels))]
    (simulate-retinas
     retinas hs-database
     :image-count image-count
     :analysis #(hebbian-analysis % X {:dampen dampen
                                       :timestep timestep
                                       :memory memory
                                       :steps-per-signal steps-per-signal
                                       :spring-strength spring-strength
                                       :distance-fn distance-fn
                                       :step-monitor step-monitor
                                       :k 0})
     :full-reduce false)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Code for writing results of simulations

(defn write-simulation
  "(write-simulation output-filename retina correlation-sim-result) writes information about the
   given retina and correlation matrix produced by the simulation to the file named by
   output-filename. If the optional :embed is passed as true, then an embedding is also written."
  [output-filename retina corr-mtx & {:keys [embed]}]
  (with-open [o (java.io.DataOutputStream. (java.io.FileOutputStream. output-filename))]
    (let [labels (get retina :labels)
          n (count labels)
          coords (get (get retina :mosaic) :coordinates)
          rmap (loop [r {}, s (seq (get retina :cones)), k 0]
                 (if s
                   (recur (assoc r (first s) k) (next s) (inc k))
                   r))]
      ;; number of cones
      (.writeInt o n)
      ;; cone types
      (doseq [l labels] (.writeInt o (get rmap l)))
      ;; coordinates
      (doseq [c coords, x c] (.writeFloat o x))
      ;; correlation matrix...
      (doseq [row (range n), col (range n)]
        (.writeFloat o (corr-mtx row col)))
      ;; embedding, if requested
      (when embed
        (let [em (embed-correlation-matrix corr-mtx 3)]
          (doseq [row em, v row]
            (.writeFloat o (float v))))))))

(defn write-kohonen-simulation
  "(write-kohonen-simulation output-filename retina kohonen-sim-result) writes information about
   the given retina and resulting data produced by the kohonen simulation to the file named by
   output-filename."
  [output-filename R K]
  (let [[memory neurons traj] (if (= (count K) 3) K (conj K nil))
        dims (count (first (first neurons)))
        nmem (count (first memory))
        labels (get R :labels)
        n (count labels)
        coords (get (get R :mosaic) :coordinates)
        rmap (loop [r {}, s (seq (get R :cones)), k 0]
               (if s
                 (recur (assoc r (first s) k) (next s) (inc k))
                 r))]
    (with-open [o (java.io.DataOutputStream. (java.io.FileOutputStream. output-filename))]
      (.writeInt o n)
      (.writeInt o nmem)
      (.writeInt o (count neurons))
      (.writeInt o dims)
      ;; cone types
      (doseq [l labels] (.writeInt o (get rmap l)))
      ;; coordinates
      (doseq [c coords, x c] (.writeFloat o (float x)))
      ;; kohonen sim stuff
      (doseq [m memory, v m] (.writeFloat o (float v)))
      (doseq [n neurons]
        (let [[pos weight] n]
          (doseq [x pos] (.writeFloat o (float x)))
          (doseq [w weight] (.writeFloat o (float w)))))
      (doseq [t traj]
        (let [[bmuCoords weights] t]
          (doseq [coord bmuCoords, x coord] (.writeFloat o (float x)))
          (doseq [nrn weights, w nrn] (.writeFloat o (float w))))))))

(defn write-hebbian-simulation
  "(write-hebbian-simulation output-filename retina hebbian-sim-result) writes information about
   the given retina and result data produced by the Hebbian simulation to the file named by
   output-filename."
  [output-filename R H]
  (let [[X V M] H
        nmem (count (first M))
        labels (get R :labels)
        n (count labels)
        coords (get (get R :mosaic) :coordinates)
        rmap (loop [r {}, s (seq (get R :cones)), k 0]
               (if s
                 (recur (assoc r (first s) k) (next s) (inc k))
                 r))]
    (with-open [o (java.io.DataOutputStream. (java.io.FileOutputStream. output-filename))]
      (.writeInt o n)
      (.writeInt o nmem)
      ;; cone types
      (doseq [l labels] (.writeInt o (get rmap l)))
      ;; coordinates
      (doseq [c coords, x c] (.writeFloat o (float x)))
      ;; hebbian sim stuff
      (doseq [x X] (.writeFloat o (float (nth x 2))))
      (doseq [v V] (.writeFloat o (float v)))
      (doseq [m M, v m] (.writeFloat o (float v))))))

