;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; core.clj
;; Code for simulating a human retina as it responds to natural images.
;; by Noah C. Benson

(ns brainardlab.nben.retina.core
  (:use brainardlab.nben.retina.constants)
  (:use clojure.contrib.generic.math-functions))

;; Retinas are basically just maps; they store all the necessary information for a retina.  All of
;; the private functions in this file assume that sanity checking has been done thus run as fast
;; as possible; the public methods sanity check input.

(defn- arg-err [& txt] (throw (IllegalArgumentException. (apply str txt))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Mosaic Code

;; The mosaic of a retina is a a sequence of 2-vectors
(defn mosaic?
  "Yields true if the give object is a mosaic and false otherwise."
  {:version "1.0"
   :added "1.0"
   :author "Noah C. Benson"}
  [m]
  (and (map? m)
       (coll? (get m :coordinates))
       (every? #(and (seq %)
                     (= (count %) 2)
                     (number? (first %))
                     (number? (fnext %)))
               (get m :coordinates))
       (fn? (get m :signal-fn))))

;; used by both signal-fn and surround-fn for blur and surround suppression
(defn- normal-2D [[x0 y0] [x y] sigma]
  (let [dx (- x0 x), dy (- y0 y)]
    (/ (Math/exp (/ (+ (* dx dx) (* dy dy))
                    (* -2.0 (* sigma sigma))))
       (Math/sqrt (double (* Math/PI 2.0 sigma))))))

;; Used by mosaic and retina to create signal functions for mosaics
(defn- signal-fn [xy size blur eps]
  (let [[height width] (or size [(ceil (- (max (map fnext xy)) eps))
                                 (ceil (- (max (map first xy)) eps))])
        deps (if blur
               (map (fn [x0-and-y0]
                      (filter #(> (nth % 2) eps)
                              (for [row (range height), col (range width)]
                                [row col (normal-2D x0-and-y0 [col row] blur)])))
                    xy)
               (map (fn [[x y]]
                    (let [xbot (floor (+ x eps)), xtop (ceil (- x eps)),
                          ybot (floor (+ y eps)), ytop (ceil (- y eps))]
                      (filter
                       #(> (nth % 2) 0)
                       (map (fn [[col row]]
                              (if (or (> col width) (> row height)
                                      (< col 1) (< row 1))
                                (arg-err "signal-fn cannot interpolate outside of [1,w] and [1,h]")
                                (let [d (+ (* (- x col) (- x col))
                                           (* (- y row) (- y row)))]
                                  [row col (if (> d 1) 0.0 (- 1.0 d))])))
                            (or (and (= xbot xtop)
                                     (or (and (= ybot ytop)
                                              (list [xbot ybot]))
                                         (list [xbot ybot] [xbot ytop])))
                                (and (= ybot ytop)
                                     (list [xbot ybot] [xtop ybot]))
                                (list [xbot ybot] [xbot ytop] [xtop ybot] [xtop ytop]))))))
                    xy))
        weight-total (map (fn [dep] (reduce #(+ %1 (nth %2 2)) 0 dep)) deps)
        scaled-deps (map (fn [dep wtot]
                           (map (fn [[r c w]] [r c (/ w wtot)]) dep))
                         deps weight-total)]
    (fn [image]
      (map (fn [dep]
             (loop [s dep, w 0.0, sig nil]
               (if s
                 (let [[row col weight] (first s)
                       spec (nth (nth image (dec row)) (dec col))]
                   (recur (next s)
                          (+ w weight)
                          (if sig (map + sig spec) spec)))
                 (vec (map #(/ % w) sig)))))
           scaled-deps))))

(defn- surround-fn [xy surround eps]
  (if surround
    (let [idx (range (count xy))
          [weight std] (cond (and (coll? surround) (= (count surround) 2)) surround
                             (= :automatic surround) [0.25 3]
                             (number? surround) [0.25 surround]
                             :else (arg-err "Invalid :surround argument"))
          surr (map (fn [x0-and-y0]
                      (let [comps (filter #(> (first %) eps)
                                          (map (fn [x-and-y cone-index]
                                                 [(if (= x-and-y x0-and-y0)
                                                    0
                                                    (normal-2D x0-and-y0 x-and-y std))
                                                  cone-index])
                                               xy
                                               idx))
                            scaling-factor (/ weight (reduce #(+ %1 (first %2)) 0 comps))]
                        (map (fn [[cone-weight idx]] [(* cone-weight scaling-factor) idx])
                             comps)))
                  xy)]
      (fn [sig]
        (map (fn [center-cone-surround center-cone-index]
               (reduce -
                       (nth sig center-cone-index)
                       (map (fn [[cone-weight cone-index]]
                              (* cone-weight (nth sig cone-index)))
                            center-cone-surround)))
             surr
             idx)))
    identity))

(defn mosaic
  "Yields a mosaic given the parameters provided.  Mosaics are stored as vectors or 2-element
   vectors, each of which specifies a cone position.  Note that when showing a retina images, the
   images' coordinates are considered to lie at (1,1) for the top left corner to (width,height) for
   the lower right corner.  Because the optics of the eye invert the image, it is sensible to
   consider a retina covering the same range of coordinates to exist on a standard axis with the y
   unit vector pointing up.  The mosaic structure that is returned is a map with three fields:
   :coordinates for the cone positions and :signal-fn for the function that, given a hyperspectral
   image tensor, returns the hyper-spectrum observed by each cone in the same order as they are
   listed in the mosaic.

   The following optional arguments may be used:
    :layout may be :rectilinear (default) or :hexagonal
    :filter may be :rectangle (default), :circular, or a function that will filter over the mosaic
       parameters
    :size may be a positive integer (default: 20) or 2-element seq of positive integers representing
       the height and width of the retina; the height and width must match the image patches shown
       to this retina.
    :separation (default 1.0) specifies how far a cone is from its nearest neighbors.
    :x0 (default: 1.0) specifies the starting x-value used in tiling the space.
    :y0 (default: 1.0) specifies the starting y-value used in tiling the space.
    :blur (default: nil) specifies the standard deviation of a blurring Gaussian, if any, to use on
       input images when calculating the response signal."
  {:version "1.0"
   :added "1.0"
   :author "Noah C. Benson"}
  [&{:keys [layout filter size separation x0 y0 pixel-epsilon blur]
     :or {layout :rectilinear filter nil size 20 separation 1.0
          x0 1.0 y0 1.0 pixel-epsilon 0.001}}]
  (let [sep (or (and (number? separation)
                     (> separation 0)
                     separation)
                (arg-err ":separation must be a number > 0"))
        eps (or (and (number? pixel-epsilon)
                     (>= pixel-epsilon 0)
                     (< pixel-epsilon 0.5)
                     pixel-epsilon)
                (arg-err ":pixel-epsilon must be a number >= 0 and < 0.5"))
        [height width] (or (and (integer? size)
                                (> size 0)
                                [size size])
                           (and (coll? size) (= (count size) 2)
                                (integer? (first size)) (integer? (fnext size))
                                (> (first size) 0) (> (fnext size) 0)
                                size)
                           (arg-err ":size must be a number or a pair of numbers > 0"))
        [cx cy] [(* 0.5 (inc width)) (* 0.5 (inc height))]
        [a b] [(* 0.5 width) (* 0.5 height)]
        xstart (or (and (number? x0) x0) (arg-err ":x0 must be a number"))
        ystart (or (and (number? y0) y0) (arg-err ":x0 must be a number"))
        filt (or (and (or (= filter :circle) (= filter :ellipse))
                      (fn [[x y]] (let [dx (/ (- x cx) a), dy (/ (- y cy) b)]
                                    (<= (+ (* dx dx) (* dy dy)) 1.0))))
                 (and (or (= filter :rect) (= filter :rectangle))
                      (fn [[x y]] (let [dx (/ (- x cx) a), dy (/ (- y cy) b)]
                                    (and (<= (* dx dx) 1.0)
                                         (<= (* dy dy) 1.0)))))
                 (and (ifn? filter) filter)
                 (and filter (arg-err ":filter option must be a function, nil, :rect, or :circle")))
        xy (case layout
             :rectilinear ;; easy case; place them every unit like in an image
             (loop [res (transient []) x xstart y ystart]
               (cond (> x width) (recur res xstart (+ y sep))
                     (> y height) (persistent! res)
                     :else (recur (if (and (>= x 1.0) (>= y 1.0)
                                           (or (nil? filt) (filt x y)))
                                    (conj! res [x y])
                                    res)
                                  (+ x sep)
                                  y)))
             :hexagonal ;; hard case; there's an offset for off rows
             (let [half-sep (* 0.5 sep) dsep (* 0.5 (sqrt 3.0))]
               (loop [res (transient []) x xstart y ystart xoffs false]
                 (cond (> y height) (persistent! res)
                       (> x width) (recur res
                                          (if xoffs xstart (+ xstart half-sep))
                                          (+ ystart dsep)
                                          (not xoffs))
                       :else (recur (if (and (>= x 1.0) (>= y 1.0)
                                             (or (nil? filt) (filt x y)))
                                      (conj! res [x y])
                                      res)
                                    (+ x dsep) y xoffs))))
             (arg-err ":layout must be :rectilinear or :hexagonal"))]
    {:params {:layout layout :filter filter :size size
              :separation separation :x0 x0 :y0 y0 :pixel-epsilon pixel-epsilon}
     :coordinates xy
     :size [height width]
     :signal-fn (signal-fn xy [height width] blur eps)}))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Spectra Code

(defn spectral-sensitivities
  "Yields a map of cone id to spectral sensitivities for the given cones and given lambda-max values."
  {:version "1.0"
   :added "1.0"
   :author "Noah C. Benson"}
  [cones lambda-max]
  (if (not (or (seq? cones) (coll? cones)))
    (arg-err "cones, passed to spectra, must be seqable"))
  (if (not (or (seq? lambda-max) (coll? lambda-max)))
    (arg-err "lambda-max, passed to spectra, must be seqable"))
  (if (not (or (map? lambda-max) (= (count cones) (count lambda-max))))
    (arg-err "cones and lambda-max should be the same length"))
  (loop [qc (seq cones), ql (seq lambda-max), res {}]
    (if (nil? qc)
      res
      (let [c (first qc)
            lm (or (and (map? lambda-max)
                        (or (contains? lambda-max c)
                            (arg-err "lambda-max, when a map, must contain all cone ids"))
                        (get lambda-max c))
                   (first ql))
            up (int (ceil lm))
            down (int (floor lm))
            wup (* 0.5 (- 1.0 (- up lm)))
            wdown (* 0.5 (- 1.0 (- lm down)))
            sup (get spectral-sensitivities-matrix up)
            sdown (get spectral-sensitivities-matrix down)]
        (if (or (nil? sup) (nil? sdown))
          (arg-err "Spectral sensitivities can only be created for values in [420 560]"))
        (recur (next qc)
               (next ql)
               (assoc res c (vec (map #(+ (* wup %1) (* wdown %2)) sup sdown))))))))

(defn cone-class-responses-to-spectra
  "Yields the given cone type's response to the given spectra of light from the given retina.  The
   option :wavelengths may be given if the wavelengths at which the spectra are sampled is not
   equivalent to (range 400 721 10)."
  [retina cone-class spectra & {:keys [wavelengths]}]
  (let [senss (or (get retina :spectral-sensitivities)
                  (arg-err "Retina must contain spectral sensitivities"))
        orig-sens (or (get senss cone-class)
                      (arg-err "Retina must contain given cone"))
        sens (if (nil? wavelengths)
               orig-sens
               (map #(let [rescaled (/ (- % 400) 10.0)
                           lower (or (and (number? %) (int (floor rescaled)))
                                     (arg-err "wavelengths must be numbers"))
                           upper (or (and (>= % 400) (<= % 720) (int (ceil rescaled)))
                                     (arg-err "wavelengths must be in [400,720]"))
                           lwght (- 1.0 (- rescaled lower))
                           uwght (or (and (= upper lower) 0)
                                     (- 1.0 (- upper rescaled)))]
                       (+ (* lwght (nth orig-sens lower))
                          (* uwght (nth orig-sens upper))))
                    wavelengths))]
    (map (fn [spectrum] (reduce + (map * sens spectrum)))
         spectra)))

(defn cone-responses
  "Yields the responses of all cones in the retina (in the order listed in both the :coordinates and
   the :labels of the retina). The image argument must be a hyperspectral image matrix (ie, a seq of
   seqs of seqs or 3D vector) that is being exposed to the retina or a map containing entries for
   :rows, :cols, and :image-fn, which must give the size of the retina and provide access to the
   hyperspectral seq of numbers when a row and column is passed to the :image-fn.  If the
   wavelengths that the image is sampled on are not equivalent to (range 400 720 10), then the
   :wavelengths argument should be passed to indicate the wavelengths at which the image was
   sampled."
  [retina image & {:keys [wavelengths]}]
  (let [mosaic (or (get retina :mosaic) (arg-err "retina must contain a :mosaic"))
        signalfn (or (get mosaic :signal-fn) (arg-err "mosaic must contain a :signal-fn"))
        cone-classes (or (get retina :cones) (arg-err "retina must have a :cones"))
        labels (or (get retina :labels) (arg-err "retina must have :labels"))
        img (cond (instance? brainardlab.nben.retina.jvm.HSImage image) (.seq image)
                  (map? image) (let [rows (get image :rows)
                                     cols (get image :cols)
                                     imgfn (get image :image-fn)]
                                 (map (fn [r] (map #(imgfn r %) (range cols))) (range rows)))
                  (coll? image) image
                  :else (arg-err "image must be a map, a seq, or an HSImage"))
        signals (vec (signalfn img))
        ;; build up the separate cone class responses
        signals-by-class
        (loop [c labels
               s signals
               res (reduce #(assoc %1 %2 (transient [])) {} cone-classes)]
          (if (nil? c)
            (reduce #(assoc %1 %2 (persistent! (get %1 %2))) res cone-classes)
            (let [cone-class (first c)
                  signal (first s)]
              (conj! (get res cone-class) signal)
              (recur (next c) (next s) res))))
        ;; next, get the responses
        responses-by-class
        (reduce #(assoc %1
                   (key %2) (cone-class-responses-to-spectra
                             retina (key %2) (val %2) :wavelengths wavelengths))
                {}
                signals-by-class)
        responses
        ;; finally, zip them back into the correct vector
        (loop [q labels
               responses responses-by-class
               r (transient [])]
          (if (nil? q)
            (persistent! r)
            (let [c (first q)
                  response-queue (get responses c)]
              (recur (next q)
                     (assoc responses c (next response-queue))
                     (conj! r (first response-queue))))))
        surrfn (get retina :surround-fn)]
    ;; if there's a surround, deal with it
    (if surrfn (surrfn responses) responses)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Cone Label Code

(defn label-cones
  "Yields a sequence of labels for the cones in the given mosaic using the makeup instructions.  The
   cones argument must be a sequable collection of cone id's; retinal-mosaic must be a mosaic map
   (ie, (mosaic? retinal-mosaic) must be true), and makeup must be a map of makeup instructions for
   each cone id.  These maps may contain the following instructions:
   If a cone's id maps to ..., then ...
    nil -> by default, when all other cone types have been placed, the remaining cones in the mosaic
      will be assigned random types from all cone types whose cone id's mapped to nil or were not
      present in the instructions map.
    <number> -> If the cone id maps to a number, this is taken to be a relative ratio part; ie, if
      {:L 5, :M 6, :S 1} then the ratio of :L to :M to :S cones in the final labelling should be
      close to 5:6:1.
    {:even-spacing <boolean>} -> any instruction map containing :spacing will be evenly spaced 
      according to a heuristic; if <spacing> is :hexagonal or :rectilinear then the spacing is such
      that the distance between these cones is approximately hexagonal or rectilinear.
    {:noise <noise>} -> for cones with :even-spacing, this specifies the fraction of the distance 
      between cones to allow for noise in the placement; <noise> is the std.dev. of a Gaussian
      centered at the ideal spacing position where 1 is scaled to be the distance between spaced
      cones (default is 0).
    {:weight <number>} -> same as just <number>, but other options may be provided.
    {:fraction <number>} -> exactly (round (* <number> (count retinal-mosaic))) cones of this type
      should be placed.
    {:count <number>} -> exactly <number> cones of this type should be placed."
  {:version "1.0"
   :added "1.0"
   :author "Noah C. Benson"}
  [cones retinal-mosaic instructions & {:keys [seed] :or {seed nil}}]
  ;; First, cleanup the makeup instructions; save the evenly spaced ones for first use and those
  ;; without instructions as remainder-fillers
  (let [inst (or (and (map? instructions) instructions)
                 (and (coll? instructions) (= (count instructions) (count cones))
                      (apply hash-map (apply concat (map list cones instructions))))
                 (arg-err "instructions arg must be a map or a coll the same size as cones"))
        xy (and (or (mosaic? retinal-mosaic) (arg-err "retinal-mosaic must be a mosaic"))
                (get retinal-mosaic :coordinates))
        ;; for convenience...
        num-cones (count xy)
        num-cone-types (count cones)
        [height width] (get retinal-mosaic :size)
        ;; sort cone instructions into types: evenly-spaced cones as well as cones whose count is
        ;; specified exactly, by weight, or not specified/random
        [spaced weighted exact remainder]
        (loop [q (seq cones), sp {}, wt {}, ex {}, rm nil]
          (if (nil? q)
            [(if (empty? sp) nil sp), (if (empty? wt) nil wt), (if (empty? ex) nil ex), rm]
            (let [cone-id (first q)
                  entry (get inst cone-id)
                  entry-spacing (get entry :even-spacing)
                  new-sp (or (and (not entry-spacing) sp)
                             (and (true? entry-spacing)
                                  (assoc sp cone-id entry))
                             (arg-err ":even-spacing must be true or false/nil"))
                  entry-weight (get entry :weight)
                  new-wt (or (and (number? entry) (or (> entry 0) (arg-err "weights must be > 0"))
                                  (assoc wt cone-id entry))
                             (and (nil? entry-weight) wt)
                             (and (number? entry-weight) (> 0 entry-weight)
                                  (assoc wt cone-id entry-weight))
                             (arg-err ":weight must be a number > 0"))
                  entry-frac (get entry :fraction)
                  entry-count (get entry :count)
                  new-ex (or (and entry-frac entry-weight
                                  (arg-err "cannot specify both :weight and :fraction"))
                             (and entry-count entry-weight
                                  (arg-err "cannot specify both :weight and :count"))
                             (and entry-count entry-frac
                                  (arg-err "cannot specify both :fraction and :count"))
                             (and entry-frac
                                  (or (and (number? entry-frac) (> entry-frac 0) (<= entry-frac 1)
                                           (assoc ex cone-id (round (* entry-frac (count xy)))))
                                      (arg-err ":fraction must be a number > 0 and <= 1")))
                             (and entry-count
                                  (or (and (integer? entry-count) (> entry-count 0)
                                           (assoc ex cone-id entry-count))
                                      (arg-err ":count must be an integer > 0")))
                             ex)
                  new-rm (or (and (= wt new-wt) (= ex new-ex)
                                  (cons cone-id rm))
                             rm)]
              (recur (next q) new-sp new-wt new-ex new-rm))))
        ;; get exact counts of each type now, excluding the remainder
        count-exacts
        (or (and weighted remainder
                 (arg-err ":weight and :remainder instructions cannot be used together"))
            (let [count-ex (reduce + (vals exact))]
              (if (> (+ count-ex (count weighted) (count remainder)) num-cones)
                (arg-err "The number of exactly specified cone labels is so high that other cone"
                         " types cannot be included")
                count-ex)))
        ;; get the remaining counts
        count-inexacts (- (count xy) count-exacts)
        ;; get the total of all weights
        weight-total (if weighted (reduce + (vals weighted)) (count remainder))
        ;; we start by making an ordered list (from least to greatest weight) of cone -> weight
        ;; where the weight is approximately the number of cones; then we solidify the number of
        ;; cones of each cone type with a minimum of 1 cone
        weighted-counts
        (loop [q (sort-by fnext
                          (if weighted
                            (map #(list (key %) (/ (* (val %) count-inexacts) weight-total))
                                 weighted)
                            (map #(list % (/ count-inexacts weight-total))
                                 remainder)))
               remaining count-inexacts
               res {}]
          (if (nil? q)
            (if (= remaining 0)
              res
              (recur nil
                     (dec remaining)
                     (let [i (rand-int (count res))
                           me (nth (seq res) i)]
                       (assoc res (key me) (inc (val me))))))
            (let [[cid w] (first q)
                  cnt (min (- remaining (dec (count q))) (max 1 (int (floor (double w)))))]
              (recur (next q) (- remaining cnt) (assoc res cid cnt)))))
        ;; and get exact counts for each cone id
        counts (merge exact weighted-counts)
        ;; this is the label assignment that we will modify
        labels (transient (vec (repeat num-cones nil)))
        ;; first, place the spaced cones; save the cones that are not yet assigned
        unassigned
        (loop [q (seq spaced)
               initlist (map list xy (range (count xy)))]
          (if (nil? q)
            (map fnext initlist)
            (let [[cone-id sp] (first q)
                  n (get counts cone-id) ;; how many of these cones?
                  area (* height width)
                  radius2 (- (ceil (/ area n)) 1)
                  ;; now, find the positions for these cones
                  labeled-cones
                  (loop [i 0, left initlist, r2 radius2, res nil, attempt 0]
                    (cond (< r2 1) (throw (Exception. "Could not draw evenly spaced cones"))
                          (> attempt 10) (recur 0 initlist (- r2 1) nil 0)
                          (= i n) res
                          (empty? left) (recur 0 initlist r2 nil (inc attempt))
                          :else (let [k (rand-int (count left))
                                      [x0 y0] (first (nth left k))
                                      idx (fnext (nth left k))]
                                  (recur (inc i)
                                         (filter #(let [[x y] (first %)]
                                                    (> (+ (* (- x x0) (- x x0))
                                                          (* (- y y0) (- y y0)))
                                                       r2))
                                                 left)
                                         r2
                                         (cons (fnext (nth left k)) res)
                                         attempt))))]
              ;; set all the labels in the transient...
              (loop [s (seq labeled-cones)]
                (when s
                  (assoc! labels (first s) cone-id)
                  (recur (next s))))
              ;; and recur...
              (recur (next q) (filter #(nil? (nth labels (fnext %))) initlist)))))]
    ;; last, place the remaining cones; these can be placed randomly
    ;; we randomize by drawing between 0 and K where K is the total number of cones to be drawn yet;
    ;; the first cone is assigned to class Ci such that the random draw is between sum(Kj for j < i)
    ;; and sum(Kj for j <= i)
    (loop [Ktotal count-inexacts
           K (loop [q spaced, res counts]
               (if q
                 (recur (next q) (dissoc res (first (first q))))
                 res))
           left unassigned]
      (when left
        ;; we assign the first cone of those left
        (let [r (rand-int Ktotal)
              [lbl new-count] (loop [q K, c 0]
                                (if (nil? q)
                                  (throw (Exception.
                                          (str "Could not assign a cone class!  " c " :: " K " :: " r " :: " Ktotal)))
                                  (let [[cid n] (first q)
                                        top (+ c n)]
                                    (if (and (>= r c) (< r top))
                                      [cid (dec n)]
                                      (recur (next q) top)))))]
          (assoc! labels (first left) lbl)
          (recur (dec Ktotal) (assoc K lbl new-count) (next left)))))
    ;; finally, return the labels!
    (persistent! labels)))

(def ^:private human-cone-makeup {:S {:even-spacing true, :fraction 0.06}})

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Retina Code

(defn retina?
  "Yields true if the given object is a retina, otherwise false."
  {:version "1.0"
   :added "1.0"
   :author "Noah C. Benson"}
  [r]
  (and (map? r) (true? (get r :retina?))))

;; These private maps are used below to interpret options to the retina function
(def human-lambda-max {:L 558.9 :M 530.3 :S 420.7})

(let [mosaicfn mosaic]
(defn retina
  "Yields a retina object using the given arguments.  Optional arguments include:
   :mosaic (default :automatic), if :automatic, this creates a rectilinear mosaic with a separation
     of 1 to fill the retina; this may alternately be a formed mosaic or a seq of mosaic
     coordinates, in which case a mosaic structure is created out of the points.
   :cones (default :trichromat), indicates the cones in the mosaic; this may be either a seq of
     cone type identifiers (e.g. [:L :M :S]) or one of :trichromat, :tetrachromat, or :dichromat.
     These map to [:L :M :S], [:L :A :M :S], and [:L :S] respectively.  As long as one of these
     tags or seq's is passed, the below arguments will auto-adjust their defaults.
   :lambda-max (default :automatic), indicates the lambda-max values of the cones to be used in
     this retina.  This may be either a map of cone id to lambda-max value or a seq in the same
     order as the :cones option.  By default uses {:L 558.9, :A 545.0, :M 530.0, :S 420.7}.
   :cone-makeup (default :automatic), indicates the parameters of the random draw for the cones.
     This may be a map of cone id to instruction or a seq in the same order as the :cones option.
     Instructions may be a number (taken as a ratio relative to other ratios) or a map which may
     contain the fields :exact (indicating, if true, that exactly this many cones must be drawn,
     within the limitations placed by using an integer number of cones), :count (indicating that
     a specific number of a cone type should be drawn), :fraction (indicating that this fraction
     of the cones should be this cone type), or :ratio (indicating that the ratio of this cone
     type to that of others should be drawn.
   :surround (default :automatic), declares the type of surround suppression; nil indicates none,
     otherwise may be a standard deviation (in which case 0.25 is the assumed weigh) or [weight
     standard-deviation].
   :type (default :trichromat), indicates which kind of default values should be used for the 
     above options.  A value of nil uses the default options; other allowed values are
     :tetrachromat and :dichromat.  Additionally, :Bayer will produce a Bayer-like rectilinear
     mosaic."
  {:version "1.0"
   :added "1.0"
   :author "Noah C. Benson"}
  [& {:keys [mosaic cones lambda-max cone-makeup type pixel-epsilon surround]
      :or {mosaic :rectilinear, cones :automatic, lambda-max :automatic,
           cone-makeup :automatic, type :trichromat, pixel-epsilon 0.001, surround :automatic}}]
  (let ;; first we parse the params into something sanity-checked; for this we make a temp
      [eps (or (and (number? pixel-epsilon)
                    (>= pixel-epsilon 0)
                    (< pixel-epsilon 0.5)
                    pixel-epsilon)
               (arg-err ":pixel-epsilon must be a number >= 0 and < 0.5"))
       base (case type :trichromat {:cones [:L :M :S]
                                    :lambda-max human-lambda-max
                                    :cone-makeup human-cone-makeup
                                    :mosaic :rectilinear}
                       :dichromat {:cones [:L :S]
                                   :lambda-max (dissoc human-lambda-max :M)
                                   :cone-makeup human-cone-makeup
                                   :mosaic :rectilinear}
                       :tetrachromat {:cones [:L :A :M :S]
                                      :lambda-max (assoc human-lambda-max :A 545.0)
                                      :cone-makeup human-cone-makeup
                                      :mosaic :rectilinear}
                       (arg-err ":type must be :trichromat, :dichromat, or :tetrachromat"))
       M ;; start by making the mosaic; if they gave us one, we need to check it against size
       (or (and (mosaic? mosaic) mosaic)
           (and (or (= mosaic :rectilinear)
                    (= mosaic :hexagonal))
                (mosaicfn :size 20 :layout mosaic))
           (and (coll? mosaic)
                (or (and (every? #(and (= (count %) 2) (number? (first %)) (number? (fnext %)))
                                 mosaic)
                         {:params nil :coordinates mosaic :signal-fn (signal-fn mosaic nil eps)})
                    (and (or (= (first mosaic) :rectilinear)
                             (= (first mosaic) :hexagonal))
                         (case (count mosaic)
                           1 (mosaicfn :size 20 :layout :rectilinear)
                           2 (and (integer? (fnext mosaic)) (> 0 (fnext mosaic))
                                  (mosaicfn :size (fnext mosaic) :layout (first mosaic)))
                           3 (and (integer? (fnext mosaic)) (> 0 (fnext mosaic))
                                  (integer? (nth mosaic 2)) (> 0 (nth mosaic 2))
                                  (mosaicfn :size (next mosaic) :layout (first mosaic)))
                           (arg-err
                            ":mosaic [<type> ...] must be followed by 1 or 2 positive integers")))))
           (arg-err ":mosaic must be a valid mosaic, :rectilinear, or :hexagonal"))
       C (or (and (= cones :automatic) (:cones base))
             (and (coll? cones) (vec cones))
             (arg-err ":cones must be :automatic or a seq of cone labels"))
       Lmax (or (and (= lambda-max :automatic) (:lambda-max base))
                (and (map? lambda-max)
                     (every? #(let [v (get lambda-max %)]
                                (and (number? v) (> v 0)))
                             C)
                     lambda-max)
                (and (coll? lambda-max)
                     (= (count lambda-max) (count C))
                     (every? #(and (number? %) (> % 0)) lambda-max)
                     (apply hash-map (apply concat (map list C lambda-max))))
                (arg-err ":lambda-max must be a map or seq of numbers >0 the same size as :cones"))
       makeup (or (and (= cone-makeup :automatic) (:cone-makeup base))
                  (and (= cone-makeup :human) human-cone-makeup)
                  (and (map? cone-makeup) cone-makeup)
                  (arg-err ":cone-makeup must be :automatic, :human, or a map"))
       lbls (label-cones C M makeup)]
    {:cones C
     :mosaic M
     :spectral-sensitivities (spectral-sensitivities C Lmax) ;; make spectra
     :labels lbls ;; make labels
     :surround-fn (surround-fn (:coordinates M) surround eps)
     :params {:type type
              :lambda-max lambda-max
              :pixel-epsilon pixel-epsilon
              :cone-makeup cone-makeup
              :cones cones
              :surround surround}})))
