;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; hyperspectral.clj
;; Hyperspectral image loaders and processors for dealing with these hyperspectral iamge databases:
;;   1) Parragas1998: http://www.cvc.uab.es/color_calibration/Bristol_Hyper/
;;   2) Chakrabarti2011: http://vision.seas.harvard.edu/hyperspec/index.html
;;   3) Foster2006: http://personalpages.manchester.ac.uk/staff/david.foster/Hyperspectral_images_of_natural_scenes_04.html
;; by Noah C. Benson

(ns brainardlab.nben.retina.hyperspectral
  (:use brainardlab.nben.retina.constants)
  (:use clojure.contrib.generic.math-functions)
  (:import java.io.File)
  (:import [java.net URL URLClassLoader]))

;; for convenience
(defn- arg-err [& txt]
  (throw (IllegalArgumentException. (apply str txt))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; These functions deal with the reading/writing of hyperspectral image caches.  These caches are
;; used for retina simulations and can be turned into an arbitrary seq of hyperspectral image
;; tensors.  They are read/written in (big-endian) binary code, so are quite fast as well.

(def ^{:dynamic true :private true} *hyperspec-add-fn*
  (fn [img]
    (throw (IllegalStateException. "No hyperspectral-cache currently being built"))))

;; actually does add a hyperspectral image to the given output stream
(defn- hyperspec-add-fn [out img]
  (let [h (count img), w (count (first img)), d (count (first (first img)))]
    (if (not (every? (fn [row] (and (= (count row) w)
                                    (every? #(= (count %) d) row)))
                     img))
      (arg-err "hyperspec-add-fn received a ragged image"))
    (.writeInt out h)
    (.writeInt out w)
    (.writeInt out d)
    (doseq [row img, col row, val col]
      (.writeFloat out val))))

(defn push-to-cache
  "Pushes the given image, absent all exclusions, to the open hyperspectral cache.  This function
   can be called only from within a build-hyperspectral-cache form.  Exclusions may be specified as
   matrices of booleans (or numbers where 0/false/nil are taken to be included), or as seqs of
   rectangles specified as [[row-ul col-ul] [row-lr col-lr]] where ul and lr indicate upper left and
   lower right, respectively (rectangles are assumed to be inclusive)."
  [hs-image]
  ;; the job here is to split the image into the largest possible squares; start by making a single
  ;; exclusion matrix
  (let [img (or (get hs-image :image-fn) (arg-err "hyperspectral image must have an :image-fn"))
        rows (or (get hs-image :rows) (arg-err "hyperspectral image must have a :rows"))
        cols (or (get hs-image :cols) (arg-err "hyperspectral image must have a :rows"))
        exclude (or (vec (map transient (get hs-image :exclusions)))
                    (arg-err "hyperspectral image must have :exclusions"))
        regions (transient [])]
    ;; okay, now we have a mask that is true for any point we should exclude
    (doseq [r (range rows), c (range cols)]
      ;; grow it as far as possible...
      (when (nil? (nth (nth exclude r) c))
        (let [[rsq csq right-next] ;; as a square first
              (loop [rr (inc r), cc (inc c)]
                ;; see if we can extend this far; right extension first, then lower
                (cond (or (= cc cols)
                          (some #(nth (nth exclude %) cc) (range r rr))) ;; excluded region
                      [(dec rr) cc false]
                      (or (= rr rows)
                          (some #(nth (nth exclude rr) %) (range c (inc cc))))
                      [rr cc true]
                      :else (recur (inc rr) (inc cc))))
              ;; now we want to extend as far as possible in whichever direction
              [rfn cfn]
              (if right-next
                (loop [cc (inc csq)]
                  (if (or (= cc cols)
                          (some #(nth (nth exclude %) cc) (range r rsq)))
                    [rsq cc]
                    (recur (inc cc))))
                (loop [rr (inc rsq)]
                  (if (or (= rr rows)
                          (some #(nth (nth exclude rr) %) (range c csq)))
                    [rr csq]
                    (recur (inc rr)))))]
          ;; okay, if this region is larger than 16 pixels, we store it and wipe the exclude region
          (if (> (* (- rfn r) (- cfn c)) 16)
            (do (doseq [rr (range r rfn), cc (range c cfn)]
                  (assoc! (nth exclude rr) cc true))
                (conj! regions [[r c] [rfn cfn]])
                (*hyperspec-add-fn* (for [row (range r rfn)]
                                      (for [col (range c cfn)]
                                        (img row col)))))))))
    (persistent! regions)))

(defn hyperspectral-image
  "Yields a hyperspectral-image structure given the hyperspectral image data specified in the
   options.  This hyperspectral-image structure can be understood by push-to-cache and
   build-hyperspectral-cache.  The arguments to this function must be as follows:
   image-fn must be a function that, given two integers, returns the pixel (a seq of numbers) at
    those coordinates
   size must be the [rows cols] of the image
   :filter may be a sequence of either masks indicating which bits should be excluded (all cells
    set to nil, false, or 0 will be excluded), rectangular regions to be excluded, specified
    by [[row-ul col-ul] [row-lr col-lr]] where ul and lr stand for upper-left and lower-right, or a
    function that takes a row and col as arguments and returns true for any pixel to be included.
   :wavelengths is by default equivalent to (range 400 721 10) and specifies the wavelengths on
    which the image is sampled."
  [image-fn [rows cols] & {:keys [filter wavelengths]}]
  (let [std-wavelengths (range 400 721 10)
        filt (or (and (ifn? filter) [filter])
                 (and (coll? filter) (= (count filter) 2)
                      (every? (fn [f] (and (= (count f) 2) (every? integer? f))) filter)
                      (cons filter nil))
                 (and (coll? filter)
                      (every? (fn [f]
                                (or (ifn? f)
                                    (and (coll? f) (= (count f) 2)
                                         (every? (fn [ff] (and (= (count ff) 2)
                                                               (every? integer? ff)))
                                                 f)))
                                    (and (= (count f) rows) (every? #(= (count %) cols) f)))
                              filter)
                      filter)
                 (and (= (count filter) rows) (every? #(= (count %) cols) filter)
                      (cons filter nil))
                 (nil? filter) nil
                 (arg-err "improperly composed :filter option"))]
    {:image-fn (cond (not (ifn? image-fn))
                     (arg-err "image-fn argument must be an IFn")
                     (or (nil? wavelengths) (= wavelengths std-wavelengths))
                     image-fn
                     :else ;; we make a new function that extrapolates
                     (let [weights
                           (map (fn [lambda]
                                  ;; start by finding the two closest wavelengths we have
                                  (let [sorted-wl (sort-by #(abs (- (first %) lambda))
                                                           (map list
                                                                wavelengths
                                                                (range (count wavelengths))))
                                        [x1 x2 x3] (map first (take 3 sorted-wl))
                                        idcs (map fnext (take 3 sorted-wl))
                                        getys (fn [spec] (map #(nth spec %) idcs))
                                        denom (* (- x1 x2) (- x1 x3) (- x2 x3))]
                                    ;; if there is this wavelength, nothing special needs to happen
                                    (if (= x1 lambda)
                                      (int (first idcs))
                                      [[x1 x2 x3] getys denom lambda])))
                                std-wavelengths)
                           ;; here's the function that interpolates
                           interp
                           (fn [spec [x1 x2 x3] getys denom x]
                             (let [[y1 y2 y3] (getys spec)
                                   a (/ (+ (* x3 (- y2 y1))
                                           (* x2 (- y1 y3))
                                           (* x1 (- y3 y2)))
                                        denom)
                                   b (/ (+ (* x3 x3 (- y1 y2))
                                           (* x1 x1 (- y2 y3))
                                           (* x2 x2 (- y3 y1)))
                                        denom)
                                   c (/ (+ (* x1 x3 (- x3 x1) y2)
                                           (* x2 x2 (- (* x3 y1) (* x1 y3)))
                                           (* x2 (- (* x1 x1 y3) (* x3 x3 y1))))
                                        denom)]
                               (+ (* a x x) (* b x) c)))]
                       (fn [r c]
                         (let [s (image-fn r c)]
                           (map (fn [w] (if (integer? w)
                                          (nth s w)
                                          (apply interp (cons s (seq w)))))
                                weights)))))
     :rows (or (and (integer? rows) (> rows 4) rows)
               (arg-err "rows must be an integer > 4"))
     :cols (or (and (integer? cols) (> cols 4) cols)
               (arg-err "cols must be an integer > 4"))
     ;; we want to mix all the exclusions into a single exclusion
     :exclusions (loop [ex (vec (repeatedly rows #(transient (vec (repeat cols nil)))))
                        q filter]
                  (if (nil? q)
                    (vec (map persistent! ex))
                    (let [f (first q)]
                      (cond (and (coll? f) (= (count f) 2)
                                 (every? #(and (= (count %) 2) (every? integer? %)) f))
                            (recur (let [[[rul cul] [rlr clr]] f]
                                     (doseq [r (range rul (inc rlr)), c (range cul (inc clr))]
                                       (assoc! (nth ex r) c true))
                                     ex)
                                   (next q))
                            (and (= (count f) rows) (every? #(= (count %) cols) f))
                            (recur (do (doseq [r (range rows), c (range cols)]
                                         (let [this-row (nth f r), px (nth this-row c)]
                                           (if (or (nil? px) (= px 0) (= px 0.0))
                                             (assoc! (nth ex r) c true))))
                                       ex)
                                   (next q))
                            (ifn? f)
                            (recur (do (doseq [r (range rows), c (range cols)]
                                         (if (not (f r c)) (assoc! (nth ex r) c true)))
                                       ex)
                                   (next q))
                            :else (arg-err "Bad exclusion to push-to-cache")))))}))
  
(defmacro build-hyperspectral-cache
  "Used to build hyperspectral image caches.  The body is executed such that any call to the
   push-to-cache function will insert the image into the cache.  After body has been executed, the
   cache is closed and the filename returned.  On error, an exception is thrown.  The file argument
   may be may be a string or a File object."
  [file & body]
  (let [flarg (gensym)
        fl (gensym)
        out (gensym)
        haddfn hyperspec-add-fn]
    `(let [~flarg ~file
           ~fl (cond (instance? java.io.File ~flarg) ~flarg
                 (string? ~flarg) (java.io.File. ~flarg)
                 :else (throw (IllegalArgumentException. "file argument must be a string or File")))
           ~out (java.io.DataOutputStream. (java.io.FileOutputStream. ~fl))]
       (try
         (binding [*hyperspec-add-fn* (fn [~'img] (~haddfn ~out ~'img))]
           ~@body)
       (finally (.close ~out))))))

(defn read-hyperspectral-cache
  "Opens the given file (either a File objerct or a filename) and reads in each hyperspectral image
   in the cache sequentially and lazilly.  Note that in order to ensure that the file is closed,
   the sequence must be entirely evaluated or the .close method must be called on the :stream field
   of the returned map.  The returned map contains three fields: :stream (the file), :index (the
   sizes of all images contained in this cache), and :images (the lazy sequence of images)."
  [file]
  (let [fl (cond (instance? java.io.File file) file
                 (string? file) (java.io.File. file)
                 :else (throw (IllegalArgumentException. "file argument must be a string or File")))
        stream (java.io.RandomAccessFile. fl "r")
        bytes-per-float 4
        bytes-per-int 4
        index (try (loop [idx nil pos 0]
                     (let [height (try (.readInt stream)
                                       (catch java.io.EOFException e
                                         (.seek stream 0)
                                         0))
                           width (if (= height 0) 0 (.readInt stream))
                           depth (if (= height 0) 0 (.readInt stream))
                           new-pos (+ pos (* 3 bytes-per-int))
                           tot (* height width depth bytes-per-float)
                           skipped (if (= tot 0) 0 (.skipBytes stream tot))]
                       (cond (= 0 skipped) (reverse idx)
                             (= skipped tot) (recur (cons [height width depth new-pos] idx)
                                                    (+ new-pos tot))
                             :else (throw (Exception. "Format error: could not read index")))))
                   (catch Exception e
                     (.close stream)
                     (throw e)))]
    (.seek stream 0)
    {:stream stream
     :index index}))

(defn images-from-cache
  "Yields a lazy seq of the hyperspectral images in the cache (read by read-hyperspectral-cache).
   An optional :delay may be set to true if delay objects are prefered (ie, if not all images will
   be used, then the unused images will not be loaded)."
  [cache & {:keys [delay]}]
  (letfn [(read-one [stream index bytes-per-float]
          (when index
            (try
              (let [load-pixels
                    (fn []
                      (let [[height width depth pos] (first index)
                            nfloats (* height width depth)
                            nbytes (* bytes-per-float nfloats)
                            buf (.asFloatBuffer
                                 (java.nio.ByteBuffer/wrap 
                                  (let [bytearr (byte-array nbytes)]
                                    (.seek stream pos)
                                    (.readFully stream bytearr 0 nbytes)
                                    bytearr)))
                            pxar (make-array Float/TYPE height width depth)]
                        (doseq [r (range height), c (range width)]
                          (.get buf (nth (nth pxar r) c) 0 depth))
                        (brainardlab.nben.retina.jvm.HSImage. pxar)))
                    pixels (if delay (clojure.core/delay (load-pixels)) (load-pixels))]
                (cons pixels
                      (lazy-seq (read-one stream (next index) bytes-per-float))))
              (catch java.io.EOFException e
                (.close stream)
                nil)
              (catch Exception e
                (.close stream)
                (throw e)))))]
    (lazy-seq (read-one (get cache :stream) (get cache :index) 4))))

(defn hyperspectral-to-RGB
  "Yields an BufferedImage object containing the RGB matrix converted from the given hyperspectral
   image.  This function uses CIE-xyz1931-matrix for conversion, which was originally obtained from
   the psychtoolbox."
  [hsimage]
  (let [[rows cols imgfn]
        (cond (instance? brainardlab.nben.retina.jvm.HSImage hsimage)
              (assoc (hsimage) 2 hsimage)
              (and (map? hsimage) (contains? hsimage :image-fn))
              [(:rows hsimage) (:cols hsimage) (:image-fn hsimage)]
              (or (coll? hsimage) (seq? hsimage))
              (let [v (or (and (vector? hsimage) hsimage) (vec hsimage))]
                [(count hsimage) (count (first hsimage)) #(nth (nth hsimage %1) %2)])
              :else (arg-err "Unrecognized image format"))
        [Rvec Gvec Bvec] CIE-xyz1931-matrix
        rgb-mx (map (fn [r]
                      (let [clrs-mx (map (fn [c]
                                           (let [hs (imgfn r c)
                                                 clr [(reduce + (map * hs Rvec))
                                                      (reduce + (map * hs Gvec))
                                                      (reduce + (map * hs Bvec))]]
                                             [clr (apply max clr)]))
                                         (range cols))]
                        [(map first clrs-mx) (apply max (map fnext clrs-mx))]))
                    (range rows))
        [rgb mx] [(map first rgb-mx) (apply max (map fnext rgb-mx))]
        bi (java.awt.image.BufferedImage. cols rows java.awt.image.BufferedImage/TYPE_INT_RGB)]
    (loop [row rgb, col (first rgb), rnum 0, cnum 0]
      (cond (= rnum rows) bi
            (= cnum cols) (recur (next row) (first (next row)) (inc rnum) 0)
            :else (let [px (first col)]
                    (.setRGB bi cnum rnum
                             (bit-or (bit-shift-left (int (* 0xff (max 0 (/ (nth px 0) mx)))) 16)
                                     (bit-shift-left (int (* 0xff (max 0 (/ (nth px 1) mx)))) 8)
                                     (bit-shift-left (int (* 0xff (max 0 (/ (nth px 2) mx)))) 0)))
                    (recur row (next col) rnum (inc cnum)))))))

(defn write-hyperspectral-as-png
  "Writes a .png file to the given file location by converting the given hyperspectral image to
   an RGB image using hyperspectral-to-RGB.  Yields the java.io.File object if successful and nil
   otherwise."
  [file hsimage]
  (let [fl (if (instance? java.io.File file) file (java.io.File. file))]
    (when (javax.imageio.ImageIO/write (hyperspectral-to-RGB hsimage) "PNG" fl)
      fl)))

(defn count-patches-in-image
  "Yields the number of patches of size patch-height x patch-width in a height x width image"
  ([height width patch-height patch-width]
     (let [rows (inc (- height patch-height))
           cols (inc (- width patch-width))]
       (if (or (<= rows 0) (<= cols 0))
         0
         (* rows cols))))
  ([image patch-height patch-width]
     (let [[height width] (if (instance? brainardlab.nben.retina.jvm.HSImage image)
                            (image)
                            [(count image) (count (first image))])
           rows (inc (- height patch-height))
           cols (inc (- width patch-width))]
       (if (or (<= rows 0) (<= cols 0))
         0
         (* rows cols)))))

(defn patch-fn-from-image
  "Yields a function that, given an index i, yields the ith image patch of size patch-height x
   patch-width from the given image.  The index i should begin at 0.  The image should be either
   an HSImage instance or a data structure that can handle nth."
  [image patch-height patch-width]
  (let [[height width] (if (instance? brainardlab.nben.retina.jvm.HSImage image)
                         (image)
                         [(count image) (count (first image))])
        rows (inc (- height patch-height))
        cols (inc (- width patch-width))]
    (if (instance? brainardlab.nben.retina.jvm.HSImage image)
      ;; an HSImage is easy to deal with...
      (fn [idx]
        (let [r (int (floor (/ idx cols)))
              c (int (floor (mod idx cols)))]
          (image r patch-height c patch-width)))
      ;; a vector or somesuch (handles nth)
      (fn [idx]
        (let [r (floor (/ idx cols))
              c (floor (mod idx cols))]
          (map (fn [row]
                 (map (fn [col]
                        (nth (nth image (+ row r)) (+ col c)))
                      (range patch-width)))
               (range patch-height)))))))

(defn count-patches-in-cache
  "Yields the total number of image caches of size patch-width by patch-height in the hyperspectral
   image cache given.  The cache must be a map, as returned by read-hyperspectral-cache."
  [cache patch-height patch-width]
  (loop [idx (or (get cache :index) (throw (IllegalArgumentException. "cache has no :index")))
         tot 0]
    (if (nil? idx)
      tot
      (recur (next idx)
             (+ tot (count-patches-in-image (first (first idx)) (fnext (first idx))
                                            patch-height patch-width))))))

(defn draw-patches
  "Yields a lazy sequence of image patches drawn from the given hyperspectral cache map, which must
   be created by read-hyperspectral-cache.  The seq will contain n image patches of size patch-width
   by patch-height sampled randomly from the database such that the first several patches come from
   the first image, the next several from the next image, etc.  The number of patches drawn from
   each is based on the fraction of possible patches that can be drawn from that image.  The
   optional argument :random-seed may be given to force a particular random seed.  If n is nil then
   all patches are returned."
  [cache n patch-height patch-width & {:keys [random-seed]}]
  (let [[total-patches-per-image total-patches]
        (loop [idx (or (get cache :index) (throw (IllegalArgumentException. "cache has no :index")))
               per-img nil
               tot 0]
          (if (nil? idx)
            [(vec (reverse per-img)) tot]
            (let [c (count-patches-in-image
                     (first (first idx)) (fnext (first idx)) patch-height patch-width)]
              (recur (next idx) (cons c per-img) (+ tot c)))))
        N (or n total-patches)
        patches-per-image (loop [s (seq total-patches-per-image), res nil, tot-drawn 0, tot-ideal 0]
                            (if (nil? s)
                              (reverse
                               (cond (< tot-drawn tot-ideal) (cons (inc (first res)) (next res))
                                     (> tot-drawn tot-ideal) (cons (dec (first res)) (next res))
                                     :else res))
                              (let [ideal (* n (/ (first s) total-patches))
                                    kideal (int (round (double ideal)))
                                    off (- tot-ideal tot-drawn)
                                    k (min (first s)
                                           (cond (> off 1) (min (inc kideal) total-patches)
                                                 (< off -1) (max (dec kideal) 0)
                                                 :else kideal))]
                                (recur (next s)
                                       (cons k res)
                                       (+ tot-drawn k)
                                       (+ tot-ideal ideal)))))
        rr (if random-seed
             (let [rndm (java.util.Random. random-seed)]
               (fn [n] (.nextInt rndm n)))
             (fn [n] (rand-int n)))]
    ;; make a lazy seq function......
    (letfn [(draw-next [img-seq count-seq total-seq count-left]
              (cond (nil? img-seq) nil
                    (= 0 count-left) (draw-next (next img-seq)
                                                (next count-seq)
                                                (next total-seq)
                                                (fnext count-seq))
                    :else (cons
                           ((first img-seq) (rr (first total-seq)))
                           (lazy-seq (draw-next img-seq count-seq total-seq (dec count-left))))))]
      (lazy-seq
       (draw-next (map #(if (= 0 %2) nil (patch-fn-from-image @%1 patch-height patch-width))
                       (images-from-cache cache :delay true)
                       patches-per-image)
                  patches-per-image
                  total-patches-per-image
                  (first patches-per-image))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Matlab functions

(def ^:private JMatIO-classname "com.jmatio.io.MatFileReader")

(let ;; here, we define the jmatio class; it can be loaded later at runtime
    [jmatio-class (ref (try (Class/forName JMatIO-classname true nil)
                            (catch ClassNotFoundException e
                              nil)))]
  (defn JMatIO-reader
    "Yields a jmatio MatFileReader (com.jmatio.io.MatFileReader) object for the given file.  This
     may throw an exception if (a) the file is not found, (b) there was a read error, ie, the file
     is not a proper Matlab v5 compatible file, or (c) the jmatio library cannot be loaded.  In the
     latter case, the library may be loaded dynamically using jmatio-load.  The file argument may be
     a (string) filename or a File object."
    [file]
    (let [reader-class @jmatio-class
          constr-file (if reader-class
                        (.getConstructor reader-class (into-array Class [java.io.File]))
                        nil)
          constr-str (if reader-class
                        (.getConstructor reader-class (into-array Class [String]))
                        nil)]
      (cond
       (nil? reader-class)
       (throw (Exception. (str "The JMatIO library has not been loaded, thus matlab .MAT files"
                               " may not currently be read.  To obtain the JMatIO library,"
                               " visit one of these pages:\n"
                               "   http://www.mathworks.com/matlabcentral/fileexchange/10759\n"
                               "   http://jmatio.sourceforge.net/\n"
                               "When the jmatio.jar file has been downloaded and unpacked, you"
                               " may use it by either re-running this java program with the"
                               " jmatio.jar file on your classpath, or by calling the"
                               " load-JMatIO function with the path of the jmatio.jar file.")))
       (string? file)
       (.newInstance constr-str (into-array String [file]))
       (instance? java.io.File file)
       (.newInstance constr-file (into-array java.io.File [file]))
       :else
       (arg-err "file must be a String or a File object"))))

  (defn load-JMatIO
    "Loads the JMatIO class and returns true if successful.  This function throws exceptions on
     error due to not wrapping a try statement around the ClassLoader and Class.forName calls.
     The path argument may be a File, String, or URL.  If the path is a string, it is assumed to be
     a file unless it begins with \"http://\".  If the JMatIO library is already loaded, this
     returns nil."
    [path]
    (when (nil? @jmatio-class)
      (dosync
       (ref-set jmatio-class
         (let [url (cond (instance? java.io.File path) (.toURL path)
                         (instance? java.net.URL path) path
                         (string? path) (if (re-find #"^http://" path)
                                          (.toURL path)
                                          (.toURL (java.io.File. path)))
                         :else (arg-err "path must be a File, URL, or string"))
               cl (URLClassLoader. (into-array java.net.URL [url]))]
           ;; just load the class!
           (.loadClass cl JMatIO-classname))))
      true))

  (defn read-matlab-file
    "Yields a map of variable name to the MLArray object saved to that name in the given MAT file.
     The given file may be either a File object or a string.  If JMatIO is not loaded, an exception
     will be thrown."
    [file]
    (clojure.lang.PersistentHashMap/create (.getContent (JMatIO-reader file)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Chakrabarti2011 Database

(defn read-Chakrabarti2011-image
  "Yields the hyperspectral image structure for the hyperspectral image specified by file. File may
   be a string or a File object.  If errors are encountered, exceptions are thrown.  All given
   filters are included in the hyperspectral image definition, along with the exclusions included
   in the Chakrabarti2011 matlab file.  If you do not wish to exclude the lbl field of the matlab
   file, then the option :filter-lbl should be false.  The contents of the calib.txt file can also
   be passed in the :calib argument; otherwise it is loaded automatically."
  [file & {:keys [filter-lbl filter calib] :or {filter-lbl true}}]
  (let [cnt (read-matlab-file file)
        cal (or (and (= 31 (count calib)) calib)
                (let [fl (if (instance? java.io.File file) file (java.io.File. file))
                      calibfl (java.io.File. (.getParent fl) "calib.txt")
                      scanner (java.util.Scanner. (slurp calibfl))]
                  (loop [res (transient [])]
                    (if (.hasNext scanner)
                      (recur (conj! res (.nextFloat scanner)))
                      (persistent! res)))))
        ref (.getArray (or (get cnt "ref") (arg-err "File does not contain a ref field")))
        lbl (.getArray (or (get cnt "lbl") (arg-err "File does not contain a lbl field")))
        fref (into-array (map float-array ref))
        [rows cols depth] (.getDimensions (get cnt "ref"))]
    (hyperspectral-image (brainardlab.nben.retina.jvm.HSImage. fref cols cal)
                         [rows cols]
                         :filter (or (and filter-lbl
                                          (or (and (coll? filter) (= (count filter) rows)
                                                   (every? #(= (count %) cols) filter)
                                                   (list filter lbl))
                                              (and (coll? filter) (= (count filter) 2)
                                                   (every? #(= (count %) 2) filter)
                                                   (list filter lbl))
                                              (and (ifn? filter) (list filter lbl))
                                              (cons lbl filter)))
                                     filter)
                         :wavelengths (range 420 721 10))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Foster2004

(defn read-Foster2004-image
  "Yields the hyperspectral image structure for the hyperspectral image specified by file. File may
   be a string or a File object.  If errors are encountered, exceptions are thrown.  All given
   filters are included in the hyperspectral image definition."
  [file & {:keys [filter]}]
  (let [dat (get (read-matlab-file file) "reflectances")
        ref (.getArray (or dat (arg-err "File does not contain a ref field")))
        [rows cols depth] (.getDimensions dat)
        fref (into-array (map float-array ref))]
    (hyperspectral-image (brainardlab.nben.retina.jvm.HSImage. fref cols)
                         [rows cols]
                         :filter filter
                         :wavelengths (range 400 (if (= depth 33) 721 711) 10))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Parragas1998

(defn read-Parragas1998-image
  "Yields the hyperspectral image structure for the hyperspectral image specified by directory. The
   directory argument may be a string or a File object.  If errors are encountered, exceptions are
   thrown.  All given filters are included in the hyperspectral image definition."
  [directory & {:keys [filter]}]
  (let [dir-file (or (and (string? directory) (java.io.File. directory))
                     (and (instance? java.io.File directory) directory)
                     (arg-err "directory must be a string or File"))
        name (or (and (.isDirectory dir-file) (.getName dir-file))
                 (arg-err "Parragas1998 images must be specified by their directory"))
        img (vec
             (map (fn [flnum]
                    (let [fl (java.io.File. (str (.getAbsolutePath dir-file) "/" name "." flnum))
                          in (or (and (.exists fl)
                                      (java.io.DataInputStream. (java.io.FileInputStream. fl)))
                                 (arg-err "directory is missing file " (.getName fl)))
                          sz (* 256 256)]
                      (try
                        (loop [to-skip 32]
                          (let [s (.skipBytes in to-skip)]
                            (cond (<= s 0) (throw (Exception. "Could not skip bytes in file!"))
                                  (< s to-skip) (recur (- to-skip s)))))
                        (loop [input (into-array Byte/TYPE (repeat sz 0))
                               off 0
                               left sz]
                          (let [s (.read in input off sz)]
                            (cond (<= s 0) (throw (Exception. "Error reading bytes from file"))
                                  (< s left) (recur input (+ off s) (- left s))
                                  :else input)))
                        (finally (.close in)))))
                  (range 400 701 10)))]
    ;; to pass the img to HSImage, we need to turn it into a row x col x depth seq
    (hyperspectral-image (brainardlab.nben.retina.jvm.HSImage.
                          (map (fn [rnum]
                                 (let [rowoff (* rnum 256)]
                                   (map (fn [col]
                                          (map #(float
                                                 (let [b (nth % (+ rowoff col))]
                                                   (if (< b 0)
                                                     (+ 256 b)
                                                     b)))
                                               img))
                                        (range 256))))
                               (range 256)))
                         [256 256]
                         :filter filter
                         :wavelengths (range 400 701 10))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Get a list of files from these sources

(def ^:dynamic *Parragas1998-file-list*
     "This variable lists all the files that were found in the Parragas1998 database (Aug 2013)"
     ["ashton2"
      "ashton2b"
      "ashton3"
      "cold1"
      "cold3"
      "fern1"
      "ferns2"
      "fort04"
      "fuschia"
      "gleaves"
      "inlab1"
      "inlab2"
      "inlab3"
      "inlab4"
      "inlab5"
      "inlab7"
      "jan10pm"
      "jan13am"
      "moss"
      "pink7"
      "plaza"
      "red1"
      "rleaves"
      "rocks"
      "rwood"
      "valley"
      "windy"
      "yellow1"
      "yleaves"])

(def ^:dynamic *Foster2004-file-list* 
     "This variable lists all the files that were found in the Foster2004 database (Aug 2013)"    
     ["/scene1/ref_crown3bb_reg1_lax.mat"
      "/scene2/ref_ruivaes1bb_reg1_lax.mat"
      "/scene3/ref_mosteiro4bb_reg1_lax.mat"
      "/scene4/ref_cyflower1bb_reg1.mat"
      "/scene5/ref_cbrufefields1bb_reg1_lax.mat"
      "/scene6/ref_braga1bb_reg1.mat"
      "/scene7/ref_ribeira1bbb_reg1.mat"
      "/scene8/ref_farme1bbbb_reg1.mat"])

(def ^:dynamic *Chakrabarti2011-file-list* 
     "This variable lists all the files that were found in the Chakrabarti2011 database (Aug 2013)"    
     ["img1.mat"
      "img2.mat"
      "imga1.mat"
      "imga2.mat"
      "imga5.mat"
      "imga6.mat"
      "imga7.mat"
      "imgb0.mat"
      "imgb1.mat"
      "imgb2.mat"
      "imgb3.mat"
      "imgb4.mat"
      "imgb5.mat"
      "imgb6.mat"
      "imgb7.mat"
      "imgb8.mat"
      "imgb9.mat"
      "imgc1.mat"
      "imgc2.mat"
      "imgc4.mat"
      "imgc5.mat"
      "imgc7.mat"
      "imgc8.mat"
      "imgc9.mat"
      "imgd2.mat"
      "imgd3.mat"
      "imgd4.mat"
      "imgd7.mat"
      "imgd8.mat"
      "imgd9.mat"
      "imge0.mat"
      "imge1.mat"
      "imge2.mat"
      "imge3.mat"
      "imge4.mat"
      "imge5.mat"
      "imge6.mat"
      "imge7.mat"
      "imgf1.mat"
      "imgf2.mat"
      "imgf3.mat"
      "imgf4.mat"
      "imgf5.mat"
      "imgf6.mat"
      "imgf7.mat"
      "imgf8.mat"
      "imgh0.mat"
      "imgh1.mat"
      "imgh2.mat"
      "imgh3.mat"])

(defn hyperspectral-file-list
  "Yields a list of files given the database locations.  For any database not given, those files are
   not included (but an empty string \"\" or a period \".\" may be used to indicate the current
   directory).  File objects are returned in a single lazy seq; only those files that exist are
   returned, unless :strings is true, in which case string objects are returned.
   Note that for the Parragas1998 database, the \"brelstaff\" directory containing the image folders
   is the expected argument.  For the Chakrabarti2011 database, CZ_hsdb directory is expected.  For
   Foster2004, the directory containing the \"scene1\", \"scene2\" etc. directories is expected."
  [& {:keys [Chakrabarti2011 Foster2004 Parragas1998 strings]}]
  ;; make sure they exist
  (when (not strings)
    (dorun
     (map #(if (and %1
                    (not (.exists (or (and (instance? java.io.File %1) %1)
                                 (java.io.File. %1)))))
             (arg-err "Hyperspectral directory for " %2 " does not exist!"))
          [Chakrabarti2011 Foster2004 Parragas1998]
          [:Chakrabarti2011 :Foster2004 :Parragas1998])))
  (let [names
        (concat
         ;; Parragas first
         (when Parragas1998
           (let [f (or (and (string? Parragas1998) (java.io.File. Parragas1998))
                       (and (instance? java.io.File Parragas1998) Parragas1998)
                       (arg-err ":Parragas1998 is neither a File nor a string"))
                 path (.getAbsolutePath f)]
             (map #(str path "/" %)
                  *Parragas1998-file-list*)))
         ;; Then Foster
         (when Foster2004
           (let [f (or (and (string? Foster2004) (java.io.File. Foster2004))
                       (and (instance? java.io.File Foster2004) Foster2004)
                       (arg-err ":Foster2004 is neither a File nor a string"))
                 path (.getAbsolutePath f)]
             (map #(str path "/" %)
                  *Foster2004-file-list*)))
         ;; Finally Chakrabarti
         (when Chakrabarti2011
           (let [f (or (and (string? Chakrabarti2011) (java.io.File. Chakrabarti2011))
                       (and (instance? java.io.File Chakrabarti2011) Chakrabarti2011)
                       (arg-err ":Chakrabarti2011 is neither a File nor a string"))
                 path (.getAbsolutePath f)]
             (map #(str path "/" %)
                  *Chakrabarti2011-file-list*))))]
    (if strings
      names
      (filter #(.exists %) (map #(java.io.File. %) names)))))
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Handy function for auto-building the database
(defn autobuild-hyperspectral-cache
  "See build-hyperspectral-cache.  Give then locations of the hyperspectral image databases, builds
   a hyperspectral image cache from all images."
  [file & {:keys [Chakrabarti2011 Foster2004 Parragas1998]}]
  (build-hyperspectral-cache file
    ;; first Chakrabarti2011
    (when Chakrabarti2011
      (doall (map #(push-to-cache (read-Chakrabarti2011-image %))
                  (hyperspectral-file-list :Chakrabarti2011 Chakrabarti2011))))
    ;; Then Foster2004
    (when Foster2004
      (doall (map #(push-to-cache (read-Foster2004-image %))
                  (hyperspectral-file-list :Foster2004 Foster2004))))
    ;; then Parrayas1998
    (when Parragas1998
      (doall (map #(push-to-cache (read-Parragas1998-image %))
                  (hyperspectral-file-list :Parragas1998 Parragas1998))))))
