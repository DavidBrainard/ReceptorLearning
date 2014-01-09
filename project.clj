(defproject retina "0.1.0-SNAPSHOT"
  :description (str "The retina library enables high quality simulation of the responses of retinal"
                    " mosaics to natural images.")
  :url "http://github.com/nben/retina"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  ;; We want to be able to use large amounts of memory
  ;; This is something that you'll likely want to set for each machine on which this library is used
  :jvm-opts ["-Xmx7G"]
  ;; dependencies
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [org.clojure/clojure-contrib "1.2.0"]
                 [net.sourceforge.jmatio/jmatio "1.0"]
                 [incanter "1.5.4"]]
  :plugins [[codox "0.6.6"]]
  ;; for generating codox documentation
  :codox {:output-dir "doc/html"}
  :resource-paths ["resources/mdsj.jar"]
  ;; location of source codes
  :source-paths ["src/main/clojure"]
  :java-source-paths ["src/main/java"]
  :test-paths ["src/test/clojure"]
  ;; the main class
  :main brainardlab.nben.retina.main
  ;; target for aot compilation
  :target-path "target/"
  :compile-path "target/classes"
  ;; and the namespaces to aot compile
  :aot [brainardlab.nben.retina.constants
        brainardlab.nben.retina.core
        brainardlab.nben.retina.hyperspectral
        brainardlab.nben.retina.simulation
        brainardlab.nben.retina.util
        brainardlab.nben.retina.main]
  ;; targets that get cleaned...
  :clean-targets [:target-path :compile-path]
  ;; jar file options...
  :jar-name "retina.jar"
  :omit-source false
  :jar-exclusions [#"(?:^|/).svn/"]
  ;; And some options for the REPL...
  :repl-options {:init (use '(brainardlab.nben.retina core hyperspectral simulation util))})
