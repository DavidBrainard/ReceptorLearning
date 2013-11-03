////////////////////////////////////////////////////////////////////////////////////////////////////
// KohonenAnalysis.java
// Definition of the class that manages Kohonen network analysis
// by Noah C. Benson

package brainardlab.nben.retina.jvm;
import clojure.lang.*;
import java.util.*;

public final class KohonenAnalysis
   extends clojure.lang.AFn
{
   // private data: this stores the neuron positions, their weights, etc.
   private class Neuron {
      public int index;
      public float[] weights;
      public float[] oldWeights;
      public float[] position;
      public Neuron[] neighbors;
      public float[] distances;
      public Neuron(int idx, float[] pos, int nweights, Random r) {
         index = idx;
         position = pos.clone();
         weights = new float[nweights];
         float tot = 0;
         int i;
         for (i = 0; i < nweights; ++i) {
            weights[i] = r.nextFloat();
            tot += weights[i]*weights[i];
         }
         tot = (float)Math.sqrt(tot);
         for (i = 0; i < nweights; ++i)
            weights[i] /= tot;
      }
   }
   // the neurons themselves
   private Neuron[] m_neurons;
   // used for calculating distances squared
   private float distance2(float[] a, float[] b) {
      float tot = 0;
      for (int i = 0; i < a.length; ++i) {
         tot += (a[i] - b[i])*(a[i] - b[i]);
      }
      return tot;
   }
   // this private class is used to order neurons by distance
   private class DistCmp implements Comparator<Neuron> {
      private float[] m_dist;
      public int compare(Neuron a, Neuron b) {
         if (m_dist[a.index] < m_dist[b.index]) return -1;
         else if (m_dist[a.index] > m_dist[b.index]) return 1;
         else return 0;
      }
      public boolean equals(Object o) {return o == this;}
      public DistCmp(float[] dist) {m_dist = dist;}
   }
   // this function initializes the neurons
   private void initNeurons(Seqable coords, int nweights, long seed) {
      ISeq c = coords.seq(), s;
      int n = c.count(); // how many neurons?
      int dims = 0;
      float tmp[] = null;
      int i, j, k;
      Random rnd = (seed < 0? new Random() : new Random(seed));
      // init the basics...
      m_neurons = new Neuron[n];
      float[][] dist = new float[n][n];
      // walk through the coords...
      for (i = 0; i < n && c != null; ++i, c = c.next()) {
         s = ((Seqable)c.first()).seq();
         if (dims == 0) {
            dims = s.count();
            tmp = new float[dims];
         }
         for (j = 0; j < dims; ++j, s = s.next())
            tmp[j] = ((Number)s.first()).floatValue();
         m_neurons[i] = new Neuron(i, tmp, nweights, rnd);
         // do the distance matrix...
         for (j = 0; j < i; ++j) {
            dist[i][j] = (float)Math.sqrt(distance2(m_neurons[i].position, m_neurons[j].position));
            dist[j][i] = dist[i][j];
         }
         dist[i][i] = 0f;
      }
      // okay, now sort each neuron's neighbors by distance...
      for (i = 0; i < n; ++i) {
         m_neurons[i].neighbors = m_neurons.clone();
         Arrays.sort(m_neurons[i].neighbors, 0, m_neurons.length, new DistCmp(dist[i]));
         m_neurons[i].distances = dist[i].clone();
         Arrays.sort(m_neurons[i].distances);
      }
      // that's it!
   }

   // The memory of the network
   private float[][] m_memory;
   // how much memory we have and where we are in the memory sequence...
   private int m_memorySize;
   private int m_memPos;
   // the iteration (which is different than memPos because it starts once memory is full)
   private int m_iter;
   // the number of values in a signal
   private int m_sigSize;
   // the functions that determine the radius and the alpha learning rate for a given iteration
   private IFn m_radiusfn;
   private IFn m_alphafn;

   /** Given a Seqable of Seqables (representing a list of coordinates for the neurons), creates a
    *  KohonenAnalysis object for use with the retina simulations in brainardlab.nben.simulation.
    *  The number of weights must be specified by memory.  Memory also specifies how many signals to
    *  keep in memory for each iteration.
    */
   public KohonenAnalysis(Seqable coords, int memory, int signalsz, IFn rfn, IFn afn)
   {
      initNeurons(coords, memory, -1);
      m_memory = new float[signalsz][memory];
      m_memorySize = memory;
      m_memPos = 0;
      m_iter = 0;
      m_sigSize = signalsz;
      m_radiusfn = rfn;
      m_alphafn = afn;
   }

   public Object invoke(Object sig)
   {
      if (!(sig instanceof clojure.lang.Seqable))
         throw new IllegalArgumentException("signal must be instance of Seqable");
      int i, j, k;
      ISeq q = ((Seqable)sig).seq();
      // start by inserting this signal into the memory...
      int p = (m_memPos % m_memorySize); // where we put this signal in the memory
      for (i = 0; i < m_sigSize && q != null; ++i, q = q.next())
         m_memory[i][p] = ((Number)q.first()).floatValue();
      if (i != m_sigSize || q != null)
         throw new IllegalArgumentException("signal was the wrong size");
      ++m_memPos;
      // we want to do nothing unless we now have a full memory already
      if (m_memPos < m_memory.length) return null;
      // Okay, we have an actual iteration! we need to find the bmus for each neuron
      float r = ((Number)m_radiusfn.invoke(new Integer(m_iter))).floatValue();
      float a = ((Number)m_alphafn.invoke(new Integer(m_iter))).floatValue();
      // put the weights into the old-weights position
      for (i = 0; i < m_neurons.length; ++i)
         m_neurons[i].oldWeights = m_neurons[i].weights.clone();
      // best matching unit for each memory-length set of the signal
      float bestMatch, tmp;
      Neuron bmu, neu;
      ISeq ret = null;
      // we do this backwards to have a forward list to return
      for (i = m_sigSize-1; i >= 0; --i) {
         // find the best matching unit
         bmu = m_neurons[0];
         bestMatch = distance2(bmu.oldWeights, m_memory[i]);
         for (j = 1; j < m_neurons.length; ++j) {
            tmp = distance2(m_neurons[j].oldWeights, m_memory[i]);
            if (tmp < bestMatch) {
               bmu = m_neurons[j];
               bestMatch = tmp;
            }
         }
         // add this bmu to the list
         ret = new Cons(new Integer(bmu.index), ret);
         // adjust weights for nearby units
         for (j = 0; j < m_neurons.length && bmu.distances[j] < r; ++j) {
            neu = bmu.neighbors[j];
            for (k = 0; k < bmu.weights.length; ++k)
               neu.weights[k] += a * (m_memory[i][k] - neu.weights[k]);
         }
      }
      // that is; we've processed all the signals!
      ++m_iter;
      return ret;
   }

   public Object invoke()
   {
      // return the weights in a seq
      ISeq ret = null;
      // we do this backwards to have a forward list to return
      for (int i = m_neurons.length - 1; i >= 0; --i)
         ret = new Cons(RT.seq(m_neurons[i].weights), ret);
      return ret;
   }
}