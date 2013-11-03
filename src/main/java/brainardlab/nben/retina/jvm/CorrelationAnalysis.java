////////////////////////////////////////////////////////////////////////////////////////////////////
// CorrelationAnalysis.java
// Java code used for correlation matrix generation during simulation of retinas.
// by Noah C. Benson

package brainardlab.nben.retina.jvm;
import clojure.lang.*;

public final class CorrelationAnalysis
   extends clojure.lang.AFn
{
   // private data: this stores the data shown so far that is necessary for eventual calculation of
   // the correlation matrix
   private float[] m_Ex;
   private float[][] m_Exy;
   private int m_cones;
   private int m_signals;

   // the private static data contains retina keywords
   private static clojure.lang.Keyword s_labels;
   static {
      s_labels = clojure.lang.Keyword.intern("labels");
   }
   
   /** Given a retina size, creates a CorrelationAnalysis object ready to be fed signals via the 
       object's invoke(signal) method.  For the result of the calculation, call invoke(). */
   public CorrelationAnalysis(int signal_size)
   {
      m_cones = signal_size;
      m_signals = 0;
      m_Ex = new float[m_cones];
      m_Exy = new float[m_cones][m_cones];
   }

   public Object invoke(Object sig)
   {
      if (!(sig instanceof clojure.lang.Seqable))
         throw new IllegalArgumentException("signal must be instance of Seqable");
      float[] tmp = new float[m_cones];
      int i = 0, j;
      for (ISeq q = ((Seqable)sig).seq(); q != null; q = q.next(), ++i) {
         tmp[i] = ((Number)q.first()).floatValue();
         m_Ex[i] += tmp[i];
         for (j = 0; j <= i; ++j)
            m_Exy[i][j] += tmp[i] * tmp[j];
      }
      ++m_signals;
      return null;
   }

   public Object invoke()
   {
      if (m_signals < 2) throw new IllegalStateException("At least 2 signals must be shown!");
      float[][] corr = new float[m_cones][m_cones];
      float[] std = new float[m_cones];
      float[] Ex = new float[m_cones];
      int i, j;
      // covariance is Exy - ExEy; corr is this divided by sqrt(VxVy); we can do this all in one
      // pair of nested loops for optimum efficiency
      for (i = 0; i < m_cones; ++i) {
         Ex[i] = m_Ex[i] / m_signals;
         std[i] = (float)Math.sqrt((double)(m_Exy[i][i]/m_signals - Ex[i]*Ex[i]));
         for (j = 0; j < i; ++j) {
            corr[i][j] = (m_Exy[i][j]/m_signals - Ex[i]*Ex[j]) / (std[i] * std[j]);
            corr[j][i] = corr[i][j];
         }
         corr[i][i] = 1f;
      }
      // that's it!
      return new Matrix2D(corr);
   }
}