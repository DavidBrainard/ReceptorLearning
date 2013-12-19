////////////////////////////////////////////////////////////////////////////////////////////////////
// HSImage.java
// This class is similar to Matrix2D except that it is designed to contain efficient representations
// of the (rather large) hyperspectral images used in the retina simulations.
// by Noah C. Benson

package brainardlab.nben.retina.jvm;
import clojure.lang.*;
import java.util.Random;

public final class HSImage
   extends AFn implements Seqable
{
   // private classes that handle kinds of data given to the HSImage object (part of the point is 
   // that we never reallocate/copy the data in order to preserve memory).
   private interface HSData {
      ISeq get(int i, int j);
      ISeq all();
      HSData duplicate(HSImage newThis);
   }
   private final HSData new_Row_Col_Depth(float[][][] data) {return new Row_Col_Depth(data);}
   private final class Row_Col_Depth implements HSData {
      float[][][] m_data;
      public Row_Col_Depth(float[][][] data) {
         m_data = data;
      }
      public HSData duplicate(HSImage newThis) {
         return newThis.new_Row_Col_Depth(m_data);
      }
      private class Iter extends ASeq {
         private int m_i;
         private float[] m_data;
         public Iter(float[] data, int i) {m_data = data; m_i = i;}
         public Object first() {return new Float(m_data[m_i]);}
         public ISeq next() {return (m_i+1 == m_data.length? null : new Iter(m_data, m_i+1));}
         public Obj withMeta(IPersistentMap m) {return this;}
      }
      private final class ColSeq extends ASeq {
         private float[][] m_dat;
         private int m_col;
         public ColSeq(float[][] dat, int col) {
            m_dat = dat;
            m_col = col;
         }
         public Object first() {return new Iter(m_dat[m_col0 + m_col], 0);}
         public ISeq next() {return (m_col + 1 == m_cols? null : new ColSeq(m_dat, m_col+1));}
         public Obj withMeta(IPersistentMap m) {return this;}
      }
      private final class RowSeq extends ASeq {
         private int m_row;
         public RowSeq(int row) {
            m_row = row;
         }
         public Object first() {return new ColSeq(m_data[m_row0 + m_row], 0);}
         public ISeq next() {return (m_row + 1 == m_rows? null : new RowSeq(m_row+1));}
         public Obj withMeta(IPersistentMap m) {return this;}
      }
      public ISeq get(int i, int j) {return new Iter(m_data[i][j], 0);}
      public ISeq all() {return new RowSeq(0);}
   }
   private final HSData new_Row_ColDepth(float[][] data, int col) {
      return new Row_ColDepth(data, col);
   }
   private final class Row_ColDepth implements HSData {
      float[][] m_data;
      int m_cols;
      int m_wlens;
      public Row_ColDepth(float[][] data, int cols) {
         m_data = data;
         m_cols = cols;
         m_wlens = data[0].length / cols;
      }
      public HSData duplicate(HSImage newThis) {
         return newThis.new_Row_ColDepth(m_data, m_cols);
      }
      private final class Iter extends ASeq {
         private float[] m_data;
         private int m_col;
         private int m_i;
         public Iter(float[] data, int col, int i) {
            m_data = data;
            m_col = col;
            m_i = i;
         }
         public Object first() {return new Float(m_data[m_col + m_i*m_cols]);}
         public ISeq next() {return (m_i+1 == m_wlens? null : new Iter(m_data, m_col, m_i+1));}
         public Obj withMeta(IPersistentMap m) {return this;}
      }
      private final class ColSeq extends ASeq {
         private float[] m_dat;
         private int m_col;
         public ColSeq(float[] dat, int col) {
            m_dat = dat;
            m_col = col;
         }
         public Object first() {return new Iter(m_dat, m_col0 + m_col, 0);}
         public ISeq next() {
            return (m_col + 1 == m_cols? null : new ColSeq(m_dat, m_col+1));
         }
         public Obj withMeta(IPersistentMap m) {return this;}
      }
      private final class RowSeq extends ASeq {
         private int m_row;
         public RowSeq(int row) {
            m_row = row;
         }
         public Object first() {return new ColSeq(m_data[m_row0 + m_row], 0);}
         public ISeq next() {return (m_row + 1 == m_rows? null : new RowSeq(m_row+1));}
         public Obj withMeta(IPersistentMap m) {return this;}
      }
      public ISeq get(int i, int j) {return new Iter(m_data[i], j, 0);}
      public ISeq all() {return new RowSeq(0);}
   }

   private HSData m_data;
   private IPersistentVector m_dims;
   int m_row0;
   int m_col0;
   int m_rows;
   int m_cols;
   int m_wlens;

   /** Creates an image in which the ith row, jth column, and kth wavelength are stored at position
    *  data[i][j + cols*k].  This is the case when reading matlab hyperspectral images from the 
    *  Foster2004 or Chakrabarti2011 databases, but is probably generally true for 3d matrices read
    *  from a matlab .mat file.
    *
    *  @param cols The number of columns in the image.
    */
   public HSImage(float[][] data, int cols)
   {
      m_data = new Row_ColDepth(data, cols);
      m_dims = PersistentVector.create(data.length, cols, data[0].length / cols);
      m_row0 = 0;
      m_col0 = 0;
      m_rows = data.length;
      m_cols = cols;
      m_wlens = data[0].length / cols;
   }

   /** Creates an image in which the ith row, jth column, and kth wavelength are stored at position
    *  data[i][j + cols*k].  This is the case when reading matlab hyperspectral images from the 
    *  Foster2004 or Chakrabarti2011 databases, but is probably generally true for 3d matrices read
    *  from a matlab .mat file.  The final argument cal must be a Sequable with length equal to the
    *  number of wavelengths in the data, and will guarantee that each wavelength is divided by the
    *  appropriate calibration value.
    *
    *  @param cols The number of columns in the image.
    */
   public HSImage(float[][] data, int cols, IPersistentCollection cal)
   {
      m_row0 = 0;
      m_col0 = 0;
      m_rows = data.length;
      m_cols = cols;
      m_wlens = data[0].length / cols;
      // implement the calibration...
      int i, j, k;
      ISeq q;
      float[] c = new float[m_wlens];
      for (k = 0, q = cal.seq(); q != null; ++k, q = q.next())
         c[k] = ((Number)q.first()).floatValue();
      float[] row;
      for (i = 0; i < m_rows; ++i) {
         row = data[i];
         for (j = 0; j < m_cols; ++j) {
            for (k = 0; k < m_wlens; ++k)
               row[j + k*m_cols] /= c[k];
         }
      }
      m_data = new Row_ColDepth(data, cols);
      m_dims = PersistentVector.create(data.length, cols, m_wlens);
   }

   /** Creates an image out of a clojure Seqable; the seqable must yield an n x m x q non-ragged
    *  seq which will be interpreted as the image.  This constructor duplicates the data into a more
    *  efficient representation, so ideally the seqable should be lost for garbage collection after
    *  the HSImage is created.
    */
   public HSImage(Seqable img)
   {
      // first things first, count the first dim
      int rows = 0;
      int cols = 0;
      int wlens = 0;
      ISeq r, c, w;
      Seqable sqr, sqc;
      float[][][] data;
      for (r = img.seq(); r != null; r = r.next())
         ++rows;
      // okay, now go through each row...
      data = new float[rows][][];
      int rnum, cnum, wnum;
      for (rnum = 0, r = img.seq(); r != null; r = r.next(), ++rnum) {
         sqr = (Seqable)r.first();
         // count columns...
         if (cols == 0) {
            for (c = sqr.seq(); c != null; c = c.next())
               ++cols;
         }
         // allocate this column and fill it in...
         data[rnum] = new float[cols][];
         for (cnum = 0, c = sqr.seq(); c != null && cnum < cols; c = c.next(), ++cnum) {
            sqc = (Seqable)c.first();
            // count wavelengths...
            if (wlens == 0) {
               for (w = sqc.seq(); w != null; w = w.next())
                  ++wlens;
            }
            // allocate this pixel and fill it in...
            data[rnum][cnum] = new float[wlens];
            for (wnum = 0, w = sqc.seq(); w != null && wnum < wlens; w = w.next(), ++wnum)
               data[rnum][cnum][wnum] = ((Number)w.first()).floatValue();
            if (wnum != wlens) 
               throw new IllegalArgumentException("Seqable given to HSImage was ragged at [" 
                                                  + rnum + " " +  cnum + "]");
         }
         if (cnum != cols)
            throw new IllegalArgumentException("Seqable given to HSImage was ragged at row "+rnum);
      }
      // done!
      m_data = new Row_Col_Depth(data);
      m_dims = PersistentVector.create(rows, cols, wlens);
      m_row0 = 0;
      m_col0 = 0;
      m_rows = rows;
      m_cols = cols;
      m_wlens = wlens;
   }

   /** Creates an HSImage object out of the given data array. */
   public HSImage(float[][][] data)
   {
      m_rows = data.length;
      m_cols = data[0].length;
      m_wlens = data[0][0].length;
      m_data = new Row_Col_Depth(data);
      m_dims = PersistentVector.create(data.length, data[0].length, data[0][0].length);
   }

   private HSImage(HSData data, int wlens, int row0, int rows, int col0, int cols)
   {
      m_dims = PersistentVector.create(rows, cols, wlens);
      m_row0 = row0;
      m_col0 = col0;
      m_rows = rows;
      m_cols = cols;
      m_wlens = wlens;
      m_data = data.duplicate(this);
   }
   
   /** Yields the dimensions of the hyperspectral image as a persistent vector of 
    *  [rows columns wavelengths].
    */
   public Object invoke()
   {
      return m_dims;
   }
   /** Yields an ISeq of the wavelengths at (row, col) in the hyperspectral image. */
   public Object invoke(Object row, Object col)
   {
      int i = ((Number)row).intValue();
      int j = ((Number)col).intValue();
      if (i >= m_rows || j >= m_cols) 
         throw new IndexOutOfBoundsException("HSImage: requested [" + i + " " + j + "] of image" +
                                             " with " + m_rows + " rows and " + m_cols + " cols");
      return m_data.get(m_row0 + i, m_col0 + j);
   }
   /** Yields a sub-image of the hyperspectral image.
    *
    *  @param row0 the first row of the new image (inclusive)
    *  @param rows the number of rows in the new image
    *  @param col0 the first column of the new image (inclusive)
    *  @param cols the number of columns in the new image
    */
   public Object invoke(Object row0, Object rows, Object col0, Object cols)
   {
      return new HSImage(m_data,
                         m_wlens,
                         ((Number)row0).intValue(),
                         ((Number)rows).intValue(),
                         ((Number)col0).intValue(),
                         ((Number)cols).intValue());
   }
   /** Yields a seq of the 3d array in the hyperspectral image */
   public ISeq seq()
   {
      return m_data.all();
   }

   /** Yields a clojure IPersistentMap of [pixel-distance wavelength-difference] mapped to 
    *  a correlation.
    *
    *  @param maxDist The maximum distance in terms of rows or columns to collect.
    *  @param numSamples The maximum of pairs of pixels to collect.
    */
   public IPersistentMap collectStatistics(int maxDist, int numSamples)
   {
      if (maxDist < 1 || numSamples < 2)
         throw new IllegalArgumentException("maxDist and numSamples must be > 1 and > 2");
      ITransientMap result = PersistentHashMap.EMPTY.asTransient();
      IPersistentMap pres;
      ISeq qhold1, qhold2, q1, q2;
      int r, c, mnr, mxr, mnc, mxc, rr, cc, w1, w2;
      Random rand = new Random();
      IPersistentVector v;
      float f1, f2;
      Object tmp;
      double[] dat;
      MapEntry me;
      double d;
      for (int sampleNum = 0; sampleNum < numSamples; ++sampleNum) {
         r = rand.nextInt(m_rows) + m_row0;
         c = rand.nextInt(m_cols) + m_col0;
         mnr = (r - maxDist < m_row0 ? m_row0 : r - maxDist);
         mxr = (r + maxDist > m_rows? m_rows : r + maxDist);
         mnc = (c - maxDist < m_col0? m_col0 : c - maxDist);
         mxc = (c + maxDist > m_cols? m_cols : c + maxDist);
         qhold1 = m_data.get(r, c);
         for (rr = mnr; rr < mxr; ++rr) {
            for (cc = mnc; cc < mxc; ++cc) {
               if (rr == r && cc == c) continue;
               qhold2 = m_data.get(rr, cc);
               d = (r - rr)*(r - rr) + (c - cc)*(c - cc);
               v = PersistentVector.create(new Double(d), new Integer(0));
               w1 = 0;
               q1 = qhold1;
               while (w1 < m_wlens) {
                  f1 = ((Number)q1.first()).floatValue();
                  w2 = 0;
                  q2 = qhold2;
                  while (w2 < m_wlens) {
                     f2 = ((Number)q2.first()).floatValue();
                     v = v.assocN(1, new Integer(w1 > w2? w1 - w2 : w2 - w1));
                     tmp = result.valAt(v);
                     if (tmp == null) {
                        dat = new double[6];
                        result.assoc(v, dat);
                     } else 
                        dat = (double[])tmp;
                     dat[0] = dat[0] + 1.0;
                     dat[1] = dat[1] + f1;
                     dat[2] = dat[2] + f2;
                     dat[3] = dat[3] + f1*f1;
                     dat[4] = dat[4] + f2*f2;
                     dat[5] = dat[5] + f1*f2;
                     ++w2;
                     q2 = q2.next();
                  }
                  ++w1;
                  q1 = q1.next();
               }
            }
         }
      }
      double eX, eY, eXX, eYY, eXY, varX, varY, cov, corr;
      pres = result.persistent();
      IPersistentMap finalized = PersistentHashMap.EMPTY;
      for (q1 = pres.seq(); q1 != null; q1 = q1.next()) {
         me = (MapEntry)q1.first();
         v = (IPersistentVector)me.key();
         dat = (double[])me.val();
         eX = dat[1] / dat[0];
         eY = dat[2] / dat[0];
         eXX = dat[3] / dat[0];
         eYY = dat[4] / dat[0];
         eXY = dat[5] / dat[0];
         varX = eXX - eX*eX;
         varY = eYY - eY*eY;
         cov = eXY - eX*eY;
         corr = cov / Math.sqrt(varX * varY);
         v.assocN(0, Math.sqrt(((Double)v.nth(0)).doubleValue()));
         finalized = finalized.assoc(v, new Double(corr));
      }
      return finalized;
   }
}
