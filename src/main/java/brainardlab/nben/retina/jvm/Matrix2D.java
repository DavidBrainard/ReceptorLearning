////////////////////////////////////////////////////////////////////////////////////////////////////
// Matrix2D.java
// A simple 2D Matrix implementation that is space-efficient and compatible with clojure.
// by Noah C. Benson

package brainardlab.nben.retina.jvm;
import clojure.lang.*;

/** The Matrix2D class is intended to represent a large set of 2D float array data as efficiently as
 *  possible while still being friendly to clojure.  It does this by representing the data as a 
 *  clojure function (IFn) which, when invoked with no arguments returns the matrix's [rows cols] as
 *  a clojure vector, and when invoked with two arguments returns a Float object of the value at
 *  the given (row, col).
 */
public final class Matrix2D
   extends AFn
{
   // the matrix itself
   private float m_data[][];
   // the dimensions
   private IPersistentVector m_dims;

   /** Given a rectangular array of data, represents it as a clojure-friendly matrix */
   public Matrix2D(float[][] data)
   {
      if (data == null)
         throw new IllegalArgumentException("Matrix2D cannot contain a null array");
      int rows = data.length;
      if (rows == 0) {
         m_dims = PersistentVector.create(0, 0);
      } else {
         int cols = data[0].length;
         for (int i = 1; i < rows; ++i) {
            if (data[i].length != cols)
               throw new IllegalArgumentException("Matrix2D given a ragged array");
         }
         m_dims = PersistentVector.create(rows, cols);
      }
      m_data = data;
   }
   /** Given a vector and the number of columns, represents the data as a clojure-firiendly matrix
    *  in which the first <cols> floats are the first row, the next <cols> floats are the next row,
    *  etc.
    */
   public Matrix2D(float[] data, int cols)
   {
      if (data == null)
         throw new IllegalArgumentException("Matrix2D cannot contain a null array");
      if (data.length % cols != 0)
         throw new IllegalArgumentException("Matrix2D cannot contain a ragged array");
      if (data.length == 0) {
         m_dims = PersistentVector.create(0, 0);
         m_data = new float[0][0];
      } else {
         m_data = new float[data.length / cols][];
         for (int i = 0; i < m_data.length; ++i)
            m_data[i] = java.util.Arrays.copyOfRange(data, i*cols, i*cols + cols);
         m_dims = PersistentVector.create(data.length / cols, cols);
      }
   }

   /** Yields the dimensions of the Matrix as a persistent vector */
   public Object invoke()
   {
      return m_dims;
   }

   /** Yields the entry in the ith row and jth column of the matrix */
   public Object invoke(Object i, Object j)
   {
      return new Float(m_data[((Number)i).intValue()][((Number)j).intValue()]);
   }
}
