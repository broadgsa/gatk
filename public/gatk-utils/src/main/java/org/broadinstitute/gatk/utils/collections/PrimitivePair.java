/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.collections;


/** This class is used to group together multiple Pair classes for
 *  primitive types (thanks to generics shortcomings, these implementations
 *  are more efficient then generic ones). This class contains no methods and
 *  no fields, but only declarations of inner classes.
 */
 
public class PrimitivePair {

   /** Pair of two integers */
  public static class Int {
    // declare public, STL-style for easier and more efficient access:
    public int first; 
    public int second;

    public Int(int x, int y) { first = x; second = y; }
    public Int() { first = second = 0; }

    public void set(int x, int y) { first = x; second = y; }

    /** Java-style getter; note that we currently allow direct access to 
        the member field.
    */
    public int getFirst() { return first; }

    /** Java-style getter; note that we currently allow direct access to 
        the member field.
    */
    public int getSecond() { return second; }

       /** Increments the elements of this pair by the
        * corresponding elements of the pair <code>p</code> and returns this
        * pair (modified). This method does not allocate a new pair, but changes
        * in place the values stored in the object the method is invoked from. The
        * method is unsafe: if p is null, a runtime exception will be thrown.
        * @param p
        * @return
        */
    public PrimitivePair.Int add(PrimitivePair.Int p) {
        first += p.first;
        second += p.second;
        return this;
    }

       /** Decrements the elements of this pair by the
        * corresponding elements of the pair <code>p</code> and returns this
        * pair (modified). This method does not allocate a new pair, but changes
        * in place the values stored in the object the method is invoked from. The
        * method is unsafe: if p is null, a runtime exception will be thrown.
        * @param p
        * @return
        */
    public PrimitivePair.Int subtract(PrimitivePair.Int p) {
        first -= p.first;
        second -= p.second;
        return this;
    }

       /** Copies values from the argument <code>p</code> into the corresponding
        * elements of this pair and returns this pair (modified).
        * @param p
        * @return
        */
    public PrimitivePair.Int assignFrom(PrimitivePair.Int p ) {
        first = p.first;
        second = p.second;
        return this;
    }


  }

    public static class Long {
      // declare public, STL-style for easier and more efficient access:
      public long first;
      public long second;

      public Long(long x, long y) { first = x; second = y; }
      public Long() { first = second = 0; }

      public void set(long x, long y) { first = x; second = y; }

      /** Java-style getter; note that we currently allow direct access to
          the member field.
      */
      public long getFirst() { return first; }

      /** Java-style getter; note that we currently allow direct access to
          the member field.
      */
      public long getSecond() { return second; }

        /** Increments the elements of this pair by the
         * corresponding elements of the pair <code>p</code> and returns this
         * pair (modified). This method does not allocate a new pair, but changes
         * in place the values stored in the object the method is invoked from. The
         * method is unsafe: if p is null, a runtime exception will be thrown.
         * @param p
         * @return
         */
     public PrimitivePair.Long add(PrimitivePair.Int p) {
         first += p.first;
         second += p.second;
         return this;
     }

        /** Increments the elements of this pair by the
         * corresponding elements of the pair <code>p</code> and returns this
         * pair (modified). This method does not allocate a new pair, but changes
         * in place the values stored in the object the method is invoked from. The
         * method is unsafe: if p is null, a runtime exception will be thrown.
         * @param p
         * @return
         */
     public PrimitivePair.Long add(PrimitivePair.Long p) {
         first += p.first;
         second += p.second;
         return this;
     }

        /** Decrements the elements of this pair by the
         * corresponding elements of the pair <code>p</code> and returns this
         * pair (modified). This method does not allocate a new pair, but changes
         * in place the values stored in the object the method is invoked from. The
         * method is unsafe: if p is null, a runtime exception will be thrown.
         * @param p
         * @return
         */
     public PrimitivePair.Long subtract(PrimitivePair.Int p) {
         first -= p.first;
         second -= p.second;
         return this;
     }

        /** Decrements the elements of this pair by the
         * corresponding elements of the pair <code>p</code> and returns this
         * pair (modified). This method does not allocate a new pair, but changes
         * in place the values stored in the object the method is invoked from. The
         * method is unsafe: if p is null, a runtime exception will be thrown.
         * @param p
         * @return
         */
     public PrimitivePair.Long subtract(PrimitivePair.Long p) {
         first -= p.first;
         second -= p.second;
         return this;
     }

     /** Copies values from the argument <code>p</code> into the corresponding
       * elements of this pair and returns this pair (modified).
       * @param p
       * @return
     */
     public PrimitivePair.Long assignFrom(PrimitivePair.Long p ) {
         first = p.first;
         second = p.second;
         return this;
     }

     /** Copies values from the argument <code>p</code> into the corresponding
       * elements of this pair and returns this pair (modified).
       * @param p
       * @return
     */
     public PrimitivePair.Long assignFrom(PrimitivePair.Int p ) {
            first = p.first;
            second = p.second;
            return this;
     }

    }

}
