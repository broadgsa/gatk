package org.broadinstitute.sting.utils;


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
