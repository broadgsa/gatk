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
  }
}
