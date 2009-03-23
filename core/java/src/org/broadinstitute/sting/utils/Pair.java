package org.broadinstitute.sting.utils;


public class Pair<X,Y> {
    // declare public, STL-style for easier and more efficient access:
    public X first; 
    public Y second;

    public Pair(X x, Y y) { first = x; second = y; }

    public void set(X x, Y y) { first = x; second = y; }

    /** Java-style getter; note that we currently allow direct access to 
        the member field.
    */
    public X getFirst() { return first; }

    /** Java-style getter; note that we currently allow direct access to 
        the member field.
    */
    public Y getSecond() { return second; }
}
