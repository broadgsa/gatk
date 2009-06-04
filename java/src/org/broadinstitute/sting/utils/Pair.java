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

    /**
     * Calculate whether this pair object is equal to another object.
     * @param o The other object (hopefully a pair).
     * @return True if the two are equal; false otherwise.
     */
    @Override
    public boolean equals( Object o ) {
        if( o == null )
            return false;
        if( !(o instanceof Pair) )
            return false;

        Pair other = (Pair)o;

        // Check to see whether one is null but not the other.
        if( this.first == null && other.first != null ) return false;
        if( this.second == null && other.second != null ) return false;

        // Check to see whether the values are equal.
        //  If the param of equals is null, it should by contract return false.
        if( this.first != null && !this.first.equals(other.first) ) return false;
        if( this.second != null && !this.second.equals(other.second) ) return false;        

        return true;
    }

    /**
     * Basic hashcode function.  Assume hashcodes of first and second are
     * randomly distributed and return the XOR of the two.
     * @return Randomly distributed hashcode of the pair.
     */
    @Override
    public int hashCode() {
        if( second == null && first == null )
            return 0;
        if( second == null )
            return first.hashCode();
        if( first == null )
            return second.hashCode();
        return first.hashCode() ^ second.hashCode();
    }

    public String toString() {
        return first+","+second;
    }
}
