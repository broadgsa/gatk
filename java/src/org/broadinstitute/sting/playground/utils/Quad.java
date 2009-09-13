package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.utils.Pair;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 11, 2009
 * Time: 5:12:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class Quad<W,X,Y,Z> {
    public W first;
    public X second;
    public Y third;
    public Z fourth;

    public Quad() {
        first = null;
        second = null;
        third = null;
        fourth = null;
    }

    public Quad(W w, X x, Y y, Z z) {
        first = w;
        second = x;
        third = y;
        fourth = z;
    }

    public Quad(Pair<W,X> a, Pair<Y,Z> b) {
        first = a.getFirst();
        second = a.getSecond();
        third = b.getFirst();
        fourth = b.getSecond();
    }

    public boolean equals(Object o) {
        if(o == null) {
            return false;
        } else if (! (o instanceof Quad) ) {
            return false;
        }

        Quad other = (Quad) o;

        return ( equalToNotNull(this.first,other.first) && equalToNotNull(this.second,other.second)
                 && equalToNotNull(this.third,other.third) && equalToNotNull(this.fourth,other.fourth));
    }

    public int hashCode() {
        return getHash(first) ^ getHash(second) ^ getHash(third) ^ getHash(fourth);
    }

    public W getFirst() { return first; }
    public X getSecond() { return second; }
    public Y getThird() { return third; }
    public Z getFourth() { return fourth; }
    
    public Pair<W,X> getFirstPair() { return new Pair<W,X>(first,second); }
    public Pair<Y,Z> getSecondPair() { return new Pair<Y,Z>(third,fourth); }

    public void setFirst(Object o) { first = (W) o; }
    public void setSecond(Object o) { second = (X) o; }
    public void setThird(Object o) { third = (Y) o; }
    public void setFourth(Object o) { fourth = (Z) o; }

    private int getHash(Object o) {
        int hash = 0;
        if(o != null) {
            hash = o.hashCode();
        }
        return hash;
    }

    private boolean equalToNotNull(Object a, Object b) {
        boolean areEqual = false;
        if ( a != null && b != null ) {
            areEqual = a.equals(b);
        } else if (a == null && b == null ) {
            areEqual = true;
            // todo -- make sure we don't want to check for instanceOf here...
            // todo -- maybe this statement should be eliminated
        }

        return areEqual;
    }

}
