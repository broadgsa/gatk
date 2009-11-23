package org.broadinstitute.sting.gatk.walkers.recalibration;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 16, 2009
 */
public class Dinuc implements Comparable<Dinuc>{
    private byte first;
    private byte second;

    public Dinuc() {
        first = 0;
        second = 0;
    }

    public Dinuc(final byte _first, final byte _second) {
        first = _first;
        second = _second;
    }

    public final void setValues(final byte _first, final byte _second) {
        first = _first;
        second = _second;
    }

    public int compareTo(final Dinuc that) {
        if( this.first > that.first ) { return 1; }
        else if( this.first < that.first ) { return -1; }
        else { //this.first equals that.first
            if( this.second > that.second ) { return 1; }
            else if( this.second < that.second ) { return -1; }
            else { return 0; }
        }

    }

    public static int hashBytes(final byte byte1, final byte byte2) {
        return byte1 << 16 + byte2 << 4;
    }

    public String toString() { // This method call is how the Dinuc will get written out to the table recalibration file
        byte[] byteArray = {first,second};
        return new String(byteArray);
    }
}
