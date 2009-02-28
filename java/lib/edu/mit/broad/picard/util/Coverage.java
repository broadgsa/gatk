package edu.mit.broad.picard.util;

/**
 * A simple class that is used to store the coverage information about an interval.
 *
 * @author Tim Fennell 
 */
public class Coverage {
    private Interval interval;
    private short[] depths;

    /** Constructs a new coverage object for the provided mapping with the desired padding either side. */
    public Coverage(Interval i, int padding) {
        this.interval = i;
        this.depths = new short[interval.length() + 2*padding];
    }

    /** Adds a single point of depth at the desired offset into the coverage array. */
    public void addBase(int offset) {
        if (offset >= 0 && offset < this.depths.length) {
            this.depths[offset] += 1;
        }
    }

    /** Returns true if any base in the range has coverage of > 1 */
    public boolean hasCoverage() {
        for (short s : depths) {
            if (s > 1) return true;
        }

        return false;
    }

    /** Gets the coverage depths as an array of shorts. */
    public short[] getDepths() { return this.depths; }
}
