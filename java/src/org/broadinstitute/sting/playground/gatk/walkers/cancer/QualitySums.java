package org.broadinstitute.sting.playground.gatk.walkers.cancer;

/**
 * Created by IntelliJ IDEA.
 * User: kcibul
 * Date: Jun 27, 2009
 * Time: 3:43:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class QualitySums {
    private int a = 0;
    private int c = 0;
    private int g = 0;
    private int t = 0;
    private int aCounts = 0;
    private int cCounts = 0;
    private int gCounts = 0;
    private int tCounts = 0;

    public int getQualitySum(final char base) {
        if (base == 'a' || base == 'A') { return a; }
        if (base == 'c' || base == 'C') { return c; }
        if (base == 'g' || base == 'G') { return g; }
        if (base == 't' || base == 'T') { return t; }
        throw new RuntimeException("Unknown base: " + base);
    }

    public int getCounts(final char base) {
        if (base == 'a' || base == 'A') { return aCounts; }
        if (base == 'c' || base == 'C') { return cCounts; }
        if (base == 'g' || base == 'G') { return gCounts; }
        if (base == 't' || base == 'T') { return tCounts; }
        throw new RuntimeException("Unknown base: " + base);
    }

    public void incrementSum(final char base, final byte qual) {
        if (base == 'a' || base == 'A')      { a += qual; aCounts++;}
        else if (base == 'c' || base == 'C') { c += qual; cCounts++;}
        else if (base == 'g' || base == 'G') { g += qual; gCounts++; }
        else if (base == 't' || base == 'T') { t += qual; tCounts++; }
        else throw new RuntimeException("Unknown base: " + base);
    }

    public int getOtherQualities(final char base) {
        int total = a + c + g + t;
        if (base == 'a' || base == 'A') { return total-a; }
        else if (base == 'c' || base == 'C') { return total-c; }
        else if (base == 'g' || base == 'G') { return total-g; }
        else if (base == 't' || base == 'T') { return total-t; }
        else throw new RuntimeException("Unknown base: " + base);
    }

    public void reset() {
        a = 0; c = 0; g = 0; t = 0;
        aCounts = 0; cCounts = 0; gCounts = 0; tCounts = 0;
    }
}

