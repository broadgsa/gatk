package org.broadinstitute.sting.utils;

/**
 * Created by IntelliJ IDEA.
 * User: fromer
 * Date: Aug 18, 2010
 * Time: 4:55:23 PM
 * To change this template use File | Settings | File Templates.
 */

/* PreciseNonNegativeDouble permits arithmetic operations on NON-NEGATIVE double values
   with precision (prevents underflow by representing in log space).
 */
public class PreciseNonNegativeDouble implements Comparable<PreciseNonNegativeDouble> {
    private static double EQUALS_THRESH = 1e-6;
    private static double INFINITY = Double.POSITIVE_INFINITY;

    private double logValue;

    public PreciseNonNegativeDouble(double d) {
        this(d, false);
    }

    public PreciseNonNegativeDouble(double d, boolean isLog) {
        if (isLog) {
            this.logValue = d;
        }
        else {
            if (d < 0)
                throw new IllegalArgumentException("non-log PreciseNonNegativeDouble argument must be non-negative");
            this.logValue = Math.log(d);
        }
    }

    public PreciseNonNegativeDouble(org.broadinstitute.sting.utils.PreciseNonNegativeDouble pd) {
        this.logValue = pd.logValue;
    }

    public double getValue() {
        return Math.exp(logValue);
    }

    public double getLogValue() {
        return logValue;
    }

    public org.broadinstitute.sting.utils.PreciseNonNegativeDouble setEqual(org.broadinstitute.sting.utils.PreciseNonNegativeDouble other) {
        logValue = other.logValue;
        return this;
    }

    public org.broadinstitute.sting.utils.PreciseNonNegativeDouble plus(org.broadinstitute.sting.utils.PreciseNonNegativeDouble other) {
        return new org.broadinstitute.sting.utils.PreciseNonNegativeDouble(this).plusEqual(other);
    }

    public org.broadinstitute.sting.utils.PreciseNonNegativeDouble times(org.broadinstitute.sting.utils.PreciseNonNegativeDouble other) {
        return new org.broadinstitute.sting.utils.PreciseNonNegativeDouble(this).timesEqual(other);
    }

    public org.broadinstitute.sting.utils.PreciseNonNegativeDouble div(org.broadinstitute.sting.utils.PreciseNonNegativeDouble other) {
        return new org.broadinstitute.sting.utils.PreciseNonNegativeDouble(this).divEqual(other);
    }

    public int compareTo(org.broadinstitute.sting.utils.PreciseNonNegativeDouble other) {
        // Since log is monotonic: e^a R e^b <=> a R b, where R is one of: >, <, ==
        double logValDiff = this.logValue - other.logValue;
        if (Math.abs(logValDiff) <= EQUALS_THRESH)
            return 0; // this.equals(other)

        return new Double(Math.signum(logValDiff)).intValue();
    }

    public boolean equals(org.broadinstitute.sting.utils.PreciseNonNegativeDouble other) {
        return (this.compareTo(other) == 0);
    }

    public boolean gt(org.broadinstitute.sting.utils.PreciseNonNegativeDouble other) {
        return (this.compareTo(other) > 0);
    }

    public boolean lt(org.broadinstitute.sting.utils.PreciseNonNegativeDouble other) {
        return (this.compareTo(other) < 0);
    }

    public org.broadinstitute.sting.utils.PreciseNonNegativeDouble plusEqual(org.broadinstitute.sting.utils.PreciseNonNegativeDouble other) {
        logValue = addInLogSpace(logValue, other.logValue);
        return this;
    }

    // If x = log(a), y = log(b), returns log(a+b)
    public static double addInLogSpace(double x, double y) {
        if (x == INFINITY || y == INFINITY) return INFINITY; //log( e^INFINITY + e^y ) = INFINITY

        if (x == -INFINITY) return y;
        if (y == -INFINITY) return x;

        double maxVal, negDiff;
        if (x > y) {
            maxVal = x;
            negDiff = y - x;
        }
        else { // x <= y
            maxVal = y;
            negDiff = x - y;
        }

        // x + log(1+e^(y-x)) = log(a) + log(1+e^(log(b)-log(a))) = log(a) + log(1+b/a) = log(a+b)
        return maxVal + Math.log(1.0 + Math.exp(negDiff));
    }

    public org.broadinstitute.sting.utils.PreciseNonNegativeDouble timesEqual(org.broadinstitute.sting.utils.PreciseNonNegativeDouble other) {
        logValue += other.logValue;
        return this;
    }

    public org.broadinstitute.sting.utils.PreciseNonNegativeDouble divEqual(org.broadinstitute.sting.utils.PreciseNonNegativeDouble other) {
        logValue -= other.logValue;
        return this;
    }

    public String toString() {
        return new StringBuilder().append(getValue()).toString();
    }
}
