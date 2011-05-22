/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
package org.broadinstitute.sting.gatk.walkers.phasing;

/* PreciseNonNegativeDouble permits arithmetic operations on NON-NEGATIVE double values
   with precision (prevents underflow by representing in log10 space).
 */
public class PreciseNonNegativeDouble implements Comparable<PreciseNonNegativeDouble> {
    private static final double EQUALS_THRESH = 1e-6;
    private static final double INFINITY = Double.POSITIVE_INFINITY;

    private double log10Value;

    public PreciseNonNegativeDouble(double d) {
        this(d, false);
    }

    public PreciseNonNegativeDouble(double d, boolean isLog10) {
        if (isLog10) {
            this.log10Value = d;
        }
        else {
            if (d < 0)
                throw new IllegalArgumentException("non-log PreciseNonNegativeDouble argument must be non-negative");
            this.log10Value = Math.log10(d);
        }
    }

    public PreciseNonNegativeDouble(PreciseNonNegativeDouble pd) {
        this.log10Value = pd.log10Value;
    }

    public double getValue() {
        return Math.pow(10, log10Value);
    }

    public double getLog10Value() {
        return log10Value;
    }

    public PreciseNonNegativeDouble setEqual(PreciseNonNegativeDouble other) {
        log10Value = other.log10Value;
        return this;
    }

    public PreciseNonNegativeDouble plus(PreciseNonNegativeDouble other) {
        return new PreciseNonNegativeDouble(this).plusEqual(other);
    }

    public PreciseNonNegativeDouble times(PreciseNonNegativeDouble other) {
        return new PreciseNonNegativeDouble(this).timesEqual(other);
    }

    public PreciseNonNegativeDouble div(PreciseNonNegativeDouble other) {
        return new PreciseNonNegativeDouble(this).divEqual(other);
    }

    public PreciseNonNegativeDouble absDiff(PreciseNonNegativeDouble other) {
        return new PreciseNonNegativeDouble(absSubLog(this.log10Value, other.log10Value), true);
    }

    public int compareTo(PreciseNonNegativeDouble other) {
        // Since log is monotonic: e^a R e^b <=> a R b, where R is one of: >, <, ==
        double logValDiff = this.log10Value - other.log10Value;
        if (Math.abs(logValDiff) <= EQUALS_THRESH)
            return 0; // this.equals(other)

        return new Double(Math.signum(logValDiff)).intValue();
    }

    public boolean equals(PreciseNonNegativeDouble other) {
        return (this.compareTo(other) == 0);
    }

    public boolean gt(PreciseNonNegativeDouble other) {
        return (this.compareTo(other) > 0);
    }

    public boolean lt(PreciseNonNegativeDouble other) {
        return (this.compareTo(other) < 0);
    }

    public PreciseNonNegativeDouble plusEqual(PreciseNonNegativeDouble other) {
        log10Value = addInLogSpace(log10Value, other.log10Value);
        return this;
    }

    public PreciseNonNegativeDouble timesEqual(PreciseNonNegativeDouble other) {
        log10Value += other.log10Value;
        return this;
    }

    public PreciseNonNegativeDouble divEqual(PreciseNonNegativeDouble other) {
        log10Value -= other.log10Value;
        return this;
    }

    // If x = log(a), y = log(b), returns log(a+b)
    private static double addInLogSpace(double x, double y) {
        if (x == INFINITY || y == INFINITY) return INFINITY; // log(e^INFINITY + e^y) = INFINITY

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
        return maxVal + Math.log10(1.0 + Math.pow(10, negDiff));
    }

    // If x = log(a), y = log(b), returns log |a-b|
    private double absSubLog(double x, double y) {
      if (x == -INFINITY && y == -INFINITY) {
        // log |e^-INFINITY - e^-INFINITY| = log |0-0| = log(0) = -INFINITY
        return -INFINITY;
      }
      else if (x >= y) // x + log(1-e^(y-x)) = log(a) + log(1-e^(log(b)-log(a))) = log(a) + log(1-b/a) = a - b = |a-b|, since x >= y
        return x + Math.log10(1 - Math.pow(10, y-x));
      else
        return y + Math.log10(1 - Math.pow(10, x-y));
    }    

    public String toString() {
        return new StringBuilder().append(getValue()).toString();
    }
}