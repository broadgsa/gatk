/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.activeregion;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.gatk.utils.GenomeLoc;

/**
 * The state of an active region walker's isActive call at a specific locus in the genome
 *
 * User: rpoplin
 * Date: 7/27/12
 */
public class ActivityProfileState {
    final private GenomeLoc loc;
    public double isActiveProb;
    public Type resultState;
    public Number resultValue;

    public enum Type {
        NONE,
        HIGH_QUALITY_SOFT_CLIPS
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of isActiveProb
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param isActiveProb the probability of being active (between 0 and 1)
     */
    @Requires({"loc != null", "isActiveProb >= 0.0 && isActiveProb <= 1.0"})
    public ActivityProfileState(final GenomeLoc loc, final double isActiveProb) {
        this(loc, isActiveProb, Type.NONE, null);
    }

    /**
     * Create a new ActivityProfileState at loc with probability of being active of isActiveProb that maintains some
     * information about the result state and value
     *
     * The only state value in use is HIGH_QUALITY_SOFT_CLIPS, and here the value is interpreted as the number
     * of bp affected by the soft clips.
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param isActiveProb the probability of being active (between 0 and 1)
     */
    @Requires({"loc != null", "isActiveProb >= 0.0 && isActiveProb <= 1.0"})
    public ActivityProfileState(final GenomeLoc loc, final double isActiveProb, final Type resultState, final Number resultValue) {
        // make sure the location of that activity profile is 1
        if ( loc.size() != 1 )
            throw new IllegalArgumentException("Location for an ActivityProfileState must have to size 1 bp but saw " + loc);
        if ( resultValue != null && resultValue.doubleValue() < 0 )
            throw new IllegalArgumentException("Result value isn't null and its < 0, which is illegal: " + resultValue);

        this.loc = loc;
        this.isActiveProb = isActiveProb;
        this.resultState = resultState;
        this.resultValue = resultValue;
    }

    /**
     * The offset of state w.r.t. our current region's start location
     * @param regionStartLoc the start of the region, as a genome loc
     * @return the position of this profile relative to the start of this region
     */
    public int getOffset(final GenomeLoc regionStartLoc) {
        return getLoc().getStart() - regionStartLoc.getStart();
    }


    /**
     * Get the genome loc associated with the ActivityProfileState
     * @return the location of this result
     */
    @Ensures("result != null")
    public GenomeLoc getLoc() {
        return loc;
    }

    @Override
    public String toString() {
        return "ActivityProfileState{" +
                "loc=" + loc +
                ", isActiveProb=" + isActiveProb +
                ", resultState=" + resultState +
                ", resultValue=" + resultValue +
                '}';
    }
}
