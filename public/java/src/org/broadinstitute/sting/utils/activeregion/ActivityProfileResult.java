package org.broadinstitute.sting.utils.activeregion;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 7/27/12
 */
public class ActivityProfileResult {
    private GenomeLoc loc;
    public double isActiveProb;
    public ActivityProfileResultState resultState;
    public Number resultValue;

    public enum ActivityProfileResultState {
        NONE,
        HIGH_QUALITY_SOFT_CLIPS
    }

    /**
     * Create a new ActivityProfileResult at loc with probability of being active of isActiveProb
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param isActiveProb the probability of being active (between 0 and 1)
     */
    @Requires({"loc != null", "isActiveProb >= 0.0 && isActiveProb <= 1.0"})
    public ActivityProfileResult( final GenomeLoc loc, final double isActiveProb ) {
        this(loc, isActiveProb, ActivityProfileResultState.NONE, null);
    }

    /**
     * Create a new ActivityProfileResult at loc with probability of being active of isActiveProb that maintains some
     * information about the result state and value (TODO RYAN -- what do these mean?)
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param isActiveProb the probability of being active (between 0 and 1)
     */
    @Requires({"loc != null", "isActiveProb >= 0.0 && isActiveProb <= 1.0"})
    public ActivityProfileResult( final GenomeLoc loc, final double isActiveProb, final ActivityProfileResultState resultState, final Number resultValue ) {
        // make sure the location of that activity profile is 1
        if ( loc.size() != 1 )
            throw new IllegalArgumentException("Location for an ActivityProfileResult must have to size 1 bp but saw " + loc);

        this.loc = loc;
        this.isActiveProb = isActiveProb;
        this.resultState = resultState;
        this.resultValue = resultValue;
    }

    /**
     * Get the genome loc associated with the ActivityProfileResult
     * @return the location of this result
     */
    @Ensures("result != null")
    public GenomeLoc getLoc() {
        return loc;
    }

    @Override
    public String toString() {
        return "ActivityProfileResult{" +
                "loc=" + loc +
                ", isActiveProb=" + isActiveProb +
                ", resultState=" + resultState +
                ", resultValue=" + resultValue +
                '}';
    }
}
