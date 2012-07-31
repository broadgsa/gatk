package org.broadinstitute.sting.utils.activeregion;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 7/27/12
 */

public class ActivityProfileResult {
    public double isActiveProb;
    public ActivityProfileResultState resultState;
    public Number resultValue;

    public enum ActivityProfileResultState {
        NONE,
        HIGH_QUALITY_SOFT_CLIPS
    }

    public ActivityProfileResult( final double isActiveProb ) {
        this.isActiveProb = isActiveProb;
        this.resultState = ActivityProfileResultState.NONE;
        this.resultValue = null;
    }

    public ActivityProfileResult( final double isActiveProb, final ActivityProfileResultState resultState, final Number resultValue ) {
        this.isActiveProb = isActiveProb;
        this.resultState = resultState;
        this.resultValue = resultValue;
    }

}
